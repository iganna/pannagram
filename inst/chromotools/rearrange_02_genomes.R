# Wind a match between chromosomes and get the blocks of interest

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(pannagram)
  library(crayon)
})

source(system.file("chromotools/rearrange_func.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("--path.aln"),       type = "character", default = NULL, help = "Path to the output directory with alignments"),
  make_option(c("--path.chr"),       type = "character", default = NULL, help = "Path to the output directory with chromosomes"),
  make_option(c("--path.processed"), type = "character", default = NULL, help = "Path to the output directory with processed chromosomes"),
  
  make_option(c("--ref"),            type = "character", default = NULL, help = "Name of the reference genome"),
  
  make_option(c("--cores"),          type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option(c("--path.log"),       type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"),      type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)


# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- Values of parameters ----

# Number of cores
num.cores <- opt$cores

path.aln       <- ifelse(!is.null(opt$path.aln), opt$path.aln, stop('Folder with Alignments is not specified'))
path.chr       <- ifelse(!is.null(opt$path.chr), opt$path.chr, stop('Folder with Chromosomes is not specified'))
path.processed <- ifelse(!is.null(opt$path.processed), opt$path.processed, stop('Output folder with processed genomes is not specified'))
ref            <- ifelse(!is.null(opt$ref), opt$ref, stop('Reference genome is not specified'))

checkDir(path.aln)
checkDir(path.chr)
checkDir(path.processed)

pokaz('Output folder with processed genomes:', path.processed)

# Create folders for the alignment results
if(!dir.exists(path.processed)) dir.create(path.processed)

# Accessions
pokaz('Folder with alignments:', path.aln)
files.aln <- list.files(path.aln, pattern = ".*_maj\\.rds$", full.names = F)
files.aln = sub("_maj.rds", "", files.aln)
accessions = c()
for(f in files.aln){
  s = strsplit(f, '_')[[1]]
  s = s[-(length(s)-(0:1))]
  s = paste0(s, collapse = '_')
  accessions = c(accessions, s)
  accessions = unique(accessions)
}

pokaz('Accessions:', accessions)

# Accessions without correspondence (skipped in rearrange_01) are not rearranged
files.corresp = paste0(path.processed, 'corresp_', accessions, '_to_', ref, '.rds')
if(any(!file.exists(files.corresp))){
  pokazAttention('No correspondence file, accessions are not rearranged:', accessions[!file.exists(files.corresp)])
  accessions = accessions[file.exists(files.corresp)]
  if(length(accessions) == 0) stop('No accessions with correspondence to the reference')
}

# ***********************************************************************
# ---- Test ----

# ref = 'GCA_002079055.1'
# path.aln = "~/Library/CloudStorage/OneDrive-Personal/iglab/projects/pannagram_meta/yeast_wild/alignment/alignments_GCA_002079055.1/"
# path.chr = "~/Library/CloudStorage/OneDrive-Personal/iglab/projects/pannagram_meta/yeast_wild/alignment/chromosomes/"
# acc = "GCA_002079175.1"
# 
# 
# ref = 'MN47'
# acc = '1741'
# path.aln = "~/Library/CloudStorage/OneDrive-Personal/iglab/projects/pannagram_meta/arabidopsis/alignment/"
# path.chr = "~/Library/CloudStorage/OneDrive-Personal/iglab/projects/pannagram_meta/arabidopsis/chromosomes/"


# ***********************************************************************
# ---- Variables ----

file.ref.len = paste0(path.chr, ref, '_chr_len.txt', collapse = '')
if(!file.exists(file.ref.len)) stop('File', file.ref.len, 'does not exist.')
ref.len = read.table(file.ref.len, header = 1)
min.overlap = 0.01
min.overlap.fragment = 0.01


# ***********************************************************************
# ---- Define the chromosomes, which have the correspondence in all accessions ----
ref.chromosomes = c()
for(acc in accessions){
  pokaz('Accession', acc)
  
  file.acc.len = paste0(path.chr, acc, '_chr_len.txt', collapse = '')
  acc.len = read.table(file.acc.len, header = 1)
  
  corresp.acc2ref = readRDS(paste0(path.processed, 'corresp_',acc,'_to_',ref, '.rds'))
  if(acc == accessions[1]){
    ref.chromosomes = unique(corresp.acc2ref$i.ref)
  } else {
    ref.chromosomes = intersect(ref.chromosomes, corresp.acc2ref$i.ref)
  }
  pokaz(ref.chromosomes)
}
ref.chromosomes = sort(ref.chromosomes)
pokaz('Reference chromosomes:', ref.chromosomes)


# Cut out every piece that comes from one source chromosome of `acc`.
# The chromosome is read once here. Previously it was re-read inside the
# per-piece loop, so the cost was `pieces x chromosome length` rather than the
# sum of the piece lengths -- on a heavily rearranged genome (e.g. Mytilus,
# 39-83 pieces per chromosome) that is 786 reads of a 118 Mb FASTA per accession
# where 14 are enough.
# Pieces are cut with substr() on the string instead of seq2nt()+indexing, which
# avoids materialising a 100M-element character vector per chromosome. Case is
# preserved on forward pieces and uppercased on reversed ones, exactly as before.
cut.acc.chr <- function(i.acc, corresp, acc, path.chr){
  file.chr = paste0(path.chr, acc, '_chr', i.acc, '.fasta')
  checkFile(file.chr)
  s.chr = readFasta(file.chr)[[1]]

  idx = which(corresp$i.acc == i.acc)
  pieces = lapply(idx, function(irow){
    s.tmp = substr(s.chr, corresp$beg[irow], corresp$end[irow])
    if(corresp$pos[irow] < 0) s.tmp = unname(revCompl(s.tmp))
    return(s.tmp)
  })
  names(pieces) = as.character(idx)       # position of the piece in the genome
  return(pieces)
}

# ---- Per-accession worker: build every chromosome of one genome ----
# Work is organised by source chromosome (one read each), not by reference
# chromosome, so no source FASTA is ever read twice.
build.acc.genome <- function(acc, path.chr, parallel.chr = F){
  corresp = readRDS(paste0(path.processed, 'corresp_',acc,'_to_',ref, '.rds'))
  corresp = corresp[corresp$i.ref %in% ref.chromosomes,]
  if(nrow(corresp) == 0) return(NULL)
  # Reference order: chromosome by chromosome, pieces by their position within it
  corresp = corresp[order(corresp$i.ref, abs(corresp$pos)),]

  i.acc.uniq = unique(corresp$i.acc)
  if(parallel.chr){
    cut.res = foreach(i.acc = i.acc.uniq,
                      .packages = c('pannagram', 'crayon')) %dopar% {
                cut.acc.chr(i.acc, corresp, acc, path.chr)
              }
  } else {
    cut.res = lapply(i.acc.uniq, function(i.acc) cut.acc.chr(i.acc, corresp, acc, path.chr))
  }

  pieces = unlist(cut.res, recursive = F)               # named by row of `corresp`
  rm(cut.res)
  pieces = pieces[as.character(1:nrow(corresp))]        # back to reference order

  # Glue each reference chromosome in one pass instead of growing a vector
  genome.ref = sapply(split(1:nrow(corresp), corresp$i.ref),
                      function(idx) paste0(pieces[idx], collapse = ''))
  rm(pieces)
  names(genome.ref) = paste0(acc, '_Chr', names(genome.ref))
  return(genome.ref)
}

process.acc.02 <- function(acc, parallel.chr = F){
  pokaz('Accession', acc)
  genome.ref = build.acc.genome(acc, path.chr, parallel.chr = parallel.chr)
  if(is.null(genome.ref)) return(invisible(NULL))
  writeFasta(genome.ref, paste0(path.processed, acc, '.fasta'))
  return(invisible(NULL))
}

# ---- Parallel backend ----
# FORK workers are a copy of this process as it is now, so the cluster is
# created after all the functions above exist -- otherwise a worker cannot see
# a function that is not named in the %dopar% expression itself (cut.acc.chr).
if(num.cores > 1){
  myCluster <- makeCluster(num.cores, type = "FORK")
  registerDoParallel(myCluster)
}

# ---- Main loop ----
# Spend the cores along whichever dimension is larger, without nesting:
#   * many accessions    -> one genome per worker (memory bounded to one genome);
#   * a single accession -> parallel over that genome's source chromosomes.
if(num.cores > 1 && length(accessions) > 1){
  foreach(acc = accessions, .packages = c('pannagram', 'crayon')) %dopar% process.acc.02(acc)
} else if(num.cores > 1){
  for(acc in accessions) process.acc.02(acc, parallel.chr = T)
} else {
  for(acc in accessions) process.acc.02(acc)
}


if(num.cores > 1) stopCluster(myCluster)

# Also copy the reference genome
genome = c()
for(i.chr.ref in ref.chromosomes){
  file.chr = paste0(path.chr, ref, '_chr', i.chr.ref, '.fasta')
  checkFile(file.chr)
  genome[paste0(ref, '_Chr',i.chr.ref)] = readFasta(file.chr)
}
writeFasta(genome, paste0(path.processed, ref, '.fasta'))

