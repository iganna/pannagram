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


# ---- Parallel backend (comb-style; FORK inherits all vars/functions) ----
if(num.cores > 1){
  myCluster <- makeCluster(num.cores, type = "FORK")
  registerDoParallel(myCluster)
}

# Build one reference-ordered chromosome for accession `acc`.
# Fast path: a single source chromosome taken whole and forward is copied as a
# string (no per-nucleotide seq2nt/nt2seq) -> avoids building multi-gigabase
# character vectors. Only reversed/spliced pieces are expanded nucleotide-wise.
build.ref.chr <- function(i.chr.ref, corresp.acc2ref, acc, acc.len, path.chr){
  corresp.tmp = corresp.acc2ref[corresp.acc2ref$i.ref == i.chr.ref,]
  if(nrow(corresp.tmp) == 0) return(NULL)
  corresp.tmp = corresp.tmp[order(abs(corresp.tmp$pos)),]

  read.chr <- function(i.acc){
    file.chr = paste0(path.chr, acc, '_chr', i.acc, '.fasta')
    checkFile(file.chr)
    readFasta(file.chr)[[1]]
  }

  if(nrow(corresp.tmp) == 1 &&
     corresp.tmp$beg[1] == 1 &&
     corresp.tmp$end[1] == acc.len$len[corresp.tmp$i.acc[1]] &&
     corresp.tmp$pos[1] >= 0){
    out = read.chr(corresp.tmp$i.acc[1])                 # whole chromosome, forward
  } else {
    s = c()
    for(irow in 1:nrow(corresp.tmp)){
      s.tmp = seq2nt(read.chr(corresp.tmp$i.acc[irow]))[corresp.tmp$beg[irow]:corresp.tmp$end[irow]]
      if(corresp.tmp$pos[irow] < 0) s.tmp = revCompl(s.tmp)
      s = c(s, s.tmp)
    }
    out = nt2seq(s)
  }
  names(out) = paste0(acc, '_Chr', i.chr.ref)
  return(out)
}

# ---- Per-accession worker: build every chromosome of one genome and write it ----
process.acc.02 <- function(acc){
  pokaz('Accession', acc)
  file.acc.len = paste0(path.chr, acc, '_chr_len.txt', collapse = '')
  acc.len = read.table(file.acc.len, header = 1)
  corresp.acc2ref = readRDS(paste0(path.processed, 'corresp_',acc,'_to_',ref, '.rds'))
  seqs = lapply(ref.chromosomes, function(i.chr.ref)
                build.ref.chr(i.chr.ref, corresp.acc2ref, acc, acc.len, path.chr))
  genome.ref = unlist(seqs)                              # named vector; NULLs dropped
  writeFasta(genome.ref, paste0(path.processed, acc, '.fasta'))
  return(invisible(NULL))
}

# ---- Main loop ----
# Spend the cores along whichever dimension is larger, without nesting:
#   * many accessions   -> one genome per worker (memory bounded to one genome);
#   * a single accession -> parallel over that genome's chromosomes.
if(num.cores > 1 && length(accessions) > 1){
  foreach(acc = accessions, .packages = c('pannagram', 'crayon')) %dopar% process.acc.02(acc)
} else if(num.cores > 1){
  for(acc in accessions){
    pokaz('Accession', acc)
    file.acc.len = paste0(path.chr, acc, '_chr_len.txt', collapse = '')
    acc.len = read.table(file.acc.len, header = 1)
    corresp.acc2ref = readRDS(paste0(path.processed, 'corresp_',acc,'_to_',ref, '.rds'))
    seqs = foreach(i.chr.ref = ref.chromosomes,
                   .packages = c('pannagram', 'crayon')) %dopar% {
             build.ref.chr(i.chr.ref, corresp.acc2ref, acc, acc.len, path.chr)
           }
    genome.ref = unlist(seqs)
    writeFasta(genome.ref, paste0(path.processed, acc, '.fasta'))
  }
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

