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
# ---- Main loop ----
for(acc in accessions){
  pokaz('Accession', acc)
  
  file.acc.len = paste0(path.chr, acc, '_chr_len.txt', collapse = '')
  acc.len = read.table(file.acc.len, header = 1)
  
  
  print(corresp.acc2ref)
  corresp.acc2ref = readRDS(paste0(path.processed, 'corresp_',acc,'_to_',ref, '.rds'))
  
  n.acc = max(corresp.acc2ref$i.acc)
  n.ref = max(corresp.acc2ref$i.ref)

  genome = c()
  for(i.chr.acc in 1:n.acc){
    file.chr = paste0(path.chr, acc, '_Chr', i.chr, '.fasta')
    checkFile(file.chr)
    genome[i.chr.acc] = seq2nt(readFasta(file.chr))
  }
  

  genome.ref = c()
  for(i.chr.ref in 1:n.ref){
    corresp.tmp = corresp.acc2ref[corresp.acc2ref$i.ref == i.chr.ref,]
    if(nrow(corresp.tmp) == 0) next
    corresp.tmp = corresp.tmp[order(abs(corresp.tmp$pos)),]
    s = c()
    for(irow in 1:nrow(corresp.tmp)){
      s.tmp = genome[corresp.tmp$i.acc[irow]][corresp.tmp$beg[irow]:corresp.tmp$end[irow]]
      if(corresp.tmp$pos[irow] < 0){
        s.tmp = revCompl(s.tmp)
      }
      s = c(c, s.tmp)
    }
    genome.ref[i.chr.ref] = nt2seq(s)
  }
  
  writeFasta(genome.ref, paste0(path.processed, acc, '.fasta'))
  
}


