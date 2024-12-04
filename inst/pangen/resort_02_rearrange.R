# Alignment-1. Remaining syntenic (major) matches

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(pannagram)
  library(crayon)
})

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("--path.genomes"),      type = "character", default = NULL, help = "Path to the output directory with genomes"),
  make_option(c("--path.new"),      type = "character", default = NULL, help = "Path to the output directory with rearranged genomes"),
  make_option(c("--path.resort"),   type = "character", default = NULL, help = "Path to the output directory with alignments"),
  
  make_option(c("--ref"),           type = "character", default = NULL, help = "Name of the reference genome"),
  
  make_option(c("--cores"),         type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option(c("--path.log"),      type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"),     type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Variables ----

max.len = 10^6
len.blast = 50

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- Values of parameters ----

# Number of cores
num.cores <- opt$cores

path.genomes  <- ifelse(!is.null(opt$path.genomes), opt$path.genomes, stop('Folder with genomes is not specified'))
path.resort   <- ifelse(!is.null(opt$path.resort), opt$path.resort, stop('Folder combinations to resort the genomes'))
base.acc      <- ifelse(!is.null(opt$ref), opt$ref, stop('Reference genome is not specified'))

path.new      <- opt$path.new
if(is.null(path.new)){
  path.new = paste0(path.genomes, 'rearranged/')
}
if(!dir.exists(path.new)) dir.create(path.new)

# Accessions
files.aln <- list.files(path.resort, pattern = ".*\\.rds$", full.names = F)
accessions = sub(".rds", "", files.aln)

if(length(unique(accessions)) != length(accessions)) stop('Something is wring with accession names')
accessions = setdiff(accessions, base.acc)

pokaz('Accessions:', accessions)

# ***********************************************************************
# ---- Correspondence ----

for(acc in accessions){
  pokaz('Accession', acc)

  file.resort = paste0(path.resort, acc, '.rds')
  corresp = readRDS(file.resort)
  corresp = corresp[order(corresp[,1]),]
  print(corresp)

  pokaz('Reading...')
  file.genome = paste0(paste0(path.genomes, acc, ".fasta"))
  pokaz(file.genome)
  genome = readFasta(file.genome)

  genome.new = rep(NA, nrow(corresp))
  for(irow in 1:nrow(corresp)){
    idx = corresp[irow,2]
    pokaz(idx)
    if(idx < 0){
      idx = abs(idx)
      genome.new[irow] = revComplSeq(genome[idx])
    } else {
      genome.new[irow] = genome[idx]
    }
  }
  names(genome.new) = paste0('Chr', 1:length(genome.new), "_reodered")

  pokaz('Writing...')
  file.genome.new = paste0(path.new, acc, '.fasta')
  writeFasta(genome.new, file.genome.new)

}

file.genome = paste0(path.genomes, base.acc, ".fasta")
file.genome.new = paste0(path.new, base.acc, '.fasta')

# pokaz(file.genome, file.genome.new)
res.copy = file.copy(file.genome, file.genome.new, overwrite = TRUE, copy.mode = TRUE)
if(!res.copy){
  pokazAttention("Somethig is wrong with copying the reference genome file")
}


