# Alignment-1. Remaining syntenic (major) matches

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(pannagram)
})

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("--path.aln"),      type = "character", default = NULL, help = "Path to the output directory with alignments"),
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

path.aln      <- ifelse(!is.null(opt$path.aln), opt$path.aln, stop('Folder with Alignments is not specified'))
path.resort   <- ifelse(!is.null(opt$path.resort), opt$path.resort, stop('Folder combinations to resort the genomes'))
base.acc      <- ifelse(!is.null(opt$ref), opt$ref, stop('Reference genome is not specified'))

pokaz('path resort', path.resort)

# Create folders for the alignment results
if(!dir.exists(path.resort)) dir.create(path.resort)

# Accessions
files.aln <- list.files(path.aln, pattern = ".*\\.rds$", full.names = F)
files.aln = sub("_maj.rds", "", files.aln)
accessions = c()
for(f in files.aln){
  s = strsplit(f, '_')[[1]]
  s = s[-(length(s)-(0:1))]
  s = paste0(s, collapse = '_')
  accessions = c(accessions, s)
  accessions = unique(accessions)
}

# ***********************************************************************
# ---- Correspondence ----

for(acc in accessions){
  pokaz('Accession', acc)
  files.aln <- list.files(path.aln, pattern = paste0(acc,".*\\.rds$"), full.names = F)
  # pokaz(files.aln)
  files.sizes <- file.size(paste0(path.aln, files.aln))
  combinations = sub(paste0(acc, "_"), "", files.aln)
  combinations = sub("_maj.rds", "", combinations)
  
  result <- (do.call(rbind, strsplit(combinations, "_")))
  result <- data.frame(Col1 = as.numeric(result[, 1]), Col2 = as.numeric(result[, 2]))
  
  n.acc = max(result[,1])
  n.ref = max(result[,2])
  
  chr.acc.blocked = c()
  corresp = c()
  for(i.chr.ref in 1:n.ref){
    # pokaz(i.chr.ref)
    idx.tmp = grep(paste0("_", i.chr.ref, "_maj\\.rds$"), files.aln)
    idx.corr = idx.tmp[which.max(files.sizes[idx.tmp])]
    i.chr.acc = result[idx.corr,1]
    
    if(length(i.chr.acc) == 0){
      pokazAttention("Cannot get the correspondence for chromosome", i.chr.ref)
      next
    }
    file.aln = paste0(path.aln, files.aln[idx.corr])
    # pokaz(file.aln)
    x = readRDS(file.aln)
    pos.plus = sum((x$V3 - x$V2 + 1) * (1 - x$dir))
    
    pos.minus = sum((x$V3 - x$V2 + 1) * x$dir)
    
    if(pos.minus > pos.plus){
      i.chr.acc = -i.chr.acc
    }
    
    corresp = rbind(corresp, 
                    c(i.chr.ref, i.chr.acc))
  }
  
  if(length(unique(abs(corresp[,2]))) != length(corresp[,2])){
    pokazAttention("Cannot get the correspondence for some chromosomes for the accession", acc)
    next
  }
  
  if(n.acc != length(corresp[,2])){
    pokazAttention("Not all of the achromosomes have a correspondence, accession", acc)
    next
  }
  
  file.resort = paste0(path.resort, acc, '.rds')  
  saveRDS(corresp, file.resort)
  
}




