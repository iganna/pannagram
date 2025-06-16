#' How the output files look like:
#'     group           name       otype   dclass      dim
#'         /           accs   H5I_GROUP                  
#'     /accs              0 H5I_DATASET    FLOAT 28940631
#'     /accs          10002 H5I_DATASET    FLOAT 28940631
#'     /accs          10015 H5I_DATASET    FLOAT 28940631
#'     
#' In combo-files the reference accession column will not have zeros, because it will participate further lika a function.

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/synteny_func.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("--path.cons"),     type = "character", default = NULL, help = "Path to consensus directory"),
  make_option(c("--path.aln"),      type = "character", default = NULL, help = "Path to the output directory with alignments"),
  make_option(c("--path.chr"),      type = "character", default = NULL, help = "Path to the reference file"),
  
  make_option(c("--ref"),          type = "character", default = NULL, help = "Prefix of the reference file"),
  make_option(c("--accessions"),   type = "character", default = NULL, help = "File containing accessions to analyze"),
  make_option(c("--combinations"), type = "character", default = NULL, help = "File containing combinations to analyze"),
  
  make_option(c("--cores"),        type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option(c("--path.log"),     type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"),    type = "character", default = NULL, help = "Level of log to be shown on the screen")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

# ***********************************************************************
# ---- Values of parameters ----

# print(opt)

# Number of cores
num.cores <- opt$cores

# Reference accessions
if (!is.null(opt$ref)) base.acc.ref <- opt$ref
pokaz('Reference genome', base.acc.ref, file=file.log.main, echo=echo.main)

# Paths
if (!is.null(opt$path.chr)) path.chr <- opt$path.chr  # to know the chromosomal lengths
if (!is.null(opt$path.aln)) path.aln <- opt$path.aln # Path with alignments
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons # Path with the consensus output

# Create directories
if(!dir.exists(path.cons)) system(paste0('mkdir ', path.cons))

# ***********************************************************************

# Accessions
file.acc <- ifelse(!is.null(opt$accessions), opt$accessions, stop("File with accessions are not specified"))
tmp <- read.table(file.acc, stringsAsFactors = F)
accessions <- as.character(tmp[,1])
pokaz('Names of genomes for the analysis:', accessions, 
      file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- Combinations ----

# Available combinations

files.aln <- list.files(path = path.aln, pattern = "_full.rds$")
pokaz('Full files:', files.aln, file=file.log.main, echo=echo.main)
files.aln <- files.aln[sapply(files.aln, function(x) any(sapply(accessions, function(a) startsWith(x, a))))]
files.aln = gsub("_full.rds", "", files.aln)

chromosome.pairs = unique(sapply(files.aln, function(s){
  info = pref2info(s)
  return(paste(info$query.chr, info$base.chr, sep = '_'))
}))

file.combinations <- ifelse(!is.null(opt$combinations), opt$combinations, stop("File with combinations are not specified"))
if (length(readLines(file.combinations)) != 0) {
  combinations.filter = read.table(file.combinations)
  
  # String format
  combinations.filter.str <- apply(combinations.filter, 1, paste, collapse = "_")
  
  # Intersect
  chromosome.pairs <- intersect(chromosome.pairs, combinations.filter.str)
}

pokaz('Combinations:', chromosome.pairs, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- Length of reference chromosomes ----

file.chr.len = paste0(path.chr, base.acc.ref, '_chr_len.txt')
if(!file.exists(file.chr.len)) stop('File with chromosomes does not exist')
chr.len = read.table(file.chr.len, stringsAsFactors = F, header = 1)
chr.len = chr.len$len

pokaz('Chromosomal lengths:', chr.len, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(s.comb, 
                          echo.loop=T){
  
  initial.vars <- ls()
  
  s.comb = strsplit(s.comb, '_')[[1]]
  query.chr = as.numeric(s.comb[1])
  base.chr = as.numeric(s.comb[2])
  
  # Log files
  file.log.loop = paste0(path.log, 'loop_file_', 
                         query.chr, '_', base.chr, '_', base.acc.ref,
                         '.log')
  if(!file.exists(file.log.loop)){
    invisible(file.create(file.log.loop))
  }
  
  # ---- Check log Done ----
  if(checkDone(file.log.loop)){
    return()
  }
  
  pokaz('Combination:', query.chr, base.chr, file=file.log.loop, echo=echo.loop)
  pokaz('Chromosomal length', chr.len, file=file.log.loop, echo=echo.loop)
  base.len = chr.len[base.chr]
  
  file.comb = paste0(path.cons, aln.type.ref, query.chr, '_', base.chr, '_', base.acc.ref,'.h5')
  if (file.exists(file.comb)) file.remove(file.comb)
  h5createFile(file.comb)
  
  # Path to accessions chunks
  # TODO: Check the availability of the group before creating it
  h5createGroup(file.comb, gr.accs.e)
  
  for(acc in accessions){
    
    pokaz('Accession', acc, 'qchr', query.chr, 'bchr', base.chr, 
                   file=file.log.loop, echo=echo.loop)
    
    pref.comb = paste0(acc, '_', query.chr, '_', base.chr, collapse = '')
    file.aln.full <- paste(path.aln, paste0(pref.comb,  '_full.rds', collapse = ''), sep = '')
    if(!file.exists(file.aln.full)) next
    
    pokaz('Alignment file:', file.aln.full, file=file.log.loop, echo=echo.loop)
    
    # Reading the alignment
    x = readRDS(file.aln.full)
    
    pokaz('Base len', base.len, file=file.log.loop, echo=echo.loop)
    # saveRDS(x, 'tmp.rds')
    
    # Get query coordinates in base order
    x.corr = getCorresp2BaseSign(x, base.len)
    
    if(sum(duplicated(x.corr[x.corr != 0])) > 0) stop('DUPLICSTIONS', sum(duplicated(x.corr[x.corr != 0])))
    
    # Write into file
    suppressMessages({
      h5write(x.corr, file.comb, paste0(gr.accs.e, '', acc))
    })
    
    
    rmSafe(x.corr)
    rmSafe(x)
    rmSafe(v)
    rmSafe(idx.tmp.acc)
    
  }
  
  suppressMessages({
    
    h5write(base.acc.ref, file.comb, v.ref.name)
    h5write(base.len, file.comb, v.len)
    
    h5write(1:base.len, file.comb, paste0(gr.accs.e, '', base.acc.ref))
  })
  
  pokaz('Done.', file=file.log.loop, echo=echo.loop)
  
  H5close()
  
  final.vars <- ls()
  new.vars <- setdiff(final.vars, initial.vars)
  rm(list = new.vars)
  gc()
  
  return(NULL)
}

# ***********************************************************************
# ---- Loop  ----

if(num.cores == 1){
  # file.log.loop = paste0(path.log, 'loop_all.log')
  # invisible(file.create(file.log.loop))
  for(s.comb in chromosome.pairs){
    loop.function(s.comb,
                  # file.log.loop = file.log.loop, 
                  echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(s.comb = chromosome.pairs, 
                .packages=c('rhdf5', 'crayon'))  %dopar% { 
                  loop.function(s.comb,
                                echo.loop=echo.loop)
                }
  stopCluster(myCluster)
}

warnings()

pokaz('Done.',
      file=file.log.main, echo=echo.main)

