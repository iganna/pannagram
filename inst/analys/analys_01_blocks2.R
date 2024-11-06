# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
})

source(system.file("utils/utils.R", package = "pannagram"))

# Define blocks in the alignemnt

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.cons", type = "character", default = NULL, help = "Path to consensus directory"),
  make_option("--ref",       type = "character", default = "",   help = "Prefix of the reference file"),
  make_option("--cores",     type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--aln.type",  type = "character", default = NULL, help = "Prefix for the output file")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

max.len.gap = 20000

# print(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

aln.type <- ifelse(is.null(opt$aln.type), aln.type.msa, opt$aln.type)


# ***********************************************************************
# ---- Values of parameters ----

# Set the number of cores for parallel processing
num.cores <- opt$cores

# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop('Consensus folder doesn’t exist')

# Reference genome
ref.name <- opt$ref

# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type, ".*h5")
s.combinations <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
s.combinations = gsub(aln.type, "", s.combinations)
s.combinations = gsub(".h5", "", s.combinations)

if(ref.name != ""){
  ref.suff = paste0('_', ref.name)
  
  pokaz('Reference:', ref.name)
  s.combinations <- s.combinations[grep(ref.suff, s.combinations)]
  s.combinations = gsub(ref.suff, "", s.combinations)
  
} else {
  ref.suff = ''
}

if(length(s.combinations) == 0){
  stop('No Combinations found.')
} else {
  pokaz('Combinations', s.combinations)  
}
# ***********************************************************************
# ---- MAIN program body ----


# path.cons = '/Volumes/Samsung_T5/vienn/test/a600/'
# aln.type = aln.type.ref
# ref.suff = '_0'

loop.function <- function(s.comb, echo = T){
  
  i.chr = as.numeric(strsplit(s.comb, '_')[[1]][1])
  pokaz('Chromosome', i.chr)
  # --- --- --- --- --- --- --- --- --- --- ---
  
  file.comb.in = paste0(path.cons, aln.type, s.comb, ref.suff,'.h5')

  groups = h5ls(file.comb.in)
  accessions = groups$name[groups$group == gr.accs.b]
  
  idx.breaks = c()
  for(acc in accessions){
    initial.vars <- ls()
    pokaz('Accession', acc)
    s.acc = paste0(gr.accs.e, acc)
    v = h5read(file.comb.in, s.acc)
    
    v.init = v
    # Define blocks
    v.idx = 1:length(v)
    
    v.idx = v.idx[v != 0]
    v = v[v != 0]
    v.r = rank(abs(v))
    v.r[v < 0] = v.r[v < 0] * (-1)
    v.b = findRuns(v.r)
    
    v.b$v.beg = v[v.b$beg]
    v.b$v.end = v[v.b$end]
    
    v.b$i.beg = v.idx[v.b$beg]
    v.b$i.end = v.idx[v.b$end]
    
    df.res = data.frame(pan.b = v.b$i.beg,
                        pan.e = v.b$i.end,
                        own.b = abs(v.b$v.beg),
                        own.e = abs(v.b$v.end),
                        acc = acc,
                        chr = i.chr,
                        dir = (1 - sign(v.b$v.beg)) / 2)
    
    idx.breaks = rbind(idx.breaks,
                       df.res)
    
    # Cleanup variables
    final.vars <- ls()
    new.vars <- setdiff(final.vars, initial.vars)
    rm(list = new.vars)
    gc()
  }
  
  file.blocks = paste0(path.cons, aln.type, 'bl_', s.comb, ref.suff, '.rds')
  saveRDS(idx.breaks, file.blocks, compress = F)
  
  rm(idx.breaks)
  
}

# ***********************************************************************
# ---- Loop  ----


if(num.cores == 1){
  for(s.comb in s.combinations){
    loop.function(s.comb)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  list.blocks = foreach(s.comb = s.combinations, .packages=c('rhdf5', 'crayon'))  %dopar% { 
    loop.function(s.comb)
    return(tmp)
  }
  stopCluster(myCluster)
}


warnings()
