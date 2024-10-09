# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
})

source(system.file("utils/utils.R", package = "pannagram"))
# source("pangen/synteny_funcs.R")

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to consensus directory", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("--path.log"), type = "character", default = NULL,
              help = "Path for log files", metavar = "character"),
  make_option(c("--log.level"), type = "character", default = NULL,
              help = "Level of log to be shown on the screen", metavar = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

max.len.gap = 100000  # This should be a parameter!!!

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

# ***********************************************************************
# ---- Values of parameters ----

# Set the number of cores for parallel processing
num.cores.max = 10
num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))

# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop('Consensus folder doesn’t exist')

# ***********************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type.comb, ".*")
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.comb, "", files)
# pref.combinations <- sub("_ref.*$", "", pref.combinations)
pref.combinations <- sub(".h5", "", pref.combinations)

if(length(pref.combinations) == 0) {
  stop('No files with the ref-based alignments are found')
}

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(s.comb,
                          echo.loop=T, 
                          file.log.loop=NULL){
  
  # Log files
  if (is.null(file.log.loop)){
    file.log.loop = paste0(path.log, 'loop_file_', 
                           s.comb,
                           '.log')
    invisible(file.create(file.log.loop))
  }
  
  # ---- Check log Done ----
  if(checkDone(file.log.loop)){
    return()
  }
  
  # --- --- --- --- --- --- --- --- --- --- ---
  
  file.comb = paste0(path.cons, aln.type.comb, s.comb,'.h5')
  
  if(!file.exists(file.comb)){
    stop('File with combined references doesn’t exist')
  }
  
  suppressMessages({
  h5createGroup(file.comb, gr.blocks)
  })
  
  idx.trust = h5read(file.comb, v.idx.trust)
  idx.trust = idx.trust != 0
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  
  idx.breaks = c()
 
  for(acc in accessions){
    pokaz(acc)
    s.acc = paste0(gr.accs.e, acc)
    v = h5read(file.comb, s.acc)
    v = v[idx.trust]
      
    suppressMessages({
      h5delete(file.comb, s.acc)
      h5write(v, file.comb, s.acc)
    })
    
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
    
    blocks.acc = rep(0, max(v))
    for(irow in 1:nrow(v.b)){
      blocks.acc[abs(v.b$v.beg[irow]):abs(v.b$v.end[irow])] = irow
    }
    
    suppressMessages({
      s.acc = paste0(gr.blocks, acc)
      h5write(blocks.acc, file.comb, s.acc)
    })
    
    
    # Find breaks
    i.br.acc = which(abs(diff(v)) != 1)
    df = data.frame(val.beg = v[i.br.acc],
                    val.end = v[i.br.acc+1],
                    idx.beg = v.idx[i.br.acc],
                    idx.end = v.idx[i.br.acc+1])
    
    df = df[blocks.acc[abs(df$val.beg)] == blocks.acc[abs(df$val.end)],]
    
    df$acc = acc
    df$len.acc = abs(df$val.end - df$val.beg) - 1
    df$len.comb = abs(df$idx.end - df$idx.beg) - 1
    
    idx.breaks = rbind(idx.breaks, df)
    
  }
  
  # Update index to trust
  suppressMessages({
    idx.trust = h5read(file.comb, v.idx.trust)
    idx.trust = idx.trust[idx.trust != 0]
    
    h5delete(file.comb, v.idx.trust)
    h5write(idx.trust, file.comb, v.idx.trust)
  })
  
  
  file.breaks = paste0(path.cons, 'breaks_', s.comb,'.rds')
  saveRDS(idx.breaks, file.breaks)
  
  rmSafe(idx.breaks)
  
  H5close()
  gc()
  
  pokaz('Done.', file=file.log.loop, echo=echo.loop)
  return(NULL)
}


# ***********************************************************************
# ---- Loop  ----

if(num.cores == 1){
  # file.log.loop = paste0(path.log, 'loop_all.log')
  # invisible(file.create(file.log.loop))
  for(s.comb in pref.combinations){
    loop.function(s.comb,
                  # file.log.loop = file.log.loop, 
                  echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(s.comb = pref.combinations, 
                .packages=c('rhdf5', 'crayon'))  %dopar% { 
                  loop.function(s.comb,echo.loop=echo.loop)
                }
  stopCluster(myCluster)
}

warnings()

pokaz('Done.',
      file=file.log.main, echo=echo.main)


# ***********************************************************************
# ---- Testing ----
if(F){
  library(rhdf5)
  path.cons = './'
  options("width"=200, digits=10)
}

