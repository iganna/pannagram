# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
})

source(system.file("utils/utils.R", package = "pannagram"))
# source(system.file("pangen/synteny_funcs.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.cons",  type = "character", default = NULL, help = "Path to consensus directory"),
  make_option("--cores",      type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",   type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",  type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

max.len.gap = 100000  # This should be a parameter!!!

# ***********************************************************************
# ---- Logging ----
source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----
source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

aln.type.in = aln.type.comb
aln.type.out = aln.type.clean

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

s.pattern <- paste0("^", aln.type.in, ".*\\.*h5$")
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.in, "", files)
pref.combinations <- sub(".h5", "", pref.combinations)

if(length(pref.combinations) == 0) {
  stop('No files with the ref-based alignments are found')
}

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(s.comb,
                          echo.loop=T){
  
  # Log files
  file.log.loop = paste0(path.log, 'loop_file_', 
                         s.comb,
                         '.log')
  if(!file.exists(file.log.loop)){
    invisible(file.create(file.log.loop))
  }
  
  # ---- Check log Done ----
  if(checkDone(file.log.loop)){
    return()
  }
  
  # --- --- --- --- --- --- --- --- --- --- ---
  pokaz('Combination', s.comb)
  file.comb.in = paste0(path.cons, aln.type.in, s.comb,'.h5')
  file.comb.out = paste0(path.cons, aln.type.out, s.comb,'.h5')
  
  if(!file.exists(file.comb.in)){
    stop('File with combined references doesn’t exist')
  }
  
  # Create the output file
  suppressMessages({
    file.remove(file.comb.out)
    h5createFile(file.comb.out)
    h5createGroup(file.comb.out, gr.blocks)
    h5createGroup(file.comb.out, gr.accs.e)
  })
  
  idx.trust = h5read(file.comb.in, v.idx.trust)
  idx.trust = idx.trust != 0
  groups = h5ls(file.comb.in)
  accessions = groups$name[groups$group == gr.accs.b]
  
  # ---- Cleanup ----
  pokaz('Cleanup..')
  idx.nonzero = 0
  for(acc in accessions){
    pokaz('Accession', acc)
    
    s.acc = paste0(gr.accs.e, acc)
    v = h5read(file.comb.in, s.acc)
    v = v[idx.trust]
    
    v.init = v
    
    # Define blocks
    for(i in 1:2){
      v = v.init
    
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
      
      v.b.remove = v.b[v.b$len <= 100,]
      if(nrow(v.b.remove) == 0) break
      for(irow in 1:nrow(v.b.remove)){
        v.init[v.b.remove$i.beg[irow]:v.b.remove$i.end[irow]] = 0
      }
    }
    
    suppressMessages({
      h5write(v.init, file.comb.out, s.acc)
    })
  
    idx.nonzero = idx.nonzero + (abs(v.init) > 0) * 1
  }
  
  # ---- Remove zeros ----
  pokaz('Remove zeros..')
  idx.nonzero = idx.nonzero > 0
  # pokaz(length(idx.nonzero), sum(idx.nonzero))
  
  for(acc in accessions){
    pokaz('Accession', acc)
    
    s.acc = paste0(gr.accs.e, acc)
    v = h5read(file.comb.out, s.acc)
    
    v = v[idx.nonzero]
    
    # Rewrite  
    suppressMessages({
      h5delete(file.comb.out, s.acc)
      h5write(v, file.comb.out, s.acc)
    })
  }
  
  # ---- Breaks ----
  pokaz('Find breaks..')
  idx.breaks = c()
  for(acc in accessions){
    pokaz('Accession', acc)
    
    s.acc = paste0(gr.accs.e, acc)
    v = h5read(file.comb.out, s.acc)

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
    
    v.b = v.b[order(abs(v.b$v.beg)),]
    
    blocks.acc = rep(0, max(abs(v)))
    for(irow in 1:nrow(v.b)){
      blocks.acc[abs(v.b$v.beg[irow]):abs(v.b$v.end[irow])] = irow
    }
    
    suppressMessages({
      s.acc = paste0(gr.blocks, acc)
      h5write(blocks.acc, file.comb.out, s.acc)
    })
    
    # Find breaks
    i.br.acc = which(abs(diff(v)) != 1)
    if(length(i.br.acc) == 0) next
    
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
  for(s.comb in pref.combinations){
    loop.function(s.comb, echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(s.comb = pref.combinations, 
                .packages=c('rhdf5', 'crayon'))  %dopar% { 
                  loop.function(s.comb, echo.loop=echo.loop)
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

