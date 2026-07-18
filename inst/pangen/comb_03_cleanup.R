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
  make_option("--path.features.msa",  type = "character", default = NULL, help = "Path to msa directory (features)"),
  make_option("--path.inter.msa",     type = "character", default = NULL, help = "Path to msa directory (internal)"),
  make_option("--accessions",         type = "character", default = NULL, help = "File containing accessions to analyze"),
  make_option("--cores",              type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",           type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",          type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

max.len.gap = 100000  # This should be a parameter!!!
min.block.len = 100

# ***********************************************************************
# ---- Logging ----
source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----
source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

aln.type.in = paste0(aln.type.comb, '_')
aln.type.out = paste0(aln.type.clean, '_')

# ***********************************************************************
# ---- Accessions ----

file.acc <- ifelse(!is.null(opt$accessions), opt$accessions, stop("File with accessions are not specified"))
tmp <- read.table(file.acc, stringsAsFactors = F)
accessions.specified <- as.character(tmp[,1])

# ***********************************************************************
# ---- Values of parameters ----

# Set the number of cores for parallel processing
num.cores.max = 10
num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))

# Path with the consensus output
if (!is.null(opt$path.features.msa)) path.features.msa <- opt$path.features.msa
if(!dir.exists(path.features.msa)) stop('path_features_msa directory does not exist')

if (!is.null(opt$path.inter.msa)) path.inter.msa <- opt$path.inter.msa
if(!dir.exists(path.inter.msa)) stop('path_inter_msa folder does not exist')

# ***********************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type.in, ".*\\.*h5$")
files <- list.files(path = path.inter.msa, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.in, "", files)
pref.combinations <- sub(".h5", "", pref.combinations)

if(length(pref.combinations) == 0) {
  stop('No files with the ref-based alignments are found')
}

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(s.comb,
                          done.set = character(0),
                          echo.loop=T){

  # ---- Checkpoint: skip already completed combinations ----
  s.comb.id <- s.comb
  if(s.comb.id %in% done.set){
    return(NULL)
  }

  # One log file per worker (bounded number of files); also the checkpoint ledger
  file.log.loop = initLoopLog(path.log)

  # --- --- --- --- --- --- --- --- --- --- ---
  pokaz('Combination', s.comb, file=file.log.loop, echo=echo.loop)
  file.comb.in = paste0(path.inter.msa, aln.type.in, s.comb,'.h5')
  file.comb.out = paste0(path.features.msa, aln.type.out, s.comb,'.h5')
  
  if(!file.exists(file.comb.in)){
    stop('File with combined references does not exist')
  }
  
  # Create the output file
  # suppressMessages({
  #   file.remove(file.comb.out)
  #   h5createFile(file.comb.out)
  #   h5createGroup(file.comb.out, gr.blocks)
  #   h5createGroup(file.comb.out, gr.accs.e)
  # })
  if(!file.exists(file.comb.out)){
      h5createFile(file.comb.out)
      h5createGroup(file.comb.out, gr.blocks)
      h5createGroup(file.comb.out, gr.accs.e)
  }
  
  idx.trust = h5read(file.comb.in, v.idx.trust)
  idx.trust = idx.trust != 0
  groups = h5ls(file.comb.in)
  accessions = groups$name[groups$group == gr.accs.b]
  accessions = intersect(accessions, accessions.specified)
  
  # ---- Cleanup ----
  pokaz('Cleanup..', file=file.log.loop, echo=echo.loop)
  idx.nonzero = 0
  for(acc in accessions){
    
    pokaz('Accession cleanup', acc, file=file.log.loop, echo=echo.loop)
    
    s.acc = paste0(gr.accs.e, acc)
    v = h5read(file.comb.in, s.acc)
    v = v[idx.trust]
    
    v.init = v

    clean.id = paste0(s.comb.id, '_clean_', acc)
    if(!(clean.id %in% done.set)){
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
        
        v.b.remove = v.b[v.b$len <= min.block.len,]
        if(nrow(v.b.remove) == 0) break
        idx.rm = sequence(v.b.remove$i.end - v.b.remove$i.beg + 1L,
                          from = v.b.remove$i.beg)
        v.init[idx.rm] = 0
      }
      
      suppressMessages({
        try(h5delete(file.comb.out, s.acc), silent = TRUE)
        h5write(v.init, file.comb.out, s.acc)
      })

      markDone(clean.id, file=file.log.loop, echo=echo.loop)
    } else {
      # Cleaning was done in a prior run: read the CLEANED vector back from the
      # output (already trust-subset + cleaned) so idx.nonzero -- the survive-mask
      # for the whole combination -- is built from cleaned, not raw, data on resume.
      v.init = h5read(file.comb.out, s.acc)
    }

    idx.nonzero = idx.nonzero + (abs(v.init) > 0) * 1
  }
  
  # ---- Remove zeros ----
  pokaz('Remove zeros..', file=file.log.loop, echo=echo.loop)
  idx.nonzero = idx.nonzero > 0
  # pokaz(length(idx.nonzero), sum(idx.nonzero))
  
  
  for(acc in accessions){
    pokaz('Accession', acc, file=file.log.loop, echo=echo.loop)
    
    zero.id = paste0(s.comb.id, '_zero_', acc)
    if(!(zero.id %in% done.set)) {
      s.acc = paste0(gr.accs.e, acc)
      v = h5read(file.comb.out, s.acc)

      v = v[idx.nonzero]

      # Rewrite (idempotent)
      suppressMessages({
        try(h5delete(file.comb.out, s.acc), silent = TRUE)
        h5write(v, file.comb.out, s.acc)
      })
      markDone(zero.id, file=file.log.loop, echo=echo.loop)
    }
  }
  
  # ---- Breaks ----
  pokaz('Find breaks..', file=file.log.loop, echo=echo.loop)
  # NOTE: kept as an iterative rbind on purpose. do.call(rbind,list) is faster but
  # R disambiguates duplicate data.frame rownames differently there, so the saved
  # breaks_*.rds would not be byte-identical to the original (values are the same).
  idx.breaks = c()
  for(acc in accessions){
    pokaz('Accession', acc, file=file.log.loop, echo=echo.loop)
    
    s.acc = paste0(gr.accs.e, acc)
    v = h5read(file.comb.out, s.acc)
    v.idx = 1:length(v)
    v.idx = v.idx[v != 0]
    v = v[v != 0]
    
    break.id = paste0(s.comb.id, '_break_', acc)
    if(!(break.id %in% done.set)){
      # Define blocks
      
      v.r = rank(abs(v))
      v.r[v < 0] = v.r[v < 0] * (-1)
      v.b = findRuns(v.r)
      
      v.b$v.beg = v[v.b$beg]
      v.b$v.end = v[v.b$end]
      
      v.b$i.beg = v.idx[v.b$beg]
      v.b$i.end = v.idx[v.b$end]
      
      v.b = v.b[order(abs(v.b$v.beg)),]
      
      if(length(v) == 0){
        # Fully-absent accession: no blocks. Placeholder so gr.blocks exists for
        # comb_04 (it is never indexed for an absent accession).
        blocks.acc = 0L
      } else {
        blocks.acc = rep(0, max(abs(v)))
        # abs(v.beg) may exceed abs(v.end) on minus-strand runs; the original loop
        # used beg:end (descending) which covers the same index set -> use lo/hi.
        b.lo = pmin(abs(v.b$v.beg), abs(v.b$v.end))
        b.hi = pmax(abs(v.b$v.beg), abs(v.b$v.end))
        b.lens = b.hi - b.lo + 1L
        idx.fill = sequence(b.lens, from = b.lo)
        blocks.acc[idx.fill] = rep(seq_len(nrow(v.b)), b.lens)
      }
      
      suppressMessages({
        s.acc = paste0(gr.blocks, acc)
        h5write(blocks.acc, file.comb.out, s.acc)
      })
      markDone(break.id, file=file.log.loop, echo=echo.loop)
    } else {
      s.acc = paste0(gr.blocks, acc)
      blocks.acc = h5read(file.comb.out, s.acc)
    }
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

  file.breaks = paste0(path.inter.msa, 'breaks_', s.comb,'.rds')
  saveRDS(idx.breaks, file.breaks)
  
  rmSafe(idx.breaks)
  
  H5close()
  gc()

  # ---- Checkpoint marker: combination fully processed ----
  markDone(s.comb.id, file=file.log.loop, echo=echo.loop)
  return(NULL)
}

# ***********************************************************************
# ---- Loop  ----

# Rebuild the set of completed combinations from all worker logs (any core count)
done.set <- getDoneSet(path.log)
if(length(done.set) > 0){
  pokaz('Skip already done:', length(done.set), file=file.log.main, echo=echo.main)
}

if(num.cores == 1){
  assign('.worker.id', 1, envir = .GlobalEnv)   # stable single file core_1.log
  for(s.comb in pref.combinations){
    loop.function(s.comb, done.set = done.set, echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK")
  registerDoParallel(myCluster)

  # Assign a stable worker id (1..N) to each worker -> bounded, reused file names
  parallel::clusterApply(myCluster, seq_len(num.cores),
                         function(i) assign('.worker.id', i, envir = .GlobalEnv))

  tmp = foreach(s.comb = pref.combinations,
                .packages=c('rhdf5', 'crayon'))  %dopar% {
                  loop.function(s.comb, done.set = done.set, echo.loop=echo.loop)
                }
  stopCluster(myCluster)
}


pokaz('Done.',
      file=file.log.main, echo=echo.main)



