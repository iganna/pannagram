suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.features.msa", type = "character", default = NULL, help = "Path to msa directory (features)"),
  make_option("--path.inter.msa",    type = "character", default = NULL, help = "Path to msa directory (internal)"),
  make_option("--accessions",        type = "character", default = NULL, help = "File containing accessions to analyze"),
  
  make_option("--max.len.gap",       type = "integer",   default = NULL, help = "Max length of the gap"),
  
  make_option("--cores",             type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",          type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",         type = "character", default = NULL, help = "Level of log to be shown on the screen"),
  make_option("--aln.type.in",       type = "character", default = NULL, help = "Name of the alignment file")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser, args = args)

# ***********************************************************************
# ---- Logging ----
source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----
source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

aln.type.in <- ifelse(is.null(opt$aln.type.in), aln.type.clean, opt$aln.type.in)
aln.type.in = paste0(aln.type.in, '_')

# ***********************************************************************
# ---- Accessions ----

file.acc <- ifelse(!is.null(opt$accessions), opt$accessions, stop("File with accessions are not specified"))
accessions.specified <- as.character(read.table(file.acc, stringsAsFactors = FALSE)[, 1])


# ***********************************************************************
# ---- Values of parameters ----

# Max len gap
if (is.null(opt$max.len.gap)) {
  stop("Error: max.len.gap is NULL")
} else {
  len.large <- opt$max.len.gap
}

# Number of cores for parallel processing
num.cores = opt$cores
pokaz('Number of cores', num.cores, file=file.log.main, echo=echo.main)

# Path with the MSA output (features)
path.features.msa <- opt$path.features.msa
path.inter.msa <- opt$path.inter.msa

if (is.null(path.features.msa) || is.null(path.inter.msa)) {
  stop("Error: both --path.features.msa and --path.inter.msa must be provided")
}

if (!dir.exists(path.features.msa)) stop('Features MSA directory does not exist')
if (!dir.exists(path.inter.msa)) stop('Internal MSA directory does not exist')

# ***************************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type.in, ".*")
files <- list.files(path = path.features.msa, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.in, "", files)
pref.combinations <- sub(".h5", "", pref.combinations)

if(length(pref.combinations) == 0) {
  stop('No files with the ref-based alignments are found')
}

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

# Rebuild the set of completed combinations from all worker logs (any core count)
done.set <- getDoneSet(path.log)
if(length(done.set) > 0){
  pokaz('Skip already done:', length(done.set), file=file.log.main, echo=echo.main)
}
assign('.worker.id', 1, envir = .GlobalEnv)   # sequential outer loop -> single core_1.log

for(s.comb in pref.combinations){

  # ---- Checkpoint: skip already completed combinations ----
  s.comb.id <- s.comb
  if(s.comb.id %in% done.set) next

  # One log file per worker (bounded number of files); also the checkpoint ledger
  file.log.loop = initLoopLog(path.log)

  pokaz('* Combination', s.comb, file=file.log.loop, echo=echo.loop)
  q.chr = strsplit(s.comb, '_')[[1]][1]
  
  file.comb <- file.path(path.features.msa, paste0(aln.type.in, s.comb, ".h5"))
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  accessions = intersect(accessions, accessions.specified)
  n.acc = length(accessions)
  
  # ---- Read Breaks ----
  file.breaks <- file.path(path.inter.msa,    paste0("breaks_", s.comb, ".rds"))
  if(!file.exists(file.breaks)) {
    pokaz('File', file.breaks, 'does not exist')
    stop()
  }
  breaks = readRDS(file.breaks)
  n.init = nrow(breaks)
  
  breaks$in.anal = (breaks$len.acc <= len.large) & (breaks$len.comb <= len.large)
  breaks.extra = breaks[!breaks$in.anal,]  # what should be aligned later
  
  breaks = breaks[breaks$in.anal,]
  breaks.init = breaks
  
  # ---- Merge coverages ----
  pokaz('Merge coverages..', file=file.log.loop, echo=echo.loop)
  breaks <- mergeOverlapsTolerance(breaks)
  
  if(sum(breaks$cnt) != nrow(breaks.init)) stop('Checkpoint3')
  
  # ---- Solve long ----
  pokaz('Solve long..', file=file.log.loop, echo=echo.loop)
  idx.rem.init = solveLong2(breaks, breaks.init, len.large)
  if(length(idx.rem.init) > 0){
    breaks.extra = rbind(breaks.extra, breaks.init[idx.rem.init,])
    breaks.init = breaks.init[-idx.rem.init,]  
    breaks = mergeOverlapsTolerance(breaks.init)
  }
  
  if((sum(breaks$cnt) + nrow(breaks.extra)) != n.init) stop('Checkout length')
  
  pokaz('Save extra breaks..', file=file.log.loop, echo=echo.loop)
  
  file.breaks.extra <- file.path(path.inter.msa, paste0("breaks_extra_", s.comb, ".rds"))
  saveRDS(breaks.extra, file.breaks.extra)
  
  ## ---- Get begin-end positions of gaps ----
  pokaz('Get begin-end positions of gaps..', file=file.log.loop, echo=echo.loop)
  pokaz('Get begin-end positions of gaps..')
  
  # n.breaks <- nrow(breaks)
  # n.acc    <- length(accessions)
  # v.beg <- matrix(0, nrow = n.breaks, ncol = n.acc)
  # v.end <- matrix(0, nrow = n.breaks, ncol = n.acc)
  # 
  # for (i in seq_along(accessions)) {
  #   acc <- accessions[i]
  #   
  #   pokaz(acc, file=file.log.loop, echo=echo.loop)
  #   
  #   x.acc = h5read(file.comb, paste0(gr.accs.e, acc))
  #   b.acc = h5read(file.comb, paste0(gr.blocks, acc))
  #   
  #   x.beg = fillPrev(x.acc)[breaks$idx.beg]
  #   x.end = fillNext(x.acc)[breaks$idx.end]
  #   
  #   idx.no.zero = (x.beg != 0) & (x.end != 0)
  #   idx.no.zero[idx.no.zero] = b.acc[abs(x.beg[idx.no.zero])] == b.acc[abs(x.end[idx.no.zero])]
  #   
  #   x.beg[!idx.no.zero] = 0
  #   x.end[!idx.no.zero] = 0
  #   
  #   v.beg[, i] <- x.beg
  #   v.end[, i] <- x.end
  #   
  # }
  # colnames(v.beg) = accessions
  # colnames(v.end) = accessions
  # 
  # -------
  n.breaks <- nrow(breaks)
  n.acc    <- length(accessions)
  
  idx.beg <- breaks$idx.beg
  idx.end <- breaks$idx.end
  
  cl <- parallel::makeCluster(num.cores)
  registerDoParallel(cl)
  
  res <- foreach(
    i = seq_along(accessions),
    .packages = c("rhdf5", "pannagram")
  ) %dopar% {
    
    acc <- accessions[i]
    
    x.acc <- h5read(file.comb, paste0(gr.accs.e, acc))
    b.acc <- h5read(file.comb, paste0(gr.blocks, acc))
    
    x.beg <- fillPrev(x.acc)[idx.beg]
    x.end <- fillNext(x.acc)[idx.end]
    
    idx.no.zero <- (x.beg != 0L) & (x.end != 0L)
    
    if (any(idx.no.zero)) {
      idx.no.zero[idx.no.zero] <-
        b.acc[abs(x.beg[idx.no.zero])] == b.acc[abs(x.end[idx.no.zero])]
    }
    
    x.beg[!idx.no.zero] <- 0L
    x.end[!idx.no.zero] <- 0L
    
    list(beg = x.beg, end = x.end)
  }
  
  stopCluster(cl)
  
  pokaz('Combine..', file=file.log.loop, echo=echo.loop)
  pokaz('Combine..')
  v.beg <- do.call(cbind, lapply(res, `[[`, "beg"))
  v.end <- do.call(cbind, lapply(res, `[[`, "end"))
  
  colnames(v.beg) <- accessions
  colnames(v.end) <- accessions
  
  # -------
  pokaz('Filter...')
  # Filter "extra" breaks. NOTE: this zeroing is load-bearing -- besides deferring
  # large insertions, it removes duplicate v.beg/v.end that arise when one
  # accession's big insertion spans several merged breaks (fillPrev/fillNext return
  # the same flanking positions), which would otherwise trip the duplicate check
  # below. A side effect is that at mixed-size overlapping loci it can leave a single
  # accession -> a phantom `single==1 & gap!=0` break (see comb_05 logging).
  for(acc in accessions){
    breaks.acc = breaks.extra[breaks.extra$acc == acc,]
    if(nrow(breaks.acc) == 0) next
    # Union of overlaps computed on the original column values (equivalent to the
    # progressive zeroing: a not-yet-removed position always keeps its original
    # value, so the first matching interval removes it -> same final set).
    vb = v.beg[,acc]; ve = v.end[,acc]
    idx.remove = logical(length(vb))
    for(irow in 1:nrow(breaks.acc)){
      idx.remove = idx.remove | ((vb <= breaks.acc$val.end[irow]) & (ve >= breaks.acc$val.beg[irow]))
    }
    v.beg[idx.remove,acc] = 0
    v.end[idx.remove,acc] = 0
  }

  # Check inversions
  if (any(sign(v.beg * v.end) < 0)) stop('Checkpoint4')
  
  # Check direction
  if (any(sign(v.end - v.beg) < 0)){
    # save(list = ls(), file = paste0("tmp_workspace_checkpoint5_", s.comb,".RData"))
    stop('Checkpoint5')
  } 
  
  # ---- Zero-positions mask ----
  pokaz("Zero-positions mask...")
  zero.mask = (v.end == 0) | (v.beg == 0)
  v.end[zero.mask] = 0
  v.beg[zero.mask] = 0
  
  # ---- Check lengths ----
  v.len = v.end - v.beg - 1
  v.len[zero.mask] = 0
  
  if (any(v.len < 0)) stop('Checkpoint6')
  
  zero.len.mask = (v.len == 0)
  v.end[zero.len.mask] = 0
  v.beg[zero.len.mask] = 0
  
  # ---- Checkups for duplicates ----
  pokaz("Checkups for duplicates...")
  for(icol in 1:ncol(v.len)){
    # Guard with anyDuplicated over non-zero entries (short-circuits); the full
    # dup value set is only recomputed to build the error message on failure.
    col.b = v.beg[,icol]; col.b = col.b[col.b != 0]
    if(anyDuplicated(col.b) > 0) {
      idx.dup = unique(v.beg[duplicated(v.beg[,icol]),icol])
      stop(paste('Duplicated in column', icol, 'in v.beg, amount:', length(idx.dup) - 1))  # WHY -1 ?!
    }
    col.e = v.end[,icol]; col.e = col.e[col.e != 0]
    if(anyDuplicated(col.e) > 0) {
      idx.dup = unique(v.end[duplicated(v.end[,icol]),icol])
      stop(paste('Duplicated in column', icol, 'in v.end, amount:', length(idx.dup) - 1))  # WHY -1 ?!
    }
  }
  
  idx.zero = which(rowSums(v.beg != 0) == 0)
  if(length(idx.zero) != 0){
    pokaz('Number of zero-breaks is', length(idx.zero), file=file.log.loop, echo=echo.loop)
    v.beg = v.beg[-idx.zero,,drop=FALSE]
    v.end = v.end[-idx.zero,,drop=FALSE]
    breaks = breaks[-idx.zero,]
    v.len = v.len[-idx.zero,,drop=FALSE]
  }
  
  # ---- Subdivide into categories ----
  pokaz("Subdivide into categories...")
  breaks$single = rowSums(v.len != 0)
  breaks$len.acc = rowMax(v.len)
  v.len[v.len == 0] <- NA
  breaks$len.mean = rowMeans(v.len, na.rm = TRUE)

  # ---- Defer residual single-with-gap breaks to extra ----
  # A break left with exactly one accession (single==1) but a non-zero reference gap
  # is a phantom singleton: overlapping long insertions were deferred to extra,
  # leaving one accession whose interior is occupied by the deferred ones. It is not a
  # true singleton (which has gap==0), so defer it to extra too -- recorded for the
  # later extra pass -- instead of leaving it unclassifiable in comb_05.
  gap.res = breaks$idx.end - breaks$idx.beg - 1
  idx.single.extra = which((breaks$single == 1) & (gap.res != 0))
  if(length(idx.single.extra) > 0){
    pokaz('Defer residual single-with-gap breaks to extra:', length(idx.single.extra),
          file=file.log.loop, echo=echo.loop)
    file.single.extra <- file.path(path.inter.msa, paste0("breaks_single_extra_", s.comb, ".rds"))
    saveRDS(breaks[idx.single.extra, , drop=FALSE], file.single.extra)
    breaks = breaks[-idx.single.extra, , drop=FALSE]
    v.beg  = v.beg[-idx.single.extra, , drop=FALSE]
    v.end  = v.end[-idx.single.extra, , drop=FALSE]
  }

  pokaz("Save...")
  all.local.objects <- c("breaks", "v.end", "v.beg", "accessions")
  file.ws <- file.path(path.inter.msa, paste0("breaks_ws_", s.comb, ".RData"))
  save(list = all.local.objects, file = file.ws)
  
  H5close()
  gc()

  # ---- Checkpoint marker: combination fully processed ----
  markDone(s.comb.id, file=file.log.loop, echo=echo.loop)
}

