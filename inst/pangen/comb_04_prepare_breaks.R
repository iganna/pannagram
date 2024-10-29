# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
  library(muscle)  #BiocManager::install("muscle")
  # library(Biostrings)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.cons",        type = "character", default = NULL, help = "Path to consensus directory"),
  make_option("--path.chromosomes", type = "character", default = NULL, help = "Path to directory with chromosomes"),
  
  make_option("--max.len.gap",      type = "integer",   default = NULL, help = "Max length of the gap"),
  
  make_option("--cores",            type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",         type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",        type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Logging ----
source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----
source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

aln.type.in = aln.type.clean

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
if(is.null(num.cores)) stop('Whong number of cores: NULL')

pokaz('Number of cores', num.cores)
if(num.cores > 1){
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
}
# num.cores.max = 10
# num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))
# if(num.cores > 1){
#   myCluster <- makeCluster(num.cores, type = "PSOCK") 
#   registerDoParallel(myCluster)
# }

# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop('Consensus folder doesn???t exist')

# Path to chromosomes
if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes

# **************************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type.in, ".*")
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.in, "", files)
pref.combinations <- sub(".h5", "", pref.combinations)

if(length(pref.combinations) == 0) {
  stop('No files with the ref-based alignments are found')
}

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

echo = T
for(s.comb in pref.combinations[3]){
  
  if(echo) pokaz('* Combination', s.comb)
  q.chr = strsplit(s.comb, '_')[[1]][1]
  
  file.comb = paste0(path.cons, aln.type.in, s.comb,'.h5')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)
  
  # ---- Read Breaks  ----
  file.breaks = paste0(path.cons, 'breaks_', s.comb,'.rds')
  breaks = readRDS(file.breaks)
  n.init = nrow(breaks)
  
  breaks$in.anal = (breaks$len.acc <= len.large) & (breaks$len.comb <= len.large)
  breaks.extra = breaks[!breaks$in.anal,]  # what should be aligned later
  
  breaks = breaks[breaks$in.anal,]
  breaks.init = breaks
  
  # ---- Merge coverages ----
  breaks <- mergeOverlaps(breaks)
  
  if(sum(breaks$cnt) != nrow(breaks.init)) stop('Checkpoint3')
  
  # ---- Solve long ----
  idx.rem.init = solveLong(breaks, breaks.init, len.large)
  if(length(idx.rem.init) > 0){
    breaks.extra = rbind(breaks.extra, breaks.init[idx.rem.init,])
    breaks.init = breaks.init[-idx.rem.init,]  
    breaks = mergeOverlaps(breaks.init)
  }
  
  if((sum(breaks$cnt) + nrow(breaks.extra)) != n.init) stop('Checkout length')
  
  file.breaks.extra = paste0(path.cons, 'breaks_extra_', s.comb,'.rds')
  saveRDS(breaks.extra, file.breaks.extra)
  
  ## ---- Get begin-end positions of gaps ----
  v.beg = c()
  v.end = c()
  for(acc in accessions){
    pokaz(acc)
    
    x.acc = h5read(file.comb, paste0(gr.accs.e, acc))
    b.acc = h5read(file.comb, paste0(gr.blocks, acc))
    
    x.beg = fillPrev(x.acc)[breaks$idx.beg]
    x.end = fillNext(x.acc)[breaks$idx.end]
    
    idx.no.zero = (x.beg != 0) & (x.end != 0)
    idx.no.zero[idx.no.zero] = b.acc[abs(x.beg[idx.no.zero])] == b.acc[abs(x.end[idx.no.zero])]
    
    x.beg[!idx.no.zero] = 0
    x.end[!idx.no.zero] = 0
    
    v.beg = cbind(v.beg, x.beg)
    v.end = cbind(v.end, x.end)
    
  }
  colnames(v.beg) = accessions
  colnames(v.end) = accessions
  
  save(list = ls(), file = "tmp_workspace.RData")
  stop('Enough..')
  
  
  # Filter "extra" breaks
  for(acc in accessions){
    breaks.acc = breaks.extra[breaks.extra$acc == acc,]
    for(irow in 1:nrow(breaks.acc)){
      idx.remove = (v.beg[,acc] <= breaks.acc$val.end[irow]) & (v.end[,acc] >= breaks.acc$val.beg[irow])
      
      v.beg[idx.remove,acc] = 0
      v.end[idx.remove,acc] = 0
    }
  }
  
  # Check inversions
  if (any(sign(v.beg * v.end) < 0)) stop('Checkpoint4')
  # To test: which(rowSums(sign(v.beg * v.end) < 0) > 0)
  
  # Check direction
  if (any(sign(v.end - v.beg) < 0)) stop('Checkpoint5')
  
  # ---- Zero-positions mask ----
  
  zero.mask = (v.end == 0) | (v.beg == 0)
  v.end[zero.mask] = 0
  v.beg[zero.mask] = 0
  
  # ---- Check lengths ----
  v.len = v.end - v.beg - 1
  v.len[zero.mask] = 0
  
  if (any(v.len < 0)) stop('Checkpoint6')
  v.len[zero.mask] = 0
  
  zero.len.mask = (v.len == 0)
  v.end[zero.len.mask] = 0
  v.beg[zero.len.mask] = 0
  
  # ---- Checkups for duplicates ----
  for(icol in 1:ncol(v.len)){
    idx.dup = unique(v.beg[duplicated(v.beg[,icol]),icol])
    if(length(setdiff(idx.dup, 0)) != 0) {
      stop(paste('Duplicated in column', icol, 'in v.beg, amount:', length(idx.dup) - 1))
    }
    # v.len[v.beg[,icol] %in% idx.dup, icol] = 0
    
    idx.dup = unique(v.end[duplicated(v.end[,icol]),icol])
    if(length(setdiff(idx.dup, 0)) != 0) {
      stop(paste('Duplicated in column', icol, 'in v.end, amount:', length(idx.dup) - 1))
    }
    # v.len[v.end[,icol] %in% idx.dup, icol] = 0
  }
  
  # ---- Subdivide into categories ----
  
  breaks$single = rowSums(v.len != 0)
  breaks$len.acc = rowMax(v.len)
  v.len[v.len == 0] <- NA
  breaks$len.mean = rowMeans(v.len)
  breaks$len.mean <- rowMeans(v.len, na.rm = TRUE)
  
  all.local.objects <- c("breaks", "v.end", "v.beg", "accessions")
  file.ws = paste0(path.cons, 'breaks_ws_', s.comb, '.RData')
  save(list = all.local.objects, file = file.ws)
  
  H5close()
  gc()
}

if(num.cores > 1){
  stopCluster(myCluster)
}

warnings()


# ***********************************************************************
# ---- Manual testing ----




