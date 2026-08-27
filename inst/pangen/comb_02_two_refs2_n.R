#' How the output files look like:
#'     group           name       otype   dclass      dim
#'         /           accs   H5I_GROUP                  
#'     /accs              0 H5I_DATASET    FLOAT 28940631
#'     /accs          10002 H5I_DATASET    FLOAT 28940631
#'     /accs          10015 H5I_DATASET    FLOAT 28940631


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
  make_option("--ref0",       type = "character", default = NULL, help = "Reference file 1"),
  make_option("--ref1",       type = "character", default = NULL, help = "Reference file 2"),
  make_option("--cores",      type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",   type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",  type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Logging ----
source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----
source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files
source(system.file("utils/interval_func.R", package = "pannagram")) # interval codec + hdf5 layer

aln.type.in = paste0(aln.type.ref, '_')
aln.type.out = paste0(aln.type.comb, '_')

# ***********************************************************************
# ---- Accessions ----

file.acc <- ifelse(!is.null(opt$accessions), opt$accessions, stop("File with accessions are not specified"))
accessions.specified <- as.character(read.table(file.acc, stringsAsFactors = FALSE)[, 1])

# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores <- opt$cores

# Path with the consensus output
if (!is.null(opt$path.features.msa)) path.features.msa <- opt$path.features.msa
if(!dir.exists(path.features.msa)) stop('path_features_msa directory does not exist')

if (!is.null(opt$path.inter.msa)) path.inter.msa <- opt$path.inter.msa
if(!dir.exists(path.inter.msa)) stop('path_inter_msa folder does not exist')

# Reference genomes
ref0 <- if (is.null(opt$ref0)) stop("opt$ref0 is NULL") else opt$ref0
ref1 <- if (is.null(opt$ref1)) stop("opt$ref1 is NULL") else opt$ref1

pokaz('References:', ref0, ref1, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

pattern <- paste0("^", aln.type.in, "[0-9]+_[0-9]+_.*\\.h5$")
combo_files <- list.files(path = path.features.msa, pattern = pattern, full.names = F)
combo_files <- combo_files[grep(paste(ref1, ref0, sep = "|"), combo_files)]

extract_xy <- function(filename) {
  parts <- strsplit(filename, "_")[[1]]
  x <- parts[2]
  y <- parts[3]
  
  if(x != y) return(NA)  # CHANGE IN FUTURE
  return(paste(x, y, sep = '_'))
}

pref.combinations <- unique(sapply(combo_files, extract_xy))

if (any(is.na(pref.combinations))) pokazAttention('Alignments will only be processed in chromosome-to-chromosome mode, i.e. chromosome 1 with chromosome 1, chromosome 2 with chromosome 2, and so on.')
pref.combinations <- pref.combinations[!is.na(pref.combinations)]  # THIS IS THE FITURE

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
  file.comb0 <- file.path(path.features.msa, paste0(aln.type.in, s.comb, "_", ref0, ".h5"))
  file.comb1 <- file.path(path.features.msa, paste0(aln.type.in, s.comb, "_", ref1, ".h5"))

  # Combined file. If it exists, then use it for the growing correspondence
  file.res   <- file.path(path.inter.msa, paste0(aln.type.out, s.comb, ".h5"))
  if(file.exists(file.res)){
    if(v.idx.trust %in% h5ls(file.res)$name) {
      file.comb0 = file.res
    } else {
      # delete file
      file.remove(file.res)
    }
  }

  pokaz("Files", file.comb0, file.comb1, file=file.log.loop, echo=echo.loop)

  if(!file.exists(file.res)){
    h5createFile(file.res)
    h5createGroup(file.res, gr.accs.e)
  }

  pokaz(file.comb0, file=file.log.loop, echo=echo.loop)
  pokaz(file.comb1, file=file.log.loop, echo=echo.loop)

  # Get the corresponsing function between two references

  # The pangenome length now comes from the dataset attribute, not from its dim
  v.ref0  <- h5VecRead(file.comb0, ref0)
  len.aln <- length(v.ref0)
  h5PanLenSet(file.res, len.aln)

  f01 <- cbind(v.ref0, h5VecRead(file.comb0, ref1, len.aln))
  rmSafe(v.ref0)
  f01 = f01[f01[,1] != 0,,drop=F]
  f01 = f01[f01[,2] != 0,,drop=F]

  # Get accessions to combine
  groups0 = h5ls(file.comb0)
  groups1 = h5ls(file.comb1)

  accessions = intersect(groups0$name[groups0$group == gr.accs.b],
                         groups1$name[groups1$group == gr.accs.b])  # full name of accessions
  accessions = intersect(accessions, accessions.specified)

  pokaz('Number of accessions', length(accessions), file=file.log.loop, echo=echo.loop)
  for(acc in accessions){

    # ---- Checkpoint: skip accessions already written for this combination ----
    acc.id <- paste0(s.comb.id, '_', acc)
    if(acc.id %in% done.set){
      pokaz('Accession already done, skip:', acc, file=file.log.loop, echo=echo.loop)
      next
    }

    # --- --- --- --- --- --- --- --- --- --- ---
    pokaz('Accession', acc, file=file.log.loop, echo=echo.loop)
    s = paste0('/',gr.accs.e, acc)

    # Data from the main reference
    v0 = h5VecRead(file.comb0, acc, len.aln)
    v1 = h5VecRead(file.comb1, acc)
    v.final = v0
    v.final[f01[,1]] = 0

    v0 = v0[f01[,1]]
    v1 =  v1[abs(f01[,2])] * sign(f01[,2])
    v0[v1 != v0] = 0
    v.final[f01[,1]] = v0

    nz = v.final[v.final != 0]
    if(anyDuplicated(nz) > 0){
      dup.value = unique(nz[duplicated(nz)])
      v.final[v.final %in% dup.value] <- 0
      pokaz('Number of duplicated', length(dup.value), file=file.log.loop, echo=echo.loop)
    }

    pokaz('Length of saved vector', length(v.final), file=file.log.loop, echo=echo.loop)

    # Write into file as an interval table (h5VecWrite is idempotent)
    h5VecWrite(file.res, acc, v.final)

    # ---- Checkpoint marker: accession fully written ----
    markDone(acc.id, file=file.log.loop, echo=echo.loop)

  }

  # ---- Update Idx trust (accumulator: run exactly once, guarded by its own marker) ----
  # Unlike per-accession writes, this increments a counter and is NOT idempotent, so a
  # dedicated checkpoint marker prevents double-counting if a resume re-enters here.
  trust.id <- paste0(s.comb.id, '_trust')
  if(!(trust.id %in% done.set)){
    suppressMessages({

      if(v.idx.trust %in% h5ls(file.res)$name) {
        idx.trust = h5read(file.res, v.idx.trust)
      } else {
        idx.trust = rep(0, len.aln)
      }

      idx.trust[f01[,1]] = idx.trust[f01[,1]] + 1
      h5write(idx.trust, file.res, v.idx.trust)

    })
    markDone(trust.id, file=file.log.loop, echo=echo.loop)
  }

  H5close()

  # ---- Checkpoint marker: combination fully processed ----
  markDone(s.comb.id, file=file.log.loop, echo=echo.loop)

  gc()

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

# stop('TEST')

# ***********************************************************************
# ---- Manual testing ----



