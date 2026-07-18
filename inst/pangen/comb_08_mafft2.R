# Combine all alignments together into the final one

suppressMessages({
library('foreach')
library(doParallel)
library("optparse")
source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func_mafft_refine2.R", package = "pannagram"))
source(system.file("pangen/synteny_func.R", package = "pannagram"))
})

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.mafft.in"), type="character", default=NULL, 
              help="path to directory, where to combine fasta files for mafft runs", metavar="character"),
  make_option(c("--path.mafft.out"), type="character", default=NULL, 
              help="path to directory, where to mafft results are", metavar="character"),
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

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging


# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing
# num.cores.max = 10
# num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))
num.cores = opt$cores

if (!is.null(opt$path.mafft.in)) path.mafft.in <- opt$path.mafft.in
if (!is.null(opt$path.mafft.out)) path.mafft.out <- opt$path.mafft.out

# ***********************************************************************
# ---- Preparation ----

n.flank = 30

files.extra <- list.files(path = path.mafft.in, pattern = "\\.fasta$", full.names = F)
files.extra = files.extra[!grepl('_aligned', files.extra)]

if(length(files.extra) == 0){
  pokaz('Number of files for extra alignment')
  quit(status = 0)
} else {
  pokaz('Number of files to align', length(files.extra))
}

path.mafft.in.tmp = paste0(path.mafft.in, 'tmp/')

if (!file.exists(path.mafft.in.tmp)) {
  dir.create(path.mafft.in.tmp)
}

# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(f.in,
                          done.set = character(0),
                          echo.loop=T){
  pokaz(f.in)

  # ---- Checkpoint: skip already completed items ----
  item.id <- sub("\\.[^.]*$", "", basename(f.in))
  if(item.id %in% done.set){
    return(NULL)
  }

  # One log file per worker (bounded number of files); also the checkpoint ledger
  file.log.loop = initLoopLog(path.log)
  
  seqs = readFasta(paste0(path.mafft.in, f.in))
  seqs.clean = seq2clean(seqs, n.flank)
  
  seqs.clean = seqs.clean[nchar(seqs.clean) > 7]
  if(length(seqs.clean) < 2) {
    pokazAttention('Not enough sequences to align', file=file.log.loop, echo=echo.loop)
    markDone(item.id, file=file.log.loop, echo=echo.loop)
    return()
  }
  
  # Proportion of non-N nucleotides: drop N-heavy sequences INDIVIDUALLY
  # (a single N-rich sequence must not discard the whole locus).
  n.n = sapply(seqs.clean, function(s){
    s.tmp = seq2nt(s)
    sum((s.tmp != 'N') & (s.tmp != 'n')) / length(s.tmp)
  })
  seqs.clean = seqs.clean[n.n >= 0.5]

  if(length(seqs.clean) < 2){
    pokazAttention('Too many N to align', file=file.log.loop, echo=echo.loop)
    markDone(item.id, file=file.log.loop, echo=echo.loop)
    return()
  }
  
  path.work = paste0(path.mafft.in.tmp, sub('\\.fasta', '', basename(f.in)), '_')
  pokaz(path.work)
  res = refineAlignmentRouted(seqs.clean, path.work)
  
  alignments = res$aln
  
  alignment = alignments[[length(alignments)]]
  
  alignment.seq = mx2aln(alignment)
  
  file.out = paste0(path.mafft.out, sub('\\.fasta', '', basename(f.in)), "_aligned.fasta")
  writeFasta(alignment.seq, file.out)

  # ---- Checkpoint marker: item fully processed ----
  markDone(item.id, file=file.log.loop, echo=echo.loop)

}


# ***********************************************************************
# ---- Loop  ----

# Rebuild the set of completed items from all worker logs (any core count)
done.set <- getDoneSet(path.log)
if(length(done.set) > 0){
  pokaz('Skip already done:', length(done.set), file=file.log.main, echo=echo.main)
}

if(num.cores == 1){
  assign('.worker.id', 1, envir = .GlobalEnv)   # stable single file core_1.log
  for(f.in in files.extra){
    loop.function(f.in,
                  done.set = done.set,
                  echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK")
  registerDoParallel(myCluster)

  # Assign a stable worker id (1..N) to each worker -> bounded, reused file names
  parallel::clusterApply(myCluster, seq_len(num.cores),
                         function(i) assign('.worker.id', i, envir = .GlobalEnv))

  tmp = foreach(f.in = files.extra,
                .packages=c('crayon'),
                .verbose = F)  %dopar% {
                  loop.function(f.in,
                                done.set = done.set,
                                echo.loop=echo.loop)
                }
  stopCluster(myCluster)
}


pokaz('Done.', file=file.log.main, echo=echo.main)

