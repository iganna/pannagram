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

option_list <- list(
  make_option("--path.mafft.in",  type = "character", default = NULL, help = "Path to directory where to combine fasta files for mafft runs"),
  make_option("--path.mafft.out", type = "character", default = NULL, help = "Path to directory where mafft results are"),
  make_option("--cores",          type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",       type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",      type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

n.flank = 30

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores = opt$cores

# Paths
if (!is.null(opt$path.mafft.in)) path.mafft.in <- opt$path.mafft.in
if (!is.null(opt$path.mafft.out)) path.mafft.out <- opt$path.mafft.out

# ***********************************************************************
# ---- Preparation ----

files.in <- list.files(path = path.mafft.in, pattern = "\\.fasta$", full.names = F)
files.out <- list.files(path = path.mafft.out, pattern = "\\.fasta$", full.names = F)
files.extra = setdiff(files.in, gsub('_aligned', '', files.out))

path.mafft.in.tmp = paste0(path.mafft.in, 'tmp/')

if (!file.exists(path.mafft.in.tmp)) {
  dir.create(path.mafft.in.tmp)
}

loop.function <- function(f.in, 
                          echo.loop=T){
  
  pokaz(f.in)
  # Log files
  file.log.loop = paste0(path.log, 'loop_file_', 
                         sub("\\.[^.]*$", "", basename(f.in)),
                         '.log')
  invisible(file.create(file.log.loop))

  # ---- Check log Done ----
  if(checkDone(file.log.loop)){
    return()
  }
  
  seqs = readFasta(paste0(path.mafft.in, f.in))
  seqs.clean = seq2clean(seqs,n.flank)
  
  
  path.work = paste0(path.mafft.in.tmp, sub('.fasta', '', basename(f.in)), '_')
  pokaz(path.work)
  res = refineAlignment(seqs.clean, path.work)
  
  alignments = res$aln
  
  alignment = alignments[[length(alignments)]]
  
  alignment.seq = mx2aln(alignment)
  
  
  file.out = paste0(path.mafft.out, sub('.fasta', '', basename(f.in)), "_aligned2.fasta")
  writeFasta(alignment.seq, file.out)
  
}


# ***********************************************************************
# ---- Loop  ----

if(num.cores == 1){
  for(f.in in files.extra){
    loop.function(f.in, echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(f.in = files.extra, 
                .packages=c('crayon'), 
                .verbose = F)  %dopar% { 
                  loop.function(f.in, echo.loop=echo.loop)
                }
  stopCluster(myCluster)
}

warnings()

pokaz('Done.', file=file.log.main, echo=echo.main)


