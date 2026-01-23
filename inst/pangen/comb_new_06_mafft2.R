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
                          echo.loop=T){
  pokaz(f.in)
  # Log files
  file.log.loop = paste0(path.log, 'loop_file_', 
                         sub("\\.[^.]*$", "", basename(f.in)),
                         '.log')
  if(!file.exists(file.log.loop)){
    invisible(file.create(file.log.loop))
  }
  
  # ---- Check log Done ----
  if(checkDone(file.log.loop)){
    return()
  }
  
  seqs = readFasta(paste0(path.mafft.in, f.in))
  seqs.clean = seq2clean(seqs, n.flank)
  
  seqs.clean = seqs.clean[nchar(seqs.clean) > 7]
  if(length(seqs.clean) < 2) {
    pokazAttention('Not enough sequences to align', file=file.log.loop, echo=echo.loop)
    pokaz('Done.', file=file.log.loop, echo=echo.loop)
    return()
  }
  
  # Proportion of non-N nucleotides
  n.n = c()
  for(s in seqs.clean){
    s.tmp = seq2nt(s)
    n.n = c(n.n, sum((s.tmp != 'N') & (s.tmp != 'n')) / length(s.tmp))
  }
  
  if(min(n.n) < 0.5){
    pokazAttention('Too many N to align', file=file.log.loop, echo=echo.loop)
    pokaz('Done.', file=file.log.loop, echo=echo.loop)
    return()
  }
  
  path.work = paste0(path.mafft.in.tmp, sub('\\.fasta', '', basename(f.in)), '_')
  pokaz(path.work)
  res = refineAlignment(seqs.clean, path.work)
  
  alignments = res$aln
  
  alignment = alignments[[length(alignments)]]
  
  alignment.seq = mx2aln(alignment)
  
  file.out = paste0(path.mafft.out, sub('\\.fasta', '', basename(f.in)), "_aligned.fasta")
  writeFasta(alignment.seq, file.out)
  
  pokaz('Done.', file=file.log.loop, echo=echo.loop)
  
}


# ***********************************************************************
# ---- Loop  ----

if(num.cores == 1){
  for(f.in in files.extra){
    loop.function(f.in,
                  echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(f.in = files.extra, 
                .packages=c('crayon'), 
                .verbose = F)  %dopar% { 
                  loop.function(f.in,
                                echo.loop=echo.loop)
                }
  stopCluster(myCluster)
}

warnings()

pokaz('Done.', file=file.log.main, echo=echo.main)

