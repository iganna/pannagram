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
  make_option("--path.mafft.out", type = "character", default = NULL, help = "Path to directory where mafft results are"),
  make_option("--cores",          type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",       type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",      type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

n.flank = 30
sim.cutoff = 0.2

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores = opt$cores

# Paths
if (!is.null(opt$path.mafft.out)) path.mafft.out <- opt$path.mafft.out

# ***********************************************************************
# ---- Preparation ----

files.out <- list.files(path = path.mafft.out, pattern = "_aligned\\.fasta$", full.names = F)

# ***********************************************************************
# ---- MAIN program body ----

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
  
  seqs = readFasta(paste0(path.mafft.out, f.in))
  mx = aln2mx(seqs)
  
  len.aln = ncol(mx)
  
  mx = toupper(mx)
  pos.profile = mx2profile(mx)
  seq.cons = mx2cons(mx)
  pos.variation = (colSums(pos.profile == 0) != 3) * 1
  
  # Define blocks, were the alignment non well
  blocks.all = c()
  for(i.seq in 1:nrow(mx)){
    s = mx[i.seq,]
    blocks = findOnes((s != '-') *1)
    blocks$len = blocks$end - blocks$beg + 1
    if(nrow(blocks) == 0) next
    
    # Estimate diversity within each block
    blocks$pi = 0
    blocks$acc = rownames(mx)[i.seq]
    for(irow in 1:nrow(blocks)){
      idx.block = blocks$beg[irow]:blocks$end[irow]
      blocks$pi[irow] = sum(pos.variation[idx.block]) / blocks$len[irow]
    }
    blocks.all = rbind(blocks.all, blocks)
  }
  
  
  blocks.all = blocks.all[blocks.all$pi > sim.cutoff,, drop=F]
  
  
  if(nrow(blocks.all) > 0){
    pokaz('Bad alignment', f.in)
  }
  
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


