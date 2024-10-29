# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
  library(msa)
  library(muscle)  #BiocManager::install("muscle")
  library(Biostrings)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func.R", package = "pannagram"))
# source("synteny_funcs.R")

# pokazStage('Step 10. Prepare sequences for MAFFT')

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.cons",   type = "character", default = NULL, help = "Path to consensus directory"),
  make_option("--cores",       type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",    type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",   type = "character", default = NULL, help = "Level of log to be shown on the screen"),
  make_option("--max.len.gap", type = "integer",   default = NULL, help = "Max length of the gap")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);


#TODO: SHOULD BE PARAMATERS
len.short = 50
# len.large = 40000
n.flank = 30

# print(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

# ***********************************************************************
# ---- Values of parameters ----

# Max len gap
if (is.null(opt$max.len.gap)) {
  stop("Error: max.len.gap is NULL")
} else {
  len.large <- opt$max.len.gap
}

# Number of cores for parallel processing
# Number of cores for parallel processing
num.cores = opt$cores
if(is.null(num.cores)) stop('Whong number of cores: NULL')
pokaz('Number of cores', num.cores)
if(num.cores > 1){
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
}

# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop('Consensus folder doesn???t exist')

# Path to chromosomes
if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes

# Path to mafft input
if (!is.null(opt$path.mafft.in)) path.mafft.in <- opt$path.mafft.in

# **************************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type.comb, ".*")
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.comb, "", files)
pref.combinations <- sub(".h5", "", pref.combinations)

if(length(pref.combinations) == 0) {
  stop('No files with the ref-based alignments are found')
}

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

echo = F
for(s.comb in pref.combinations){
  
  if(echo) pokaz('* Combination', s.comb)
  
  # -------------------------------------
  
  file.ws = paste0(path.cons, 'small_ws_', s.comb, '.RData')
  load(file.ws)
  
  n.acc = length(accessions)
  
  
  path.tmp = '/Volumes/Samsung_T5/vienn/test/manuals/ecoli_out/intermediate/consensus/'
  
  # ---- Align Short sequences ----
  if(echo) pokaz('Align short seqs')
  
  # Checkup the number of sequences in alignments
  tmp = unlist(lapply(aln.seqs[idx.short], length))
  if(sum(tmp == 1) > 0) stop('Checkup-short1')
  
  # Core core for the short alignmgnets
  CODE_ALN_BATCH <- function(i.batch, echo=F){
    idx.use = idx.batches[[i.batch]]
    
    # file.in = paste0(path.tmp, 'batch_',i.batch, '.fasta')
    # file.out = paste0(path.tmp, 'batch_',i.batch, '_aligned.fasta')
    
    mx.list = vector("list", length = length(idx.use))
    
    for(i.aln in 1:length(idx.use)){
      pokaz(i.aln)
      idx.aln = idx.use[i.aln]
      seqs = aln.seqs[[idx.aln]]
      names(seqs) = aln.seqs.names[[idx.aln]]
      
      seqs <- DNAStringSet(seqs)
      # alignment <- msa(seqs, method = "ClustalOmega")
      # aln = as.character(alignment)
      
      alignment = muscle(seqs, quiet = T)
      aln = as.character(alignment)
      
      n.pos = nchar(aln[1])
      
      aln.text <- paste(names(aln), collapse = "\n")
      aln.info <- read.table(text = aln.text, sep = "|", stringsAsFactors = FALSE)

      mx.pos = matrix(0, nrow=n.pos, ncol=n.acc, dimnames = list(NULL, accessions))
      for(i.s in 1:nrow(aln.info)){
        acc = aln.info$V1[i.s]
        s = seq2nt(aln[i.s])
        mx.pos[s != '-',acc] = aln.info$V3[i.s]:aln.info$V4[i.s]
      }
      mx.list[[i.aln]] = mx.pos
      rm(aln.info)
      rm(aln)
      rm(mx.pos)
      rm(alignment)
      rm(seqs)
    }
    return(mx.list)
  }
  
  # save(list = ls(), file = "tmp_workspace_s.RData")
  # stop()
  
  # Two possible loops depending on the number of cores
  if(num.cores == 1){
    pokaz('No parallel computing: short sequences')
    # One core
    
    idx.batches <- list(idx.short)
    res.msa = CODE_ALN_BATCH(1)
    
  } else {
    # Many cores
    
    # Split the vector of shotp positions into n.core batches
    idx.batches <- split(idx.short, cut(seq_along(idx.short), num.cores, labels = FALSE))
    
    pokaz('Parallel computing: short sequences')
    res.msa <- foreach(i.batch = 1:num.cores,
                       .packages=c('Biostrings', 'crayon', 'msa', 'muscle'))  %dopar% {
                         return(CODE_ALN_BATCH(i.batch)) 
                       }
    
    res = list()
    for(i.batch in 1:num.cores){
      res[idx.batches[[i.batch]]] = res.msa[[i.batch]]
    }
    
    res.msa = res
  }
  
  saveRDS(list(aln = res.msa,
               ref.pos = data.frame(beg = breaks$idx.beg[idx.short],
                                    end = breaks$idx.end[idx.short]) ),
          paste0(path.cons, 'aln_short_',s.comb,'2.rds'), compress = F)
  
}


if(num.cores > 1){
  stopCluster(myCluster)
}

warnings()


# ***********************************************************************
# ---- Manual testing ----

if(F){

  library(rhdf5)
  path.cons = '/Volumes/Samsung_T5/vienn/test/symA_test_0/intermediate/consensus/'
  path.mafft.in = '/Volumes/Samsung_T5/vienn/test/symA_test_0/intermediate/mafft_in/'
  path.chromosomes = '/Volumes/Samsung_T5/vienn/test/symA_test_0/intermediate/chromosomes/'
  
}




