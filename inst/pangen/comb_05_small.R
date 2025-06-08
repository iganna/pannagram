# Get positional data for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
  library(msa)
  library(muscle)
  library(Biostrings)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.features.msa", type = "character", default = NULL, help = "Path to msa directory (features)"),
  make_option("--path.inter.msa",    type = "character", default = NULL, help = "Path to msa directory (internal)"),
  
  make_option("--cores",             type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",          type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",         type = "character", default = NULL, help = "Level of log to be shown on the screen"),
  make_option("--max.len.gap",       type = "integer",   default = NULL, help = "Max length of the gap")
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
num.cores = opt$cores
if(is.null(num.cores)) stop('Wrong number of cores: NULL')
pokaz('Number of cores', num.cores, file=file.log.main, echo=echo.main)

# Path with the MSA output (features)
path.features.msa <- opt$path.features.msa
path.inter.msa <- opt$path.inter.msa

if (is.null(path.features.msa) || is.null(path.inter.msa)) {
  stop("Error: both --path.features.msa and --path.inter.msa must be provided")
}

if (!dir.exists(path.features.msa)) stop('Features MSA directory doesn???t exist')
if (!dir.exists(path.inter.msa)) stop('Internal MSA directory doesn???t exist')

# **************************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type.comb, ".*")
files <- list.files(path = path.features.msa, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.comb, "", files)
pref.combinations <- sub(".h5", "", pref.combinations)

if(length(pref.combinations) == 0) {
  stop('No files with the ref-based alignments are found')
}
pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

for(s.comb in pref.combinations){
  # Log files
  file.log.loop = paste0(path.log, 'loop_file_', s.comb, '.log')
  if(!file.exists(file.log.loop)) invisible(file.create(file.log.loop)) 
  if(checkDone(file.log.loop)) next   # Check log Done
  
  initial.vars <- ls()
  
  pokaz('* Combination', s.comb, file=file.log.loop, echo=echo.loop)
  
  file.ws = paste0(path.inter.msa, 'small_ws_', s.comb, '.RData')
  load(file.ws)
  
  if(length(idx.short) == 0) next
  
  n.acc = length(accessions)
    
  # ---- Align Short sequences ----
  pokaz('Align short seqs', file=file.log.loop, echo=echo.loop)
  
  # Checkup the number of sequences in alignments
  tmp = unlist(lapply(aln.seqs[idx.short], length))
  if(sum(tmp == 1) > 0) stop('Checkup-short1')
  
  # Core core for the short alignmgnets
  CODE_ALN_BATCH <- function(i.batch, echo=F){
    idx.use = idx.batches[[i.batch]]
    
    mx.list = vector("list", length = length(idx.use))
    
    for(i.aln in 1:length(idx.use)){
      idx.aln = idx.use[i.aln]
      seqs = aln.seqs[[idx.aln]]
      names(seqs) = aln.seqs.names[[idx.aln]]
      
      seqs <- DNAStringSet(seqs)
      
      alignment = muscle(seqs, quiet = T)
      aln = as.character(alignment)
      
      save(list = ls(), file = "tmp_workspace_s.RData")
      
      n.pos = nchar(aln[1])
      
      aln.text <- paste(names(aln), collapse = "\n")
      aln.info <- read.table(text = aln.text, sep = "|", stringsAsFactors = FALSE)

      mx.pos = matrix(0, nrow=n.pos, ncol=n.acc, dimnames = list(NULL, accessions))
      for(i.s in 1:nrow(aln.info)){
        acc = as.character(aln.info$V1[i.s])
        s = seq2nt(aln[i.s])
        mx.pos[s != '-', acc] = aln.info$V3[i.s]:aln.info$V4[i.s]
      }
      mx.list[[i.aln]] = mx.pos
      rm(aln.info)
      rm(aln)
      rm(mx.pos)
      rm(alignment)
      rm(seqs)
      gc()
    }
    
    gc()
    
    return(mx.list)
  }
  
  
  # Two possible loops depending on the number of cores
  if(num.cores == 1){
    pokaz('No parallel computing: short sequences', file=file.log.loop, echo=echo.loop)
    # One core
    
    idx.batches <- list(idx.short)  # do not remove
    res.msa = CODE_ALN_BATCH(1)
    
  } else {
    # Many cores
    myCluster <- makeCluster(num.cores, type = "PSOCK") 
    registerDoParallel(myCluster) 
    
    # Split the vector of shotp positions into n.core batches
    idx.batches <- split(idx.short, cut(seq_along(idx.short), num.cores, labels = FALSE))
    num.cores = min(num.cores, length(idx.batches))
    
    pokaz('Parallel computing: short sequences', file=file.log.loop, echo=echo.loop)
    res.msa <- foreach(i.batch = 1:num.cores,
                       .packages=c('Biostrings', 'crayon', 'msa', 'muscle'))  %dopar% {
                         return(CODE_ALN_BATCH(i.batch)) 
                       }
    
    res = list()
    for(i.batch in 1:num.cores){
      res[idx.batches[[i.batch]]] = res.msa[[i.batch]]
    }
    
    res.msa = res.msa[idx.short]
    
    # Stop clusters
    stopCluster(myCluster)
    
  }
  
  ref.pos = data.frame(beg = breaks$idx.beg[idx.short],
                       end = breaks$idx.end[idx.short]) 
  
  if(length(res.msa) != nrow(ref.pos)){
    pokaz(length(res.msa), nrow(ref.pos), file=file.log.loop, echo=echo.loop)
    file.ws = "tmp_workspace_x.RData"
    all.local.objects <- ls()
    save(list = all.local.objects, file = file.ws)
    stop('Lengths dont match')
  }
    
  saveRDS(list(aln = res.msa,
               ref.pos = data.frame(beg = breaks$idx.beg[idx.short],
                                    end = breaks$idx.end[idx.short]) ),
          paste0(path.inter.msa, 'aln_short_',s.comb,'.rds'), compress = F)
  
  
  # Done
  pokaz('Done.', file=file.log.loop, echo=echo.loop)
  
  # Cleanup variables
  final.vars <- ls()
  new.vars <- setdiff(final.vars, initial.vars)
  pokaz('Variables to remove', new.vars, file=file.log.loop, echo=echo.loop)
  rm(list = new.vars)
  gc()
}
warnings()