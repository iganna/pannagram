# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
  # library(muscle)  #BiocManager::install("muscle")
  library(pannagram)
  library(igraph)
  # library(Biostrings)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func.R", package = "pannagram"))
source(system.file("pangen/comb_func_extra2.R", package = "pannagram"))
source(system.file("pangen/synteny_func.R", package = "pannagram"))
# source("synteny_funcs.R")

# pokazStage('Step 10. Prepare sequences for MAFFT')

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.cons"),       type = "character", default = NULL, help = "path to consensus directory"),
  make_option(c("--path.extra"),      type = "character", default = NULL, help = "path to directory, where to combine fasta files for mafft runs"),
  make_option(c("-c", "--cores"),     type = "integer",   default = 1,    help = "number of cores to use for parallel processing"),
  make_option(c("--path.log"),        type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"),       type = "character", default = NULL, help = "Level of log to be shown on the screen"),
  make_option(c("--aln.type.in"),      type = "character", default = NULL, help = "Alignment type of the input file"),
  make_option(c("--aln.type.out"),     type = "character", default = NULL, help = "Alignment type of the output file")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

# ***********************************************************************


aln.type.in <- opt$aln.type.in
if(is.null(aln.type.in)){
  aln.type.in = aln.type.add  
}

aln.type.out <- opt$aln.type.out
if(is.null(aln.type.out)){
  aln.type.out = aln.type.extra1  
}

pokazAttention("Prefs of in and out files:", aln.type.in, aln.type.out)

# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores = opt$cores
if(is.null(num.cores)) stop('Whong number of cores: NULL')

# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop('Consensus folder does not exist')

# Path to chromosomes
if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes

# Path to mafft input
if (!is.null(opt$path.extra)) path.extra <- opt$path.extra

path.work = paste0(path.extra, 'tmp/')
if (!dir.exists(path.work)) {
  dir.create(path.work)
}

# **************************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type.in, "\\d+_\\d+\\.h5$")
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
for(s.comb in pref.combinations){
  
  if(echo) pokaz('* Combination', s.comb)
  q.chr = strsplit(s.comb, '_')[[1]][1]
  
  file.comb = paste0(path.cons, aln.type.in, s.comb,'.h5')
  file.out = paste0(path.cons, aln.type.out, s.comb,'.h5')
  
  # Create the output file
  suppressMessages({
    h5createFile(file.out)
    h5createGroup(file.out, gr.blocks)
    h5createGroup(file.out, gr.accs.e)
  })
  
  # Accessions
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  
  len.comb = as.numeric(unique(groups$dim[groups$group == gr.accs.b]))
  
  file.breaks.info = paste0(path.extra, "breaks_info_",s.comb,".RData")
  load(file.breaks.info)
  
  if(sum(breaks$idx.beg == breaks$idx.end) > 0) stop('Breaks are wrong')
  
  # Define Lengths
  idx.gain = rep(1, len.comb)
  
  breaks$len = breaks$idx.end - breaks$idx.beg - 1  # MINUS! because these are not posisiotns, but positions around
  breaks$len.new = 0
  idx.skip = c()
  for (i.b in 1:nrow(breaks)) {
    file.br.len = paste0(path.extra, breaks$id.s[i.b], '_len.RData')
    file.br.out = paste0(path.extra, breaks$id.s[i.b], '_out.RData')
    if(file.exists(file.br.len) & file.exists(file.br.out)){
      load(file.br.len)
      breaks$len.new[i.b] = len.aln
    } else {
      idx.skip = c(idx.skip, i.b)
    }
  }
  
  if(length(idx.skip) > 0){
    pokazAttention('Number of breaks to skip is', length(idx.skip))
    breaks = breaks[-idx.skip,,drop=F]
    pokazAttention('Number of breaks remained', nrow(breaks))
  }
  
  breaks$len.plus = breaks$len.new - breaks$len
  idx.minus = which(breaks$len.plus < 0)
  if(length(idx.minus) > 0){
    stop('Reduce of the alignment')
  }
  
  idx.gain[breaks$idx.beg+1] = idx.gain[breaks$idx.beg+1] + breaks$len.plus
  new.idx = cumsum(idx.gain)
  len.new = max(new.idx)
  
  if(len.new != sum(length(idx.gain) + sum(breaks$len.plus))) stop('Something if wrong with positions')
  breaks$new.beg = new.idx[breaks$idx.beg] + 1
  breaks$new.end = new.idx[breaks$idx.end] - 1 
  
  len.to.check = breaks$new.end - breaks$new.beg + 1
  
  if(sum(len.to.check != breaks$len.new) > 0) stop('Calculation of positions is wrong')
  for(i.b in 1:nrow(breaks)){
    pos = breaks$idx.beg[i.b]:breaks$idx.end[i.b]
    pos = pos[-c(1, length(pos))]
    if(length(pos) > 0){
      new.idx[pos] = 0
    }
  }
  
  if((sum(new.idx != 0) + sum(breaks$len.new)) != len.new) stop('Somwthing wrong with lengths 2')
  
  pos.transfer = cbind(1:length(new.idx), new.idx)
  pos.transfer = pos.transfer[pos.transfer[,2] != 0,, drop = F]
  
  for(acc in accessions){
    pokaz('Accessions', acc)
    
    s.acc = paste0(gr.accs.e, acc)
    v = h5read(file.comb, s.acc)
    
    v.new = rep(0, len.new)
    v.new[pos.transfer[,2]] = v[pos.transfer[,1]]
    
    if (num.cores == 1) {
      for (i.b in 1:nrow(breaks)) {
        pokaz(i.b)
        file.br.out <- paste0(path.extra, breaks$id.s[i.b], '_out.RData')
        load(file.br.out) # idx.new and msa.new
        v.new[breaks$new.beg[i.b]:breaks$new.end[i.b]] <- idx.new[acc, ]
        
        rm(idx.new)
        rm(msa.new)
      }
    } else {
      # Number of cores
      myCluster <- makeCluster(num.cores, type = "PSOCK")
      registerDoParallel(myCluster)
      
      # Loop
      tmp <- foreach(i.b = 1:nrow(breaks), 
                     .packages = c('utils')) %dopar% {
                       file.br.out <- paste0(path.extra, breaks$id.s[i.b], '_out.RData')
                       load(file.br.out) # idx.new and msa.new
                       list(start = breaks$new.beg[i.b], 
                            end = breaks$new.end[i.b], 
                            values = idx.new[acc, ])
                     }
      
      # Save results
      for (res in tmp) {
        v.new[res$start:res$end] <- res$values
      }
      
      stopCluster(myCluster)
    }
    
    gc()
    
    if((sum(v.new != 0) + 1) != length(unique(v.new))){
      save(list = ls(), file ="tmp_workspace_dup.RData")
      stop("Duplicates are found")
    } 
    
    pokaz('Save new...')
    suppressMessages({
      h5write(v.new, file.out, s.acc)
    })
    
  }
  
}