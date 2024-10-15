# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
})

source(system.file("utils/utils.R", package = "pannagram"))

# Define blocks in the alignemnt

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.cons", type = "character", default = NULL, help = "Path to consensus directory"),
  make_option("--ref",       type = "character", default = "", help = "Prefix of the reference file"),
  make_option("--cores",     type = "integer",   default = 1,    help = "Number of cores to use for parallel processing")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

max.len.gap = 20000

# print(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

aln.type = aln.type.msa

# ***********************************************************************
# ---- Values of parameters ----

# Set the number of cores for parallel processing
num.cores <- opt$cores
if(num.cores > 1){
  myCluster <- makeCluster(num.cores, type = "PSOCK")
  registerDoParallel(myCluster)  
}

# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop('Consensus folder doesn’t exist')

# Reference genome
ref.name <- opt$ref

# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type)
s.combinations <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
s.combinations = gsub(aln.type, "", s.combinations)
s.combinations = gsub(".h5", "", s.combinations)

if(ref.name != ""){
  ref.suff = paste0('_', ref.name)
  
  pokaz('Reference:', ref.name)
  s.combinations <- s.combinations[grep(ref.suff, s.combinations)]
  s.combinations = gsub(ref.suff, "", s.combinations)
  
} else {
  ref.suff = ''
}

if(length(s.combinations) == 0){
  stop('No Combinations found.')
} else {
  pokaz('Combinations', s.combinations)  
}

# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(s.comb, echo = T){
  
  file.comb = paste0(path.cons, aln.type, s.comb, ref.suff,'.h5')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  
  # Create group for blocks
  suppressMessages({
    h5createGroup(file.comb, gr.blocks)
  })

  # ---- Define coverage of blocks for every position  and put NA between blocks----
  idx.blocks = 0
  for(acc in accessions){
    
    pokaz('Accession', acc, 'combination', s.comb)
  
    v.init = h5read(file.comb, paste0(gr.accs.e, acc))
    v = v.init
    
    # ----  Find breaks  ----
    
    # Find blocks of additional breaks
    v = cbind(v, 1:length(v))                 # 2 - in ref-based coordinates
    v = v[(v[,1] != 0) & !(is.na(v[,1])),]    # 1 - existing coordinates of accessions
    
    idx.block.tmp = which(abs(diff(v[,1])) != 1)
    idx.block.beg = v[c(1, idx.block.tmp+1), 2]
    idx.block.end = v[c(idx.block.tmp, nrow(v)), 2]
  
    pokaz('Number of blocks init:', length(idx.block.beg))

    # Add NA
    v.block = rep(0, length(v.init))
    for(i.bl in 1:length(idx.block.beg)){
      v.block[idx.block.beg[i.bl]:idx.block.end[i.bl]] = i.bl
    }
    v.init[v.block == 0] = NA
    
    # Save
    suppressMessages({
      h5write(v.init, file.comb, paste0(gr.accs.e, acc))
    })
    
    # Remember for common breaks
    idx.blocks = idx.blocks + (v.block != 0) * 1
    
    rmSafe(v)
    rmSafe(v.init)
    rmSafe(v.block)
  }
  
  # Find all blocks
  df.blocks.tmp = c()
  idx.no.blocks = idx.blocks <= 1
  for(acc in accessions){
    
    pokaz('Accession', acc, 'combination', s.comb)
    
    # Read correspondence
    v.init = h5read(file.comb, paste0(gr.accs.e, acc))
    
    # Add NA as common breaks
    v.init[idx.no.blocks] = NA
    suppressMessages({
      h5write(v.init, file.comb, paste0(gr.accs.e, acc))
    })
    
    pokaz('Content NA:', sum(is.na(v.init)), 'zeros:', sum(!is.na(v.init) & (v.init ==0)), 'pos:', sum(!is.na(v.init) & (v.init !=0)))
    
    # ----  Find blocks of non-NA  ----
    # not.na <- !is.na(v.init)
    # changes <- diff(c(FALSE, not.na, FALSE))
    # idx.block.beg <- which(changes == 1)
    # idx.block.end <- which(changes == -1) - 1
    
    # Find blocks of additional breaks
    v = v.init
    v = cbind(v, 1:length(v))                       # 2 - in ref-based coordinates
    v = v[(v[,1] != 0) & !(is.na(v[,1])),]                                   # 1 - existing coordinates of accessions
    
    idx.block.tmp = which(abs(diff(v[,1])) != 1)
    idx.block.beg = v[c(1, idx.block.tmp+1), 2]
    idx.block.end = v[c(idx.block.tmp, nrow(v)), 2]
    
    
    pokaz('Number of blocks', length(idx.block.beg))
    
    
    # Blocks for the current accession
    df.blocks.acc = data.frame(pan.b = idx.block.beg,
                               pan.e = idx.block.end,
                               own.b = abs(v.init[idx.block.beg]),
                               own.e = abs(v.init[idx.block.end]),
                               acc = acc,
                               chr = strsplit(s.comb, '_')[[1]][1],
                               dir = (1 - sign(v.init[idx.block.beg])) / 2)
    df.blocks.acc[df.blocks.acc$dir == 1, 1:4] = df.blocks.acc[df.blocks.acc$dir == 1, c(2,1,4,3)]
    
    # Save for the whole dataset
    df.blocks.tmp = rbind(df.blocks.tmp, df.blocks.acc)
    
    # Define block IDs
    v.block = rep(0, length(v.init))
    for(i.bl in 1:length(idx.block.beg)){
      v.block[idx.block.beg[i.bl]:idx.block.end[i.bl]] = i.bl
    }
    
    # Save blocks
    suppressMessages({
      h5write(v.block, file.comb, paste0(gr.blocks, acc))
    })
    
    rmSafe(x.corr)
    rmSafe(x)
    rmSafe(v)
    rmSafe(v.init)
    rmSafe(idx.tmp.acc)
  }
  
  H5close()
  gc()
  
  return(df.blocks.tmp)

}


# ***********************************************************************
# ---- Loop  ----


if(num.cores == 1){
  list.blocks = list()
  for(s.comb in s.combinations){
    list.blocks[[s.comb]] = loop.function(s.comb)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  list.blocks = foreach(s.comb = s.combinations, .packages=c('rhdf5', 'crayon'))  %dopar% { 
    tmp = loop.function(s.comb)
    return(tmp)
  }
  stopCluster(myCluster)
}


# columns:   pan.b    pan.e    own.b    own.e   acc chr dir

df.blocks = do.call(rbind, list.blocks)
file.blocks = paste0(path.cons, aln.type,'blocks',ref.suff,'.rds')
saveRDS(df.blocks, file.blocks, compress = F)

warnings()




