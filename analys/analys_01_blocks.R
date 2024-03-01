# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
})

source("utils/utils.R")

pokazStage('Define blocks in the alignemnt')

args = commandArgs(trailingOnly=TRUE)


option_list = list(
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to consensus directory", metavar="character"),
  make_option(c("-p", "--ref.pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("--aln.type"), type="character", default="default", 
              help="type of alignment ('msa_', 'comb_', 'v_', etc)", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# Set the number of cores for parallel processing
num.cores.max = 10
num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))
myCluster <- makeCluster(num.cores, type = "PSOCK")
registerDoParallel(myCluster)


# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop('Consensus folder doesn’t exist')

# Reference genome
if (is.null(opt$ref.pref)) {
  stop("ref.pref is NULL")
} else {
  ref.pref <- opt$ref.pref
}


# Alignment prefix
if (!is.null(opt$aln.type)) {
  aln.type = opt$aln.type
} else {
  aln.type = 'msa_'
}

# ---- Combinations of chromosomes query-base to create the alignments ----

# path.cons = './'
# ref.pref = '0'
# library(rhdf5)
# source('../../../pannagram/utils/utils.R')


s.pattern <- paste("^", aln.type, ".*", '_ref_', ref.pref, sep = '')
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type, "", files)
pref.combinations <- sub("_ref.*$", "", pref.combinations)
pref.combinations <- pref.combinations[grep("^[0-9]+_[0-9]+$", pref.combinations)]

pokaz('Reference:', ref.pref)
pokaz('Combinations', pref.combinations)

# ----  Combine correspondence  ----

gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'
max.len.gap = 20000

gr.blocks = 'blocks/'


# ---- Main ----
flag.for = F
list.blocks = foreach(s.comb = pref.combinations, .packages=c('rhdf5', 'crayon'))  %dopar% {

# flag.for = T
# list.blocks = list()
# for(s.comb in pref.combinations){

  
  file.comb = paste(path.cons, aln.type, s.comb,'_ref_',ref.pref,'.h5', sep = '')
  
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
  
    v.init = h5read(file.comb, paste(gr.accs.e, acc, sep = ''))
    v = v.init
    
    # ----  Find breaks  ----
    
    # Find blocks of additional breaks
    v = cbind(v, 1:length(v))                       # 2 - in ref-based coordinates
    v = v[(v[,1] != 0) & !(is.na(v[,1])),]                                   # 1 - existing coordinates of accessions
    
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
      h5write(v.init, file.comb, paste(gr.accs.e, acc, sep = ''))
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
    v.init = h5read(file.comb, paste(gr.accs.e, acc, sep = ''))
    
    # Add NA as common breaks
    v.init[idx.no.blocks] = NA
    suppressMessages({
      h5write(v.init, file.comb, paste(gr.accs.e, acc, sep = ''))
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
      h5write(v.block, file.comb, paste(gr.blocks, acc, sep = ''))
    })
    
    rmSafe(x.corr)
    rmSafe(x)
    rmSafe(v)
    rmSafe(v.init)
    rmSafe(idx.tmp.acc)
  }
  
  H5close()
  gc()
  
  if(flag.for){
    list.blocks[[s.comb]] = df.blocks.tmp
  } else {
    return(df.blocks.tmp)
  }
  
}

# print(list.blocks)

# columns:   pan.b    pan.e    own.b    own.e   acc chr dir

df.blocks = do.call(rbind, list.blocks)
file.blocks = paste(path.cons, aln.type,'blocks_ref_',ref.pref,'.rds', sep = '')
saveRDS(df.blocks, file.blocks, compress = F)




