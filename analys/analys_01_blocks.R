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
              help = "number of cores to use for parallel processing", metavar = "integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

print(opt)

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

# ---- Combinations of chromosomes query-base to create the alignments ----

# path.cons = './'
# ref.pref = '0'
# library(rhdf5)
# source('../../../pannagram/utils.R')

s.pattern <- paste("^", 'msa_', ".*", '_ref_', ref.pref, sep = '')
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub("msa_", "", files)
pref.combinations <- sub("_ref.*$", "", pref.combinations)

pokaz('Reference:', ref.pref)
pokaz('Combinations', pref.combinations)

# ----  Combine correspondence  ----

gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'
max.len.gap = 20000

gr.blocks = 'blocks/'


# columns:   pan.b    pan.e    own.b    own.e   acc chr dir
df.blocks = c()

#flag.for = F
#tmp = foreach(s.comb = pref.combinations, .packages=c('rhdf5', 'crayon'))  %dopar% {  # which accession to use
flag.for = T
for(s.comb in pref.combinations){

  
  file.comb = paste(path.cons, 'msa_', s.comb,'_ref_',ref.pref,'.h5', sep = '')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  
  # Create group for blocks
  h5createGroup(file.comb, gr.blocks)
  
  
  # ---- Define positions without blocks ----
  idx.blocks = 0
  for(acc in accessions){
    
    pokaz('Accession', acc, 'combination', s.comb)
  
    v.init = h5read(file.comb, paste(gr.accs.e, acc, sep = ''))
    v = v.init
    
    # ----  Find breaks  ----
    
    # Find blocks of additional breaks
    v = cbind(v, 1:length(v))                       # 2 - in ref-based coordinates
    v = v[v[,1] != 0,]                                   # 1 - existing coordinates of accessions
    v = cbind(v, 1:nrow(v))                       # 3 - ranked order in ref-based coordinates
    v = cbind(v, rank(abs(v[,1])) * sign(v[,1]))  # 4 - signed-ranked-order in accessions coordinates 
    
    # Save blocks
    idx.block.tmp = which(abs(diff(v[,4])) != 1)
    idx.block.beg = v[c(1, which(abs(diff(v[,4])) != 1)+1), 2]
    idx.block.end = v[c(which(abs(diff(v[,4])) != 1), nrow(v)), 2]
    pokaz('Number of blocks', length(idx.block.beg))
    v.block = rep(0, length(v.init))
    for(i.bl in 1:length(idx.block.beg)){
      v.block[idx.block.beg[i.bl]:idx.block.end[i.bl]] = i.bl
    }
    
    # suppressMessages({
    #   h5write(v.block, file.comb, paste(gr.blocks, acc, sep = ''))
    # })
    
    idx.blocks = idx.blocks + (v.block != 0) * 1
    
    rmSafe(v)
    rmSafe(v.init)
    rmSafe(v.block)
  }
  
  # Find all blocks
  idx.no.blocks = idx.blocks <= 1
  for(acc in accessions){
    
    pokaz('Accession', acc, 'combination', s.comb)
    
    v.init = h5read(file.comb, paste(gr.accs.e, acc, sep = ''))
    v = v.init
    v[idx.no.blocks] = 0
    
    # ----  Find breaks  ----
    
    # Find blocks of additional breaks
    v = cbind(v, 1:length(v))                       # 2 - in ref-based coordinates
    v = v[v[,1] != 0,]                                   # 1 - existing coordinates of accessions
    v = cbind(v, 1:nrow(v))                       # 3 - ranked order in ref-based coordinates
    v = cbind(v, rank(abs(v[,1])) * sign(v[,1]))  # 4 - signed-ranked-order in accessions coordinates 
    
    # Save blocks
    idx.block.tmp = which(abs(diff(v[,4])) != 1)
    idx.block.beg = v[c(1, which(abs(diff(v[,4])) != 1)+1), 2]
    idx.block.end = v[c(which(abs(diff(v[,4])) != 1), nrow(v)), 2]
    pokaz('Number of blocks', length(idx.block.beg))
    
    
    df.blocks.acc = data.frame(pan.b = idx.block.beg,
                               pan.e = idx.block.end,
                               own.b = abs(v.init[idx.block.beg]),
                               own.e = abs(v.init[idx.block.end]),
                               acc = acc,
                               chr = strsplit(s.comb, '_')[[1]][1],
                               dir = (1 - sign(v.init[idx.block.beg])) / 2)
    
    df.blocks.acc[df.blocks.acc$dir == 1, 1:4] = df.blocks.acc[df.blocks.acc$dir == 1, c(2,1,4,3)]
    
    df.blocks = rbind(df.blocks, df.blocks.acc)
    
    v.block = rep(0, length(v.init))
    for(i.bl in 1:length(idx.block.beg)){
      v.block[idx.block.beg[i.bl]:idx.block.end[i.bl]] = i.bl
    }
    
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
  
}

file.blocks = paste(path.cons, 'msa_blocks_ref_',ref.pref,'.rds', sep = '')
saveRDS(df.blocks, file.blocks, compress = F)




