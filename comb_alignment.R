library('foreach')
library(doParallel)
library("optparse")
source("synteny_infer.R")

args = commandArgs(trailingOnly=TRUE)


option_list = list(
  make_option(c("--path.work"), type="character", default=NULL, 
              help="path to working directory", metavar="character"),
  make_option(c("--path.base"), type="character", default=NULL, 
              help="path to base chromosomes", metavar="character"),
  make_option(c("--path.aln.pref"), type="character", default=NULL, 
              help="path with alignments", metavar="character"),
  
  make_option(c("--file.chr.len.ref"), type="character", default=NULL, 
              help="file with lengths of chromosomes of reference accessions", metavar="character"),
  make_option(c("--n.cores"), type="character", default=NULL, 
              help="numer of cores: 10 max", metavar="character"),
  
  make_option(c("--accs"), type="character", default=NULL, 
              help="accessions to create the alignment", metavar="character"),
  
  
  make_option(c("--n.chr.ref"), type="character", default=NULL, 
              help="number of chromosomes in the reference genome", metavar="character"),
  make_option(c("--n.chr.acc"), type="character", default=NULL, 
              help="number of chromosomes in the accessions", metavar="character"),
  
  make_option(c("--ref.acc"), type="character", default=NULL, 
              help="reference accessions", metavar="character"),
  
  make_option(c("--all.vs.all"), type="character", default=NULL, 
              help="alignment of all chromosomes vs all or not: T/F", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

print(opt)
# return()

if (!is.null(opt$path.work)) path.work <- opt$path.work

if (!is.null(opt$path.base)) path.base <- opt$path.base
if (!is.null(opt$ref.acc)) base.acc.ref <- opt$ref.acc

if (!is.null(opt$path.aln.pref)) path.pref <- opt$path.aln.pref

if (!is.null(opt$file.chr.len.ref)) file.chr.len.ref <- opt$file.chr.len.ref

if (!is.null(opt$n.chr.ref)) n.chr.ref <- opt$n.chr.ref
if (!is.null(opt$n.chr.acc)) n.chr.acc <- opt$n.chr.acc
if (!is.null(opt$all.vs.all)) all.vs.all <- as.logical(opt$all.vs.all)

accs <- opt$accs

if (!is.null(opt$n.cores)) {
  n.cores <- min(10, as.numeric(opt$n.cores))
} else {
  n.cores = 10
}

# ===========================================================================
# ===========================================================================

# 10 CORES MAXIMUM!!!!!!!!!!!
myCluster <- makeCluster(n.cores, # number of cores to use           # 10 CORES MAXIMUM!!!!!!!!!!!
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

# ------------------------------------------------------------------

# acc.ref = c('0', '10002', '10015','10024', '1741', '22005', '6046', '6124', '6244', '8236', '9537') 
# 0,10002,10015,10024,1741,22005,6046,6124,6244,8236,9537
print(path.pref)

path.aln.ref =  paste(path.pref, base.acc.ref, '/', sep = '')

# ------------------------------------------------------------------
# Get accession names
aln.pattern = '_postgap3.rds'
accessions = c()
for(p in path.aln.ref){
  accessions = c(accessions, unique(sapply(list.files(p, pattern=paste('*', aln.pattern, sep = '')), function(s) strsplit(s, '_')[[1]][1])))
}
acc.tbl = table(accessions)
accessions = names(acc.tbl)

if(!is.null(accs)){
  accessions = intersect(accessions, acc.tbl)
}

print(accessions)

# -------------------
# Combinations of chromosomes query-base to create the alignments
chromosome.pairs = c()
for(i.query in 1:length(accessions)){
  for(query.chr in 1:n.chr.acc){
    if(!all.vs.all){
      if(query.chr > n.chr.ref) next
      chromosome.pairs = rbind(chromosome.pairs, c(i.query, query.chr, query.chr))
      next
    }
    for(base.chr in 1:n.chr.ref){
      chromosome.pairs = rbind(chromosome.pairs, c(i.query, query.chr, base.chr))
    }
  }
}
chromosome.pairs = unique(chromosome.pairs[,c(2,3)])

# -------------------
# Length of reference chromosomes
if(!file.exists(file.chr.len.ref)){
  chr.len = c()
} else {
  chr.len = readRDS(file.chr.len.ref)
}

add.flag = F
for(i.chr in 1:n.chr.ref){
  acc = base.acc.ref
  # print(c(i.chr, acc))
  if(!is.null(chr.len)){
    idx = (chr.len$acc == acc) & (chr.len$chr == i.chr)
    if(sum(idx) > 0) next
  }
  add.flag = T
  g = seqinr::read.fasta(paste(path.base, acc, '_chr', i.chr, '.fasta', sep = ''))[[1]]
  chr.len = rbind(chr.len, data.frame(acc = acc, chr = i.chr, len = length(g)))
}
if(add.flag) saveRDS(chr.len, file.chr.len.ref)

# ------------------------------------------------------------------

print(base.acc.ref)
print(chromosome.pairs)

if(!dir.exists(path.work)) system(paste('mkdir ', path.work, sep = ''))

# flag.for = F
# tmp = foreach(i.chr.pair = 1:nrow(chromosome.pairs))  %dopar% {  # which accession to use
flag.for = T
for(i.chr.pair in 1:nrow(chromosome.pairs)){
  
  query.chr = chromosome.pairs[i.chr.pair, 1]
  base.chr = chromosome.pairs[i.chr.pair, 2]
  
  val.all = c()
  val.names = c()
  for(acc in accessions){
    print(acc)
    
    f.acc2ref = paste(path.aln.ref, acc, '_', query.chr, '_', base.chr, aln.pattern, sep = '')
    if(!file.exists(f.acc2ref)) next
    t.acc2ref = readRDS(f.acc2ref)
    base.len.ref = chr.len[(chr.len[,1] == base.acc.ref) & (chr.len[,2] == base.chr),3]
    val = getCorresp2BaseSign(t.acc2ref, base.len.ref)
    rm(t.acc2ref)
    
    val.all = cbind(val.all, val)
    val.names = c(val.names, acc)
    rm(val)
  }
  
  if(length(val.names) == 0){
    if(flag.for){
      next
    } else {
      return(NULL)  
    }
  }
  colnames(val.all) = val.names
  val.all = cbind(1:nrow(val.all),val.all)
  colnames(val.all)[1] = base.acc.ref
  
  file.out.consensus = paste(path.work, 'consensus_', query.chr, '_', base.chr, '_', base.acc.ref, '_direct.rds', sep = '')
  saveRDS(val.all, file.out.consensus, compress = F)
  rm(val.all)
  gc()
  
}




