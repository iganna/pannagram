suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
})



source("utils.R")
source("synteny_funcs.R")

pokazStage('Combine: alignments by chromosomes')

args = commandArgs(trailingOnly=TRUE)


option_list = list(
  make_option(c("--path.chr.len"), type="character", default=NULL, 
              help="file with lengths of chromosomes", metavar="character"),
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to consensus directory", metavar="character"),
  make_option(c("-o", "--path.aln"), type="character", default=NULL, 
              help="path to the output directory with alignments", metavar="character"),
  make_option(c("-p", "--pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("-n", "--n.chr.ref"), type="character", default=NULL, 
              help="number of chromosomes in the reference genome", metavar="character"),
  make_option(c("-m", "--n.chr.acc"), type="character", default=NULL, 
              help="number of chromosomes in the accessions", metavar="character"),
  make_option(c("-a", "--all.vs.all"), type="character", default=NULL, 
              help="alignment of all chromosomes vs all or not: T/F", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("-r", "--path.ref"), type="character", default=NULL, 
              help="path to the reference file", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="type of fasta files", metavar="character")
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
if(!dir.exists(path.cons)) system(paste('mkdir ', path.cons, sep = ''))

# File with chromosomal lengths, which should be in the consensus dir
path.chr.len = ifelse(!is.null(opt$path.chr.len), opt$path.chr.len, 'chr_len/')
path.chr.len = paste(path.cons, path.chr.len, sep = '')
if(!dir.exists(path.chr.len)) system(paste('mkdir ', path.chr.len, sep = ''))

if (!is.null(opt$path.ref)) path.base <- opt$path.ref  # to know the chromosomal lengths
if (!is.null(opt$type)) base.suff <- opt$type  # to read fasta file
if (!is.null(opt$pref)) base.acc.ref <- opt$pref

# Path with alignments
if (!is.null(opt$path.aln)) path.aln <- opt$path.aln

# To create combinations
if (!is.null(opt$n.chr.ref)) n.chr.ref <- as.numeric(opt$n.chr.ref)
if (!is.null(opt$n.chr.acc)) n.chr.acc <- as.numeric(opt$n.chr.acc)
if (!is.null(opt$all.vs.all)) all.vs.all <- as.logical(opt$all.vs.all)


# ---- Get accession names ----
aln.suff = '_maj.rds'
aln.files <- list.files(path = path.aln, 
                        pattern = sub("\\.rds", "\\\\.rds$", aln.suff))
accessions <- unique(sub("_(.*)", "", aln.files))

pokaz('Accessions:', accessions)

# ---- Combinations of chromosomes query-base to create the alignments ----
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

# ---- Length of reference chromosomes ----

file.chr.len = paste(path.chr.len, base.acc.ref, '_len.rds', sep = '')
if(!file.exists(file.chr.len)){
  for(i.chr in 1:n.chr.ref){
    acc = base.acc.ref
    # print(c(i.chr, acc))
    
    # Read base chromosome
    base.file = paste0(base.acc.ref, '_chr', i.chr , '.', base.suff, collapse = '')
    pokaz('Base:', base.file)
    base.fas.fw = readFastaMy(paste(path.base, base.file, sep = ''))
    
    chr.len = c(chr.len, length(base.fas.fw))
  }
  saveRDS(chr.len, file.chr.len)
} else {
  chr.len = readRDS(file.chr.len)
}

pokaz('Chromosomal lengths', chr.len)

# ------------------------------------------------------------------

pokaz('Reference:', base.acc.ref)

# flag.for = F
# tmp = foreach(i.chr.pair = 1:nrow(chromosome.pairs))  %dopar% {  # which accession to use
flag.for = T
for(i.chr.pair in 1:nrow(chromosome.pairs)){
  
  query.chr = chromosome.pairs[i.chr.pair, 1]
  base.chr = chromosome.pairs[i.chr.pair, 2]
  
  for(acc in accessions){
    
    pref.comb = paste0(acc, '_', query.chr, '_', base.chr, collapse = '')
    file.aln.full <- paste(path.aln, paste0(pref.comb,  '_full.rds', collapse = ''), sep = '')
    if(!file.exists(file.aln.full)) next
    pokaz('Aln exists:', file.aln.full)
    
    # 
    # 
    # print(acc)
    # 
    # f.acc2ref = paste(path.aln.ref, acc, '_', query.chr, '_', base.chr, aln.suff, sep = '')
    # if(!file.exists(f.acc2ref)) next
    # t.acc2ref = readRDS(f.acc2ref)
    # base.len.ref = chr.len[base.chr]
    # val = getCorresp2BaseSign(t.acc2ref, base.len.ref)
    # rm(t.acc2ref)
    # 
    # val.all = cbind(val.all, val)
    # val.names = c(val.names, acc)
    # rm(val)
  }
  # 
  # if(length(val.names) == 0){
  #   if(flag.for){
  #     next
  #   } else {
  #     return(NULL)  
  #   }
  # }
  # colnames(val.all) = val.names
  # val.all = cbind(1:nrow(val.all),val.all)
  # colnames(val.all)[1] = base.acc.ref
  # 
  # file.out.consensus = paste(path.cons, 'consensus_', query.chr, '_', base.chr, '_', base.acc.ref, '_direct.rds', sep = '')
  # saveRDS(val.all, file.out.consensus, compress = F)
  # rm(val.all)
  # gc()
  
}




