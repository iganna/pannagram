suppressMessages({ library(Biostrings)
  library(rhdf5)
  library('foreach')
  library(doParallel)
  library("optparse")
})

source("utils/utils.R")

pokazStage('Get sequence alignments and consensus sequence')

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--ref.pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("--path.chromosomes"), type="character", default=NULL, 
              help="path to directory with chromosomes", metavar="character"),
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to directory with the consensus", metavar="character"),
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

if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons

path.seq = paste(path.cons, 'seq/', sep = '')
if (!dir.exists(path.seq)) dir.create(path.seq)


# ---- Variables ----

gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'


# ---- Testing ----
# 
# library(rhdf5)
# source('../../../pannagram/utils/utils.R')
# path.cons = './'
# path.chromosomes = '/home/anna/storage/arabidopsis/pacbio/pan_test/tom2/chromosomes/'
# ref.pref = '0'
# s.nts = c('A', 'C', 'G', 'T', '-')


# ---- Combinations of chromosomes query-base to create the alignments ----


s.pattern <- paste("^", aln.type, ".*", '_ref_', ref.pref, sep = '')
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type, "", files)
pref.combinations <- sub("_ref.*$", "", pref.combinations)
pref.combinations <- pref.combinations[grep("^[0-9]+_[0-9]+$", pref.combinations)]

pokaz('Reference:', ref.pref)
if(length(pref.combinations) == 0){
  stop('No Combinations found.')
} else {
  pokaz('Combinations', pref.combinations)  
}


s.nts = c('A', 'C', 'G', 'T', '-')


# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(s.comb, echo = T){
# tmp = foreach(s.comb = pref.combinations, .packages=c('rhdf5', 'crayon'))  %dopar% {  # which accession to use
# # for(s.comb in pref.combinations){
  
  pokaz('* Combination', s.comb)
  
  # Get accessions
  file.comb = paste(path.cons, aln.type, s.comb,'_ref_',ref.pref,'.h5', sep = '')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)
  
  # File with sequences
  file.seq = paste(path.seq, 'seq_', s.comb,'_ref_',ref.pref,'.h5', sep = '')
  if (file.exists(file.seq)) file.remove(file.seq)
  h5createFile(file.seq)
  h5createGroup(file.seq, gr.accs.e)
  
  
  mx.consensus = NULL
  idx.negative = c()
  for(acc in accessions){
    pokaz('Sequence of accession', acc)
    v = h5read(file.comb, paste(gr.accs.e, acc, sep = ''))
    v.na = is.na(v)
    v[v.na] = 0
    if(is.null(mx.consensus)){
      mx.consensus = matrix(0, nrow = length(v), ncol = length(s.nts), dimnames = list(NULL, s.nts))
    }
    
    q.chr = strsplit(s.comb, '_')[[1]][1]
    genome = readFastaMy(paste(path.chromosomes, acc, '_chr', q.chr, '.fasta', sep = ''))
    genome = seq2nt(genome)
    genome = toupper(genome)
  
    s = rep('-', length(v))
    idx.plus = (v > 0)
    idx.mins = (v < 0)
    if(sum(idx.plus) > 0){
      s[idx.plus] = genome[v[idx.plus]]
    }
    if(sum(idx.mins) > 0){
      s[idx.mins] = justCompl(genome[abs(v[idx.mins])])
    }
    
    idx.negative = c(idx.negative, which(idx.mins))
    
    for(s.nt in s.nts){
      mx.consensus[,s.nt] = mx.consensus[,s.nt] + (s == s.nt)
    }
    
    suppressMessages({
      h5write(s, file.seq, paste(gr.accs.e, acc, sep = ''))
    })
    
    rmSafe(v)
    rmSafe(v.na)
    rmSafe(genome)
    rmSafe(s)
    rmSafe(idx.plus)
    rmSafe(idx.mins)
    gc()
  }
  
  suppressMessages({
    h5write(mx.consensus, file.seq, 'matrix')
  })
  
  # ---- Consensus sequence ----
  pokaz('Prepare consensus fasta-sequence')
  i.chr = comb2ref(s.comb)
  file.seq.cons = paste(path.seq, 'seq_cons_', i.chr, '.fasta', sep = '')
  
  
  n = nrow(mx.consensus)
  s.cons = rep('N', n)
  n.nt = rep(0, n)
  for(k in 1:4){
    idx.k =  mx.consensus[,k] > n.nt
    s.cons[idx.k] = s.nts[k]
    n.nt[idx.k] = mx.consensus[idx.k, k]
  }
  
  if(sum(s.cons == 'N') > 0) pokazAttention('Some nucleotides are missed:', sum(s.cons == 'N'))
  s.cons = paste0(s.cons, collapse = '')
  
  pokaz('Saving consensus sequence...')
  names(s.cons) = paste('PanGen_Chr', i.chr, sep = '')
  writeFastaMy(s.cons, file.seq.cons)
  
  rmSafe(mx.consensus)
  rmSafe(s.cons)
  H5close()
  gc()
  
}


# ***********************************************************************
# ---- Loop  ----


if(num.cores == 1){
  list.blocks = list()
  for(s.comb in pref.combinations){
    list.blocks[[s.comb]] = loop.function(s.comb)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  list.blocks = foreach(s.comb = pref.combinations, .packages=c('rhdf5', 'crayon'))  %dopar% { 
    tmp = loop.function(s.comb)
    return(tmp)
  }
  stopCluster(myCluster)
}

warnings()














