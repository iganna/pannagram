suppressMessages({ library(Biostrings)
  library(rhdf5)
  library('seqinr')
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
              help = "number of cores to use for parallel processing", metavar = "integer")
); 



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# # Set the number of cores for parallel processing
# num.cores.max = 10
# num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))
# myCluster <- makeCluster(num.cores, type = "PSOCK")
# registerDoParallel(myCluster)

# Reference genome
if (is.null(opt$ref.pref)) {
  stop("ref.pref is NULL")
} else {
  ref.pref <- opt$ref.pref
}


if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons

n.flank = 30

gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'


# ---- Testing ----
# 
library(rhdf5)
source('../../../pannagram/utils.R')
path.cons = './'
path.chromosomes = '/home/anna/storage/arabidopsis/pacbio/pan_test/tom2/chromosomes/'
ref.pref = '0'
nts = c('A', 'C', 'G', 'T', '-')




# ---- Combinations of chromosomes query-base to create the alignments ----


s.pattern <- paste("^", 'msa_', ".*", '_ref_', ref.pref, sep = '')
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub("msa_", "", files)
pref.combinations <- sub("_ref.*$", "", pref.combinations)

pokaz('Reference:', ref.pref)
pokaz('Combinations', pref.combinations)


nts = c('A', 'C', 'G', 'T', '-')


# ------------------------------------
# ------------------------------------
#flag.for = F
#tmp = foreach(s.comb = pref.combinations, .packages=c('rhdf5', 'crayon'))  %dopar% {  # which accession to use
flag.for = T
for(s.comb in pref.combinations){
  
  pokaz('* Combination', s.comb)
  
  # Get accessions
  file.comb = paste(path.cons, 'msa_', s.comb,'_ref_',ref.pref,'.h5', sep = '')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)
  
  # File with sequences
  file.seq = paste(path.cons, 'seq_', s.comb,'_ref_',ref.pref,'.h5', sep = '')
  if (file.exists(file.seq)) file.remove(file.seq)
  h5createFile(file.seq)
  h5createGroup(file.seq, gr.accs.e)
  
  
  mx.consensus = NULL
  idx.negative = c()
  for(acc in accessions){
    pokaz('Sequence of accession', acc)
    v = h5read(file.comb, paste(gr.accs.e, acc, sep = ''))
    if(is.null(mx.consensus)){
      mx.consensus = matrix(0, nrow = length(v), ncol = length(nts), dimnames = list(NULL, nts))
    }
    
    acc.name = gsub('acc_', '', acc)
    q.chr = strsplit(s.comb, '_')[[1]][1]
    genome = readFastaMy(paste(path.chromosomes, acc.name, '_chr', q.chr, '.fasta', sep = ''))
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
    
    for(s.nt in nts){
      mx.consensus[,s.nt] = mx.consensus[,s.nt] + (s == s.nt)
    }
    
    suppressMessages({
      h5write(s, file.seq, paste(gr.accs.e, acc, sep = ''))
    })
  }
  
  suppressMessages({
    h5write(mx.consensus, file.seq, 'matrix')
  })
  
  H5close()
  gc()
  
}

# stopCluster(myCluster)

warnings()


