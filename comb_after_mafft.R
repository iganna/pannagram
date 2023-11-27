
library(Biostrings)
library('seqinr')
library('foreach')
library(doParallel)
library("optparse")
source('utils.R')

# rm -rf gaps mob no_interesting class pos solid_aln


myCluster <- makeCluster(15, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

pokazStage('Combine all alignments together')

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.mafft.in"), type="character", default=NULL, 
              help="path to directory, where to combine fasta files for mafft runs", metavar="character"),
  make_option(c("-p", "--ref.pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("--path.mafft.out"), type="character", default=NULL, 
              help="path to directory, where to mafft results are", metavar="character"),
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to directory with the consensus", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer")
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


if (!is.null(opt$path.mafft.in)) path.mafft.in <- opt$path.mafft.in
if (!is.null(opt$path.mafft.out)) path.mafft.out <- opt$path.mafft.out
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons

n.flank = 30


# ---- Combinations of chromosomes query-base to create the alignments ----

# path.cons = './'
# ref.pref = '0'

s.pattern <- paste("^", 'res_', ".*", '_ref_', ref.pref, sep = '')
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub("res_", "", files)
pref.combinations <- sub("_ref.*$", "", pref.combinations)

pokaz('Reference:', ref.pref)
pokaz('Combinations', pref.combinations)


# ------------------------------------
# ------------------------------------
# flag.for = F
# ref = foreach(i.f = 1:length(fasta.files), .packages=c('stringr','Biostrings', 'R.utils'))  %dopar% { 

for(s.comb in pref.combinations){
  
  pokaz('* Combination', s.comb)
  
  # ---- All MAFFT results for the combination ----
  pref = paste('Gap', s.comb, sep = '_')
  mafft.res = data.frame(file = list.files(path = path.mafft.out, 
                                           pattern = paste('^', pref, '.*_flank_', n.flank, '_aligned.fasta$', sep='')))
  
  
  mafft.res$comb = lapply(mafft.res$file, function(s) strsplit(s, '_')[[1]])
  
  y = lapply(mafft.res$file, function(s) strsplit(s, '_')[[1]])
  z = t(matrix(unlist(y), ncol = length(y)))
  mafft.res$comb = paste(z[,2], z[,3], sep = '_')
  mafft.res$beg = as.numeric(z[,5])
  mafft.res$end = as.numeric(z[,6])
  mafft.res$id = as.numeric(z[,3])
  
  
  # ---- Short alignments ----
  msa.res = readRDS(paste(path.cons, 'aln_short_', s.comb, '.rds', sep = ''))
  msa.res$beg = msa.res$idx.gap.pos

  
  single.res = readRDS(paste(path.cons, 'singletons_', s.comb, '.rds', sep = ''))
  
}


# 
# for.flag = T
# for(aln.file in 1:length(fasta.files)){
#   
#   pos.file = paste(path.pos, gsub("fasta", "txt", fasta.files[i.f]), sep = '')
#   if(file.exists(pos.file)){
#     if(for.flag) next
#     return(NULL)
#   }    
#   
#   # Get sequences
#   z.fasta = paste(path.fasta, fasta.files[i.f],sep='')
#   
#     
#     if(!file.exists(aln.fasta)) {
#       stop()
#       if(for.flag) next
#       return(NULL)
#     }
#   
#   
#   # Try to read the alignment
#   aln = NULL
#   tryCatch(
#     expr = {
#       aln = ape::as.alignment(seqinr::read.fasta(aln.fasta))
#     },
#     error = function(e){
#       aln = NULL
#     }
#   )
#   if(is.null(aln)){
#     if(for.flag) next
#     return(NULL)
#   }
#   
#   # ------------------------
#   # Pipeline
#   # ------------------------
#   # Get sequence matrix
#   seqs.names = aln$nam
#   seqs.acc = sapply(seqs.names, function(s) strsplit(s, '_')[[1]][2])
#   seqs.aln = aln$seq
#   seqs.mx = c()
#   seqs.n = length(seqs.aln)
#   for(i in 1:seqs.n){
#     seqs.mx = rbind(seqs.mx, strsplit(seqs.aln[i],'')[[1]])
#   }
# 
#   # Remove flanking anchors
#   for(i in 1:nrow(seqs.mx)){
#     pos.non.zero = which(seqs.mx[i,] != '-')
#     seqs.mx[i,tail(pos.non.zero, (n.flank) )] = '-'
#     seqs.mx[i,pos.non.zero[1:(n.flank) ]] = '-'
#   }
#   seqs.mx = seqs.mx[,colSums(seqs.mx != '-') != 0, drop=F]
# 
# 
#   # Get positions of the alignment
#   n.pos = ncol(seqs.mx)
#   val.acc.pos = matrix(0, nrow=n.pos, ncol=n.acc)
#   colnames(val.acc.pos) = accessions
#   for(i.acc in 1:nrow(seqs.mx)){
# 
#     tmp = strsplit(seqs.names[i.acc],'_')[[1]]
#     pos = as.numeric(tmp[3]):as.numeric(tmp[4])
#     # pos = pos[-c(1,length(pos))]
# 
#     s = seqs.mx[i.acc,]
# 
#     acc = tmp[2]
#     val.acc.pos[s != '-', acc] = pos
#   }
#   tmp = as.numeric(strsplit(fasta.files[i.f], '_')[[1]][3])
#   idx.gap.pos = tmp
#   val.acc.pos = cbind(val.acc.pos, idx.gap.pos + (1:n.pos) / (n.pos+1))
# 
#   write.table(val.acc.pos, file = pos.file, quote = F, row.names = F, col.names = F)
#   rm(val.acc.pos)
# }
# 
# 
# 
# 
