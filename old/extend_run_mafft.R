
library(Biostrings)
library('seqinr')
library('spgs')
library('foreach')
library(doParallel)
library(R.utils)
library("optparse")

# rm -rf gaps mob no_interesting class pos solid_aln


myCluster <- makeCluster(15, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

message('=====================================================')
message('  Run MAFFT  ')
message('-----------------------------------------------------')

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--n.chr.acc"), type="character", default=NULL, 
              help="number of chromosomes in the accessions", metavar="character"),
  make_option(c("--path.chromosomes"), type="character", default=NULL, 
              help="path to the chromosomes", metavar="character"),
  make_option(c("--path.work"), type="character", default=NULL, 
              help="path with working files", metavar="character")
); 



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

print(opt)


if (!is.null(opt$n.chr.acc)) n.chr.acc <- opt$n.chr.acc
if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes
if (!is.null(opt$path.work)) path.work <- opt$path.work

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# 
# 
# n.chr.acc = 5
# path.chromosomes = '../col_chromosomes/'
# path.work = '../col_consensus/'

accessions = readRDS(paste(path.work, 'accessions.rds', sep = ''))
n.acc = length(accessions)
print(accessions)

# ===========================================================================

dist.thresh = 0.80
sim.thresh = 0.2
min.len = 15
n.flank = 30
max.len = 21000


n.wnd = 15
p.mismatches = 0.2  # more mismatches - problem




# Working directodies
path.fasta = paste(path.work, 'mafft_fasta/', sep='')
path.mafft = paste(path.work, 'mafft_result/', sep='')
path.pos = paste(path.work, 'mafft_pos/', sep='')
if(!dir.exists(path.mafft)) system(paste('mkdir ', path.mafft, sep = ''))
if(!dir.exists(path.pos)) system(paste('mkdir ', path.pos, sep = ''))

# All fasta
fasta.files = list.files(path = path.fasta, pattern = paste('_flank_',n.flank,'.fasta$',sep=''))

# ------------------------------------
# ------------------------------------
flag.for = F
ref = foreach(i.f = 1:length(fasta.files), .packages=c('stringr','Biostrings', 'spgs', 'R.utils'))  %dopar% { 
  
#flag.for = T
#for(i.f in 1:length(fasta.files)){
  
  
  pos.file = paste(path.pos, gsub("fasta", "txt", fasta.files[i.f]), sep = '')
  if(file.exists(pos.file)){
    if(flag.for){
      next
    } else {
      return(NULL)
    }
  }    
  
  # Get sequences
  z.fasta = paste(path.fasta, fasta.files[i.f],sep='')
  
  print('---------------------')
  
  # Try to run the alignment
  aln.fasta = paste(path.mafft, gsub("fasta", "txt", fasta.files[i.f]), sep = '')
  print(aln.fasta)
  if(!file.exists(aln.fasta)) {
    #if(T) {
    # print(aln.fasta)
    res <- R.utils::withTimeout({
      # --genafpair
      system(paste('mafft --op 5 --quiet --maxiterate 100 ', z.fasta, '>', aln.fasta,  sep = ' '))
      # return(T)
    }, timeout = 100, onTimeout = "warning")
    
    if(!file.exists(aln.fasta)) {
      if(flag.for){
        next
      } else {
        return(NULL)  
      }
    }
  }
  
  # Try to read the alignment
  aln = NULL
  tryCatch(
    expr = {
      aln = ape::as.alignment(seqinr::read.fasta(aln.fasta))
    },
    error = function(e){ 
      aln = NULL
    }
  )
  if(is.null(aln)){
    if(flag.for){
      next
    } else {
      return(NULL)  
    }
  }
  
  # ------------------------
  # Pipeline
  # ------------------------
  # Get sequence matrix
  seqs.names = aln$nam
  seqs.acc = sapply(seqs.names, function(s) strsplit(s, '_')[[1]][2])
  seqs.aln = aln$seq
  seqs.mx = c()
  seqs.n = length(seqs.aln)
  for(i in 1:seqs.n){
    seqs.mx = rbind(seqs.mx, strsplit(seqs.aln[i],'')[[1]])
  }
  
  # Remove flanking anchors
  for(i in 1:nrow(seqs.mx)){
    pos.non.zero = which(seqs.mx[i,] != '-')
    seqs.mx[i,tail(pos.non.zero, (n.flank) )] = '-'
    seqs.mx[i,pos.non.zero[1:(n.flank) ]] = '-'
  }
  seqs.mx = seqs.mx[,colSums(seqs.mx != '-') != 0, drop=F]
  
  
  # Get positions of the alignment 
  n.pos = ncol(seqs.mx)
  val.acc.pos = matrix(0, nrow=n.pos, ncol=n.acc)
  colnames(val.acc.pos) = accessions
  for(i.acc in 1:nrow(seqs.mx)){
    
    tmp = strsplit(seqs.names[i.acc],'_')[[1]]
    pos = as.numeric(tmp[3]):as.numeric(tmp[4])
    # pos = pos[-c(1,length(pos))]
    
    s = seqs.mx[i.acc,]
    
    acc = tmp[2]
    val.acc.pos[s != '-', acc] = pos
  }
  tmp = as.numeric(strsplit(fasta.files[i.f], '_')[[1]][3])
  idx.gap.pos = tmp
  val.acc.pos = cbind(val.acc.pos, idx.gap.pos + (1:n.pos) / (n.pos+1))
  
  write.table(val.acc.pos, file = pos.file, quote = F, row.names = F, col.names = F)
  rm(val.acc.pos)
}




