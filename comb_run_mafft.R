
library(Biostrings)
library('seqinr')
# library('spgs')
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
  make_option(c("--path.mafft.in"), type="character", default=NULL, 
              help="path to directory, where to combine fasta files for mafft runs", metavar="character"),
  make_option(c("--path.mafft.out"), type="character", default=NULL, 
              help="path to directory, where to mafft results are", metavar="character"),
  make_option(c("--path.mafft.pos"), type="character", default=NULL, 
              help="path to directory, where to mafft results are", metavar="character"),
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


if (!is.null(opt$path.mafft.in)) path.mafft.in <- opt$path.mafft.in
if (!is.null(opt$path.mafft.out)) path.mafft.out <- opt$path.mafft.out
if (!is.null(opt$path.mafft.pos)) path.mafft.pos <- opt$path.mafft.pos

n.flank = 30

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------

# accessions = readRDS(paste(path.work, 'accessions.rds', sep = ''))
# n.acc = length(accessions)
# print(accessions)

# ===========================================================================

dist.thresh = 0.80
sim.thresh = 0.2
min.len = 15
n.flank = 30
max.len = 21000


n.wnd = 15
p.mismatches = 0.2  # more mismatches - problem


# Working directodies
path.fasta = path.mafft.in
path.mafft = path.mafft.out
path.pos = path.mafft.pos
if(!dir.exists(path.mafft)) system(paste('mkdir ', path.mafft, sep = ''))
if(!dir.exists(path.pos)) system(paste('mkdir ', path.pos, sep = ''))

# All fasta
fasta.files = list.files(path = path.fasta, pattern = paste('_flank_',n.flank,'.fasta$',sep=''))

# ------------------------------------
# ------------------------------------
# flag.for = F
# ref = foreach(i.f = 1:length(fasta.files), .packages=c('stringr','Biostrings', 'R.utils'))  %dopar% { 
  
for.flag = T
for(i.f in 1:length(fasta.files)){
  
  
  pos.file = paste(path.pos, gsub("fasta", "txt", fasta.files[i.f]), sep = '')
  if(file.exists(pos.file)){
    if(for.flag) next
    return(NULL)
  }    
  
  # Get sequences
  z.fasta = paste(path.fasta, fasta.files[i.f],sep='')
  
  pokaz('---------------------')
  
  # Try to run the alignment
  aln.fasta = paste(path.mafft, gsub("fasta", "txt", fasta.files[i.f]), sep = '')
  pokaz(aln.fasta)
  if(!file.exists(aln.fasta)) {
    #if(T) {
    # print(aln.fasta)
    # res <- R.utils::withTimeout({
    #   # --genafpair
    #   system(paste('mafft --op 5 --quiet --maxiterate 100 ', z.fasta, '>', aln.fasta,  sep = ' '))
    #   # return(T)
    # }, timeout = 100, onTimeout = "warning")
    
    
    timeout <- 600 # 10 минут в секундах
    timeout <- 0.5
    cmd <- paste('mafft --op 5 --quiet --maxiterate 100 ', z.fasta, '>', aln.fasta,  sep = ' ')
    system(command = cmd, timeout = timeout)
    # tryCatch({
    #   system2(command = cmd, stdout = aln.fasta, timeout = timeout)
    # }, error = function(e) {
    #   if(grepl("reached elapsed time limit", e$message)) {
    #     pokaz(z.fasta)
    #     cat("Команда была прервана после 10 минут выполнения\n")
    #   } else {
    #     cat("Произошла другая ошибка: ", e$message, "\n")
    #   }
    # })
    
    
    if(!file.exists(aln.fasta)) {
      stop()
      if(for.flag) next
      return(NULL)
    }
  }
  
  # # Try to read the alignment
  # aln = NULL
  # tryCatch(
  #   expr = {
  #     aln = ape::as.alignment(seqinr::read.fasta(aln.fasta))
  #   },
  #   error = function(e){ 
  #     aln = NULL
  #   }
  # )
  # if(is.null(aln)){
  #   if(for.flag) next
  #   return(NULL)
  # }
  
  # # ------------------------
  # # Pipeline
  # # ------------------------
  # # Get sequence matrix
  # seqs.names = aln$nam
  # seqs.acc = sapply(seqs.names, function(s) strsplit(s, '_')[[1]][2])
  # seqs.aln = aln$seq
  # seqs.mx = c()
  # seqs.n = length(seqs.aln)
  # for(i in 1:seqs.n){
  #   seqs.mx = rbind(seqs.mx, strsplit(seqs.aln[i],'')[[1]])
  # }
  # 
  # # Remove flanking anchors
  # for(i in 1:nrow(seqs.mx)){
  #   pos.non.zero = which(seqs.mx[i,] != '-')
  #   seqs.mx[i,tail(pos.non.zero, (n.flank) )] = '-'
  #   seqs.mx[i,pos.non.zero[1:(n.flank) ]] = '-'
  # }
  # seqs.mx = seqs.mx[,colSums(seqs.mx != '-') != 0, drop=F]
  # 
  # 
  # # Get positions of the alignment 
  # n.pos = ncol(seqs.mx)
  # val.acc.pos = matrix(0, nrow=n.pos, ncol=n.acc)
  # colnames(val.acc.pos) = accessions
  # for(i.acc in 1:nrow(seqs.mx)){
  #   
  #   tmp = strsplit(seqs.names[i.acc],'_')[[1]]
  #   pos = as.numeric(tmp[3]):as.numeric(tmp[4])
  #   # pos = pos[-c(1,length(pos))]
  #   
  #   s = seqs.mx[i.acc,]
  #   
  #   acc = tmp[2]
  #   val.acc.pos[s != '-', acc] = pos
  # }
  # tmp = as.numeric(strsplit(fasta.files[i.f], '_')[[1]][3])
  # idx.gap.pos = tmp
  # val.acc.pos = cbind(val.acc.pos, idx.gap.pos + (1:n.pos) / (n.pos+1))
  # 
  # write.table(val.acc.pos, file = pos.file, quote = F, row.names = F, col.names = F)
  # rm(val.acc.pos)
}




