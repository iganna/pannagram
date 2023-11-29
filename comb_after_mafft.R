
library(Biostrings)
library(rhdf5)
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

gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'

max.block.elemnt = 3 * 10^ 6


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
  
  # Get accessions
  file.comb = paste(path.cons, 'res_', s.comb,'_ref_',ref.pref,'.h5', sep = '')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)
  base.len = length(h5read(file.comb, paste(gr.accs.e, accessions[1], sep = '')))
  
  
  n.rows.block = round(max.block.elemnt / n.acc)
  seq.blocks = seq(1, base.len - n.rows.block, n.rows.block)
  
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
  
  
  # ---- Get long alignment positions ----
  mafft.aln.pos = list()
  for(i in 1:nrow(mafft.res)){
    aln.seq = readFastaMy(paste(path.mafft.out, mafft.res$file[i], sep = ''))
    n.aln.seq = length(aln.seq)
    name.aln.seq = names(aln.seq)
    pos.aln = sapply(name.aln.seq, function(s) strsplit(s, '\\|')[[1]][3:4])
    
    
    len.aln.seq <- nchar(aln.seq[1])
    # aln.mx <- matrix(, nrow = n.aln.seq, ncol = len.aln.seq)
    pos.mx <- matrix(0, nrow = n.aln.seq, ncol = len.aln.seq)
    for (i.seq in 1:length(aln.seq)) {
      tmp = strsplit(aln.seq[i.seq], "")[[1]]
      tmp.nongap = which(tmp != '-')
      # tmp[tail(tmp.nongap, (n.flank) )] = '-'
      # tmp[tmp.nongap[1:(n.flank) ]] = '-'
      # aln.mx[i, ] <- tmp
      
      tmp.nongap = tmp.nongap[-(1:(n.flank))]
      tmp.nongap <- tmp.nongap[1:(length(tmp.nongap) - n.flank)]
      
      p1 = pos.aln[1, i.seq]
      p2 = pos.aln[2, i.seq]
      pos.tmp = p1:p2
      pos.mx[i.seq, tmp.nongap] = pos.tmp
    }
    # aln.mx = aln.mx[,colSums(aln.mx != '-') != 0]
    pos.mx = pos.mx[,colSums(pos.mx != 0) != 0]
    
    mafft.aln.pos[[i]] = pos.mx
  }
  mafft.res$len = unlist(lapply(mafft.aln.pos, ncol))
  mafft.res$extra = mafft.res$len - (mafft.res$end - mafft.res$beg - 1)
  if(min(mafft.res$extra) < 0) stop('Wrong lengths of alignment and gaps')
  
  # ---- Short alignments ----
  msa.res = readRDS(paste(path.cons, 'aln_short_', s.comb, '.rds', sep = ''))
  msa.res$len = unlist(lapply(msa.res$aln, nrow))
  msa.res$extra = msa.res$len - (msa.res$ref.pos$end - msa.res$ref.pos$beg - 1)
  if(min(msa.res$extra) < 0) stop('Wrong lengths of alignment and gaps')

  # ---- Singletons alignments ----
  single.res = readRDS(paste(path.cons, 'singletons_', s.comb, '.rds', sep = ''))
  single.res$len = rowSums(single.res$pos.end) - rowSums(single.res$pos.beg)  + 1
  single.res$extra = single.res$len - (single.res$ref.pos$end - single.res$ref.pos$beg - 1)
  if(min(single.res$extra) < 0) stop('Wrong lengths of alignment and gaps')
  
  
  # ---- Analysis of positions ----
  # Here I wouls like fo find function of positions corresponcences between 4 things: 
  # old coordinates, long, short and singleton coordinates
  
  n.shift = rep(0, base.len)
  
  n.shift[mafft.res$end] = mafft.res$extra  # Long extra
  n.shift[msa.res$ref.pos$end] = msa.res$extra  # Short extra
  n.shift[single.res$ref.pos$end] = single.res$extra # Singletons extra
  n.shift = cumsum(n.shift)
  
  
  
  fp.main
  
  
  
  
  
  
  
  
  pos.delete = rep(0, base.len)
  pos.delete[single.res$ref.pos$beg] = 1
  pos.delete[single.res$ref.pos$end] = -1
  pos.delete[ msa.res$ref.pos$beg] = 1
  pos.delete[ msa.res$ref.pos$end] = -1
  pos.delete[ mafft.res$beg] = 1
  pos.delete[ mafft.res$end] = -1
  
  pos.delete = cumsum(pos.delete)
  
  pos.delete[single.res$ref.pos$beg] = 0
  pos.delete[single.res$ref.pos$end] = 0
  pos.delete[ msa.res$ref.pos$beg] = 0
  pos.delete[ msa.res$ref.pos$end] = 0
  pos.delete[ mafft.res$beg] = 0
  pos.delete[ mafft.res$end] = 0
  
  pos.remain = which(pos.delete == 0)
  
  # ---- Define which reals positions they have ----
  # pos.delete = rep(F, base.len)  # which are replace by alignments
  
  # Singletons
  pos.single = list()
  for(i in 1:length(single.res$len)){
    n.pos = single.res$len[i]
    pos.single[[i]] = single.res$ref.pos$beg[i] + (1:n.pos) / (n.pos+1)
    # pos = single.res$ref.pos$beg[i]: single.res$ref.pos$end[i]
    # pos = pos[-c(1, length(pos))]
    # pos.delete[pos] = T
  }
  
  
  # Short
  pos.short = list()
  for(i in 1:length(msa.res$len)){
    n.pos = msa.res$len[i]
    pos.short[[i]] = msa.res$ref.pos$beg[i] + (1:n.pos) / (n.pos+1)
    # pos = msa.res$ref.pos$beg[i]: msa.res$ref.pos$end[i]
    # pos = pos[-c(1, length(pos))]
    # pos.delete[pos] = T
  }
  
  # Long
  pos.long = list()
  for(i in 1:length(mafft.aln.pos)){
    n.pos = ncol(mafft.aln.pos[[i]])
    pos.long[[i]] = mafft.res$beg[i] + (1:n.pos) / (n.pos+1)
    # pos = mafft.res$beg[i]:mafft.res$end[i]
    # pos = pos[-c(1, length(pos))]
    # pos.delete[pos] = T
  }
  
  pos.additional = c(unlist(pos.single), unlist(pos.short), unlist(pos.long))
  
  
  # if(sum(duplicated(pos.delete)) != 0) stop('')
  
  # define, which positions to delete
  
  
  # ---- Define blocks before the big alignments ----
  pos.beg = mafft.res$beg
  pos.beg.bins <- cut(pos.beg, breaks = c(seq.blocks, Inf), labels = FALSE)
  pos.block.end = tapply(pos.beg, pos.beg.bins, max)
  pos.block.end[length(pos.block.end)] = base.len
  
  for(acc in accessions){
    v = h5read(file.comb, paste(gr.accs.e, accessions[1], sep = ''))
    v = cbind(v, 1:nrow(v))
    v = v[!pos.delete, ]
    v.aln = c()
    v.pos = c()
    for(i in 1:length(single.res$len)){
      
      if(single.res$pos.beg[i, acc] == 0){
        
      } else {
        
      }
      
      
    }
    # add schort
    
    # add long
      
    for(i.bl in 1:length(pos.block.end)){
      # add singletons
    }
  }
  
}

