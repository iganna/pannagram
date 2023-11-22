suppressMessages(library('foreach'))
suppressMessages(library(doParallel))
suppressMessages(library(msa))
suppressMessages(library(dplyr))
suppressMessages(library("optparse"))

# myCluster <- makeCluster(30, # number of cores to use
#                          type = "PSOCK") # type of cluster
# registerDoParallel(myCluster)




message('=====================================================')
message('  Prepare files for alignment  ')
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


for(i.chr in 1:n.chr.acc){
  
  # Read genomes
  file.genomes = paste(path.work, 'all_',i.chr,'_genomes.rds', sep = '')
  if(!file.exists(file.genomes)){
    message('Read genomes in loop...')
    genomes.all = list()
    for(acc in accessions){
      print(paste('Reading genome', acc))
      f.genome = paste(path.chromosomes, acc, '_chr', i.chr, '.fasta', sep = '')
      genome.fw = seqinr::read.fasta(f.genome)
      genome.fw = genome.fw[[1]] # do not remove this 1, file contains only one sequence
      genomes.all[[acc]] = genome.fw
    }  
    message('Save genomes..')
    saveRDS(genomes.all, file.genomes, compress = F)
  } else {
    message('Read saved genomes...')
    genomes.all = readRDS(file.genomes)  
  }
  
  # Read primary alignment
  message(paste('Chromosome', i.chr))
  
  # val.common.file = paste(path.work,'consensus_',i.chr, '_', i.chr, '.rds', sep = '')
  val.common.file = paste(path.work,'val_common_chr_',i.chr,'.rds', sep = '')
  v = abs(readRDS(val.common.file))
  
  # get blocks without inversions and translocations
  
  block.file = paste(path.work,'blocks_chr_',i.chr,'.rds', sep = '')
  if(!file.exists(block.file)){
    message('Get directed blocks....')
    pos.bl = c()
    v.block = v * 0
    for(i.acc in 1:n.acc){
      message(i.acc)
      v.acc = v[,i.acc]
      v.acc = cbind(v.acc, 1:length(v.acc))
      v.acc = v.acc[v.acc[,1] != 0,]
      min.pos.acc = v.acc[1,2]
      max.pos.acc = v.acc[nrow(v.acc),2]
      v.acc = v.acc[order(v.acc[,1]),]
      v.acc = cbind(v.acc, 1:nrow(v.acc))
      v.acc = v.acc[order(v.acc[,2]),]
      idx = which(abs(v.acc[-1,3] - v.acc[-nrow(v.acc),3]) != 1)
      pos.end.blocks = c(v.acc[idx,2], max.pos.acc)
      pos.beg.blocks = c(min.pos.acc,v.acc[idx + 1,2])
      
      pos.bl = c(pos.bl, pos.beg.blocks, pos.end.blocks)
      
      for(i.bl in 1:length(pos.end.blocks)){
        if(v[pos.beg.blocks[i.bl], i.acc] < v[pos.end.blocks[i.bl], i.acc]){
          v.block[pos.beg.blocks[i.bl]:pos.end.blocks[i.bl], i.acc] = i.bl
        } else {
          v.block[pos.beg.blocks[i.bl]:pos.end.blocks[i.bl], i.acc] = -i.bl
        }
      }
    }
    pos.bl = sort(unique(pos.bl))
    length(pos.bl)
    saveRDS(v.block, block.file, compress = F)
  } else {
    message('Read blocks....')
    v.block = readRDS(block.file)
  }
  
  
  
  # Define positions, which should be filled
  max.len = 25000
  message('Get positions to fill....')
  v.gap = v * 0
  for(i.acc in 1:n.acc){
    message(i.acc)
    v.acc = v[,i.acc]
    v.acc = cbind(v.acc, 1:length(v.acc))
    v.acc = v.acc[v.acc[,1] != 0,]
    v.acc = v.acc[order(v.acc[,1]),]
    idx = which(abs(v.acc[-1,1] - v.acc[-nrow(v.acc),1]) != 1)
    
    pos.beg.blocks = v.acc[idx,2]
    pos.end.blocks = v.acc[idx + 1, 2]
    
    pos.d = abs(v[pos.end.blocks, ] - v[pos.beg.blocks, ])
    idx.remain = ((pos.end.blocks - pos.beg.blocks) == 1) & (pos.d[,i.acc] <= max.len )
    
    pos.d[v[pos.end.blocks, ] == 0] = 0
    pos.d[v[pos.beg.blocks, ] == 0] = 0
    pos.d[v.block[pos.beg.blocks, ] != v.block[pos.end.blocks, ]] = 0
    pos.d = apply(pos.d, 1, max)
    idx.remain = (pos.d <= max.len) | idx.remain
    pos.beg.blocks = pos.beg.blocks[idx.remain]
    pos.end.blocks = pos.end.blocks[idx.remain]
    
    # pos.d = rowSums(v.block[pos.beg.blocks, ] != v.block[pos.end.blocks, ])
    # pos.beg.blocks = pos.beg.blocks[pos.d == 0]
    # pos.end.blocks = pos.end.blocks[pos.d == 0]
    
    if(length(pos.end.blocks) == 0) next
    
    for(i.bl in 1:length(pos.end.blocks)){
      # if in the same block
      if(v.block[pos.beg.blocks[i.bl], i.acc] != v.block[pos.end.blocks[i.bl], i.acc]) next
      v.gap[pos.beg.blocks[i.bl]:pos.end.blocks[i.bl], i.acc] = 1
    }
  }
  
  idx.gap = (rowSums(v.gap) != 0) * 1
  n = length(idx.gap)
  idx.beg = which((idx.gap[-1] == 1) & (idx.gap[-n] == 0)) + 1
  idx.end = which((idx.gap[-1] == 0) & (idx.gap[-n] == 1))
  idx.gap = cbind(idx.beg, idx.end)
  
  x = v[idx.gap[,1],]
  y = v[idx.gap[,2],]
  
  x.bl = v.block[idx.gap[,1],]
  y.bl = v.block[idx.gap[,2],]
  
  d = abs(x - y) - 1
  d[x == 0] = 0
  d[y == 0] = 0
  d[x.bl != y.bl] = 0
  
  # d.max = apply(d, 1, max)
  # idx.len.filter = (d.max <= max.len)
  
  d[d > max.len] = 0
  idx.len.filter = (rowSums(d) != 0)
  
  
  idx.gap = idx.gap[idx.len.filter,]
  x = x[idx.len.filter,]
  y = y[idx.len.filter,]
  d = d[idx.len.filter,]
  x.bl = x.bl[idx.len.filter,]
  y.bl = y.bl[idx.len.filter,]
  
  # Positions are defined!
  # --------------------------------------------------
  
  # 
  # d.sign = (y > x) * 2 - 1
  # d.sign[d == 0] = 0
  # d.sign = d.sign * x.bl
  # d[d.sign < 0] = 0
  
  x[d == 0] = 0
  y[d == 0] = 0
  
  message(paste('LONG sequences', sum(d > max.len)))
  
  d.max = apply(d, 1, max)
  d.zero = rowSums(d == 0)
  
  idx.single = which(d.zero == (n.acc - 1))
  idx.aln = setdiff(1:length(d.max), idx.single)
  
  # Multiple alignment of short indels
  # gap.pos = idx.gap[,2] - idx.gap[,1] - 1
  
  aln.seqs = vector("list", max(idx.aln))
  aln.pos = vector("list", max(idx.aln))
  for(i.gap in idx.aln){
    # cat('.')
    if(round(i.gap/100) * 100 == i.gap) {
      cat(i.gap)
      cat(' ')
    }
    d.m = d[i.gap,]
    seqs = c()
    pos.idx = list()
    for(acc in accessions){
      if(d.m[acc] == 0) next
      
      s.idx = (x[i.gap,acc]):(y[i.gap,acc])
      s.idx = s.idx[-c(1, length(s.idx))]
      pos.idx[[acc]] = s.idx
      
      s = paste0(genomes.all[[acc]][s.idx], collapse = '')
      if(x[i.gap,acc] < y[i.gap,acc]){
        seqs = c(seqs, s)  
      } else {
        seqs = c(seqs, as.character(Biostrings::complement(DNAString(s))))  
      }
    }
    
    # message('Alignment')
    names(seqs) = names(pos.idx)
    
    aln.seqs[[i.gap]] = seqs
    aln.pos[[i.gap]] = pos.idx
  }
  cat('\n')
  
  
  file.aln.data = paste(path.work, "alignment_data_chr_", i.chr, '.Rdata', sep = '')
  save(list=c("aln.seqs", "aln.pos", "d", "x", "y", "idx.gap", "d.max", "d.zero", "idx.single", "idx.aln"),
       file = file.aln.data, compress = F)
  
  rm(v)
  rm(v.block)
  rm(v.gap)
  # rm(v.next)
  # rm(v.prev)
  rm(aln.seqs)
  rm(aln.pos)
}