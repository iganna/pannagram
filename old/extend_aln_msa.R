# This script first runs MSA for short indels
# then creates files for the maffs alignment

suppressMessages(library('foreach'))
suppressMessages(library(doParallel))
suppressMessages(library(msa))
suppressMessages(library(dplyr))
suppressMessages(library("optparse"))


myCluster <- makeCluster(30, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

message('=====================================================')
message('  Prepare files for msa  ')
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


len.short = 50
n.flank = 30

for(i.chr in 1:5){
  
  message(i.chr)
  
  file.aln.data = paste(path.work, "alignment_data_chr_", i.chr, '.Rdata', sep = '')
  load(file.aln.data)
  message('Data loaded')
  
  
  file.aln.data = paste(path.work, "alignment_data_chr_", i.chr, '_post.Rdata', sep = '')
  if(!file.exists(file.aln.data)){
    
    idx.short = setdiff(which(d.max <= len.short), idx.single)
    
    for.flag = F
    res.msa <- foreach(seqs = aln.seqs[idx.short], pos.idx = aln.pos[idx.short], 
                       idx.gap.pos = idx.gap[idx.short,1],
                       .packages=c('muscle', 'Biostrings'))  %dopar% {
                         # for.flag = T
                         # for(i.gap in idx.short){
                         
                         # seqs = aln.seqs[[i.gap]]
                         # pos.idx = aln.pos[[i.gap]]
                         # idx.gap.pos = dx.gap[i.gap,1]
                         
                         names(seqs) = names(pos.idx)
                         set = DNAStringSet(seqs)
                         aln = muscle(set, quiet = T)
                         
                         set = as.character(aln)
                         n.pos = nchar(set[1])
                         val.acc.pos = matrix(0, nrow=n.pos, ncol=n.acc)
                         colnames(val.acc.pos) = accessions
                         for(acc in names(set)){
                           s = strsplit(set[[acc]],'')[[1]]
                           val.acc.pos[s != '-',acc] = pos.idx[[acc]]
                         }
                         val.acc.pos = cbind(val.acc.pos, idx.gap.pos + (1:n.pos) / (n.pos+1))
                         
                         # if(for.flag){
                         #   res.msa[[i.gap]] = val.acc.pos
                         # } else {
                         #   return(val.acc.pos)
                         # }
                         #
                         return(val.acc.pos)
                       }
    
    
    save(list=c("res.msa"),
         file = file.aln.data, compress = F)
    rm(res.msa)
  }
  
  # ------------------------------------------------------------------------
  # MAFFT multiple alignment
  message('MAFFT sequences')
  path.fasta = paste(path.work, 'mafft_fasta/', sep='')
  if(!dir.exists(path.fasta)) system(paste('mkdir ', path.fasta, sep = ''))
  
  idx.long = which(d.max > len.short)
  idx.long = setdiff(idx.long, idx.single)
  
  res.msa <- foreach(seqs = aln.seqs[idx.long], pos.idx = aln.pos[idx.long], idx.gap.pos = idx.gap[idx.long,1], .packages=c('muscle', 'Biostrings'))  %dopar% {
    
    f.pref = paste(path.fasta, 'aln',
                   '_pos_', idx.gap.pos, '_chr_', i.chr, sep = '')
    fasta.f.flank = paste(f.pref,'_flank_', n.flank,'.fasta', sep = '')
    
    names(seqs) = names(pos.idx)
    
    for(acc in names(seqs)){  # not all accessions!!
      
      write(paste('>acc', acc, pos.idx[[acc]][1], pos.idx[[acc]][length(pos.idx[[acc]])],
                  'flank',n.flank, sep = '_'),
            file = fasta.f.flank, append = T)
      write(paste0(c(rep('A', n.flank), tolower(seqs[[acc]]), rep('C', n.flank)), collapse =''),
            file = fasta.f.flank, append = T)
    } 
  }
}