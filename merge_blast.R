#library(ggplot2)
suppressMessages(library(Biostrings))
suppressMessages(library('seqinr'))
suppressMessages(library('spgs'))  # reverseComplement("actg")
suppressMessages(library('foreach'))
suppressMessages(library(doParallel))
suppressMessages(library(stringr))
suppressMessages(library("optparse"))
source("synteny_infer.R")


message('=====================================================')
message('|  Alignment1: after first blast of parts           |')
message('-----------------------------------------------------')


args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-q", "--path.query"), type="character", default=NULL, 
              help="path to query chomosome fasta files", metavar="character"),  
  make_option(c("-i", "--path.blast"), type="character", default=NULL, 
              help="path to blast results", metavar="character"),
  make_option(c("-o", "--path.aln"), type="character", default=NULL, 
              help="path to the output directory with alignments", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="type of fasta files", metavar="character"),
  make_option(c("-p", "--pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("-g", "--path.gaps"), type="character", default=NULL, 
              help="prefix of the directory with gaps", metavar="character"),
  make_option(c("-r", "--path.ref"), type="character", default=NULL, 
              help="path to the reference file", metavar="character"),
  make_option(c("-n", "--n.chr.ref"), type="character", default=NULL, 
              help="number of chromosomes in the reference genome", metavar="character"),
  make_option(c("-m", "--n.chr.acc"), type="character", default=NULL, 
              help="number of chromosomes in the accessions", metavar="character"),
  make_option(c("-a", "--all.vs.all"), type="character", default=NULL, 
              help="alignment of all chromosomes vs all or not: T/F", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

if (!is.null(opt$path.query)) path.query <- opt$path.query
if (!is.null(opt$path.blast)) path.blast.res <- opt$path.blast
if (!is.null(opt$path.aln)) path.aln <- opt$path.aln
if (!is.null(opt$path.ref)) path.base <- opt$path.ref
if (!is.null(opt$type)) base.suff <- opt$type
if (!is.null(opt$pref)) base.acc <- opt$pref
if (!is.null(opt$path.gaps)) path.gaps <- opt$path.gaps

if (!is.null(opt$n.chr.ref)) n.chr.ref <- opt$n.chr.ref
if (!is.null(opt$n.chr.acc)) n.chr.acc <- opt$n.chr.acc
if (!is.null(opt$all.vs.all)) all.vs.all <- as.logical(opt$all.vs.all)


if(!dir.exists(path.aln)) dir.create(path.aln)
if(!dir.exists(path.gaps)) dir.create(path.gaps)

# ============================================================================

myCluster <- makeCluster(30, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

# ============================================================================

max.len = 10^6

file.log <- NULL

files.query = list.files(path = path.query, pattern = "\\.fasta$")
query.name = gsub("*.fasta","",files.query)

query.name = unique(sapply(query.name, function(s){ strsplit(s, '_')[[1]][1] }))
names(query.name) <- NULL

message(paste0(c('Accessions:', query.name), collapse = " "))
message(paste('Base accession:', base.acc))


# Combinations of chromosomes query-base to cheate the alignments
chromosome.pairs = c()
for(i.query in 1:length(query.name)){
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
# print(chromosome.pairs)


for.flag = F
tmp = foreach(i.chr.pair = 1:nrow(chromosome.pairs), .packages=c('stringr','Biostrings', 'seqinr', 'spgs'))  %dopar% {  # which accession to use
  
  #for.flag = T  
  #for(i.chr.pair in 1:nrow(chromosome.pairs)){
  
  
  
  i.query = chromosome.pairs[i.chr.pair, 1]
  
  if(query.name[i.query] == base.acc){
    if(for.flag){
      next
    } else {
      return(NULL)
    }
  }  
  
  query.chr = chromosome.pairs[i.chr.pair, 2]
  base.chr = chromosome.pairs[i.chr.pair, 3]
  
  pref.comb = paste0(query.name[i.query], '_', query.chr, '_', base.chr, collapse = '')
  
  # If the blast-result is not there -> next
  t.file <- paste(path.blast.res, pref.comb, '.txt', sep = '')
  # message(t.file)
  if(!file.exists(t.file)) {
    if(for.flag){
      next
    } else {
      return(NULL)
    }
  }  
  
  file.aln.postgap3 <- paste(path.aln, paste0(pref.comb,  '_postgap3.rds', collapse = ''), sep = '')
  if(file.exists(file.aln.postgap3)) {
    if(for.flag){
      next
    } else {
      return(NULL)
    }
  }  
  
  
  base.file = paste0(base.acc, '_chr', base.chr , '.', base.suff, collapse = '')
  print(paste('Base:', base.file, sep=' '))
  
  # Read fasta files
  base.fas.fw = read.fasta(paste(path.base, base.file, sep = ''))
  base.fas.fw = base.fas.fw[[1]] # do not remove this 1, file contains only one sequence
  base.fas.bw = reverseComplement(base.fas.fw)
  base.len = length(base.fas.fw)
  
  query.file = paste(query.name[i.query], '_chr',query.chr, '.fasta', sep = '')
  print(paste('Query:', query.file, sep=' '))
  query.fas = read.fasta(paste(path.query, query.file, sep = ''))
  query.fas.chr = query.fas[[1]] # do not remove this 1, file contains only one sequence
  query.len = length(query.fas.chr)
  
  
  file.aln.pre <- paste(path.aln, paste0(pref.comb, '_pre7.rds', collapse = ''), sep = '')
  if(!file.exists(file.aln.pre)){    
    
    message(paste('Alignment:', query.name[i.query], query.chr, base.chr, sep = ' '))
    
    # ---- Read blast results ----
    t = read.table(t.file, stringsAsFactors = F, header = F)
    t = t[(t$V3-t$V2) > 100,]
    t = removeOverlapdAfterBlast(t)
    
    message(paste('Read blast results finished, numer of rows is', nrow(t)))
    
    # ---- Get right positions + dir ----
    start.pos = as.numeric(sapply(strsplit(t[,1], "\\|"), "[", 4)) - 1
    t[,2:3] = t[,2:3] + start.pos
    t <- setDir(t, base.len)
    
    message('Glue Zero..')
    t <- glueZero(t)
    
    message('Glue Second..')
    t <- glueAlmostZero(t, threshold = 2000,
                        query.fas.chr = query.fas.chr,
                        base.fas.fw = base.fas.fw,
                        base.fas.bw = base.fas.bw)
    
    
    t = removeCompleteOverlaps2(t, base.len)
    t = removeShortOverlaps2(t, base.len)
    t <- glueZero(t)
    t <- glueAlmostZero(t, threshold = 2000,
                        query.fas.chr = query.fas.chr,
                        base.fas.fw = base.fas.fw,
                        base.fas.bw = base.fas.bw)
    
    checkCorrespToGenome(t, query.fas = query.fas.chr,
                         base.fas.fw = base.fas.fw,
                         base.fas.bw = base.fas.bw)
    
    
    t <- cleanJumps2(t, base.len)
    t <- glueAlmostZero(t, threshold = 2000,
                        query.fas.chr = query.fas.chr,
                        base.fas.fw = base.fas.fw,
                        base.fas.bw = base.fas.bw)
    
    
    
    thresholds = c(100, 500, 1000, 1500)
    t = glueByThreshold(t, thresholds, query.fas = query.fas.chr,
                        base.fas.fw = base.fas.fw,
                        base.fas.bw = base.fas.bw, file.log=NULL)
    
    
    t = removeCompleteOverlaps2(t, base.len)
    t = removeShortOverlaps2(t, base.len)
    t <- glueZero(t)
    t <- glueAlmostZero(t, threshold = 2000,
                        query.fas.chr = query.fas.chr,
                        base.fas.fw = base.fas.fw,
                        base.fas.bw = base.fas.bw)
    
    t <- cleanJumps2(t, base.len)
    print('Remove Complete and short')
    
  } else {
    print('reading')
    t = readRDS(file.aln.pre)
  } # if blast alignment exists
  
  if((nrow(t) <= 1) || (is.null(t))) {
    message('no gaps')
    
    if(for.flag){
      next
    } else {
      return(NULL)
    }
    
  }
  
  
  # ---- Get gaps ----
  print('Get gaps')
  t = t[t[,'V7'] > 10000,,drop=F]
  if((nrow(t) <= 1) || (is.null(t))) {
    message('no gaps')
    if(for.flag){
      next
    } else {
      return(NULL)
    }
  }
  
  t = cleanJumps2(t, base.len, min.len = 30000)
  saveRDS(object = t, file = file.aln.pre) 
  
  query.len <- length(query.fas.chr)
  
  pos.q.free = rep(0, query.len)  # free positions if query
  pos.b.free = rep(0, base.len)  # free positions if base
  
  t = getBase(t, base.len)
  for(irow in 1:nrow(t)){
    pos.q.free[t$V2[irow]:t$V3[irow]] <- 1
    pos.b.free[t$V4[irow]:t$V5[irow]] <- 1
    # print(sum(pos.q.free))
    # print(sum(pos.b.free))
  }
  
  print(sum(pos.q.free))
  print(sum(pos.b.free))
  
  idx.gaps.ok = c()
  for(irow in 1:(nrow(t) - 1)){
    pos.gap.q = t$V3[irow]:t$V2[irow + 1] 
    pos = sort(c(t$V4[irow], t$V5[irow], t$V4[irow+1], t$V5[irow+1]))
    pos.gap.b = pos[2]:pos[3]
    
    if(length(pos.gap.q) == 2) next
    if(length(pos.gap.b) == 2) next
    
    pos.gap.q = pos.gap.q[-c(1, length(pos.gap.q))]
    pos.gap.b = pos.gap.b[-c(1, length(pos.gap.b))]
    
    n.q.occ = sum(pos.q.free[pos.gap.q])
    n.b.occ = sum(pos.b.free[pos.gap.b])
    if(!((n.q.occ == 0) && (n.b.occ == 0))) next

    i.gap = irow
    file.gap.query = paste0(path.gaps, query.name[i.query], '_', query.chr, '_', base.chr, '_query_',i.gap,'.txt', collapse = '')
    file.gap.base = paste0(path.gaps, query.name[i.query], '_', query.chr, '_', base.chr, '_base_',i.gap,'.txt', collapse = '')
    if(file.exists(file.gap.query)) next
      
    
    p.beg = min(pos.gap.q)
    p.end = max(pos.gap.q)
    if(p.end - p.beg > max.len) next
    
    # Split region into pieces of 500nt for more effective BLAST
    if(p.end-p.beg > 500){      
      p.tmp <- seq(p.beg, p.end, 500)
      
      # print(c(p.beg, p.end))
      # print(p.tmp)
      
      p.beg = p.tmp
      p.end = p.tmp-1
      p.beg = p.beg[-length(p.beg)]
      p.end = p.end[-1]
      p.end[length(p.end)] = max(pos.gap.q)
    }
    
    s.query = c()
    s.query.names = c()
    for(j.gap in 1:length(p.beg)){
      # print(c(p.beg[j.gap], p.end[j.gap], length(query.fas.chr)))
      s.query = c(s.query, paste0(query.fas.chr[p.beg[j.gap]:p.end[j.gap]], 
                                  collapse = ''))
      s.query.names = c(s.query.names, 
                        paste0('gap', 
                               '_', i.gap, '_', i.gap + 1,
                               '|',query.name[i.query] , '|', query.chr, '|',
                               'query', '|', p.beg[j.gap], '|', p.end[j.gap], collapse = ''))
    }
    
    names(s.query) <- s.query.names      
    # Write query
    # print(file.gap.query)
    writeFasta(s = s.query, file.gap.query)
    
    s.base= c()
    s.base.names = c()

    p.beg = min(pos.gap.b)
    p.end = max(pos.gap.b)
    if(p.end - p.beg > max.len) next
    
    s.base = paste0(base.fas.fw[p.beg:p.end], collapse = '')
    s.base.names = paste0('gap', 
                          '_', i.gap, '_', i.gap + 1,
                            '|',query.name[i.query] , '|', query.chr, '|',
                            'base','|', p.beg, '|', p.end, collapse = '')

    names(s.base) <- s.base.names
    # Write base
    writeFasta(s = s.base, file.gap.base)
    
  }
}  # combinations
