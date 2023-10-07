#library(ggplot2)
suppressMessages(library(Biostrings))
suppressMessages(library('seqinr'))
suppressMessages(library('spgs'))  # reverseComplement("actg")
suppressMessages(library('foreach'))
suppressMessages(library(doParallel))
suppressMessages(library(stringr))
suppressMessages(library("optparse"))
source("utils.R")
source("synteny_funcs.R")

pokazStage('Alignment1: after first blast of parts')


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


# for.flag = F
# tmp = foreach(i.chr.pair = 1:nrow(chromosome.pairs), .packages=c('stringr','Biostrings', 'seqinr', 'spgs'))  %dopar% {  # which accession to use
  
for.flag = T
for(i.chr.pair in 1:nrow(chromosome.pairs)){
  
  
  
  i.query = chromosome.pairs[i.chr.pair, 1]
  
  if(query.name[i.query] == base.acc){
    if(for.flag) next
    return(NULL)
  }  
  
  query.chr = chromosome.pairs[i.chr.pair, 2]
  base.chr = chromosome.pairs[i.chr.pair, 3]
  
  pref.comb = paste0(query.name[i.query], '_', query.chr, '_', base.chr, collapse = '')
  
  # If the blast-result is not there -> next
  t.file <- paste(path.blast.res, pref.comb, '.txt', sep = '')
  # message(t.file)
  if(!file.exists(t.file)) {
    if(for.flag) next
    return(NULL)
  }  
  
  file.aln.postgap3 <- paste(path.aln, paste0(pref.comb,  '_postgap3.rds', collapse = ''), sep = '')
  if(file.exists(file.aln.postgap3)) {
    if(for.flag) next
    return(NULL)
  }  
  
  # # Read reference sequences
  # base.file = paste0(base.acc, '_chr', base.chr , '.', base.suff, collapse = '')
  # pokaz('Base:', base.file)
  # 
  # base.fas.fw = readFastaMy(paste(path.base, base.file, sep = ''))
  # base.fas.fw = seq2nt(base.fas.fw)
  # base.fas.bw = revCompl(base.fas.fw)
  # base.len = length(base.fas.bw)
  # 
  # # Read query sequences
  # query.file = paste(query.name[i.query], '_chr',query.chr, '.fasta', sep = '')
  # pokaz('Query:', query.file)
  # 
  # query.fas.chr = readFastaMy(paste(path.query, query.file, sep = ''))
  # query.fas.chr = seq2nt(query.fas.chr)
  # query.len = length(query.fas.chr)
  
  # Output files
  file.aln.pre <- paste(path.aln, paste0(pref.comb, '_maj.rds', collapse = ''), sep = '')
  file.maj.idx <- paste(path.aln, paste0(pref.comb, '_maj_idx.rds', collapse = ''), sep = '')
  file.raw.idx <- paste(path.aln, paste0(pref.comb, '_raw_idx.rds', collapse = ''), sep = '')
  
  if(!file.exists(file.aln.pre)){    
    
    pokaz('Alignment:', query.name[i.query], query.chr, base.chr)
    
    # ---- Read blast results ----
    x = read.table(t.file, stringsAsFactors = F, header = F)
    
    pokaz('Read blast results finished, numer of rows is', nrow(x))
    
    ## ---- Pre-processing ----
    # Save true base coordinate
    x$p.beg <- ifelse(x$V4 < x$V5, x$V4, x$V5)
    x$p.end <- ifelse(x$V4 < x$V5, x$V5, x$V4)
    
    # Set correct position
    start.pos = as.numeric(sapply(strsplit(x[,1], "\\|"), "[", 4)) - 1
    x[,2:3] = x[,2:3] + start.pos
    
    # Set direction
    x$dir = (x$V4 > x$V5) * 1
    
    # Set past id
    x$part.id = cumsum(c(T, diff(as.numeric(as.factor(x$V1))) != 0))
    
    # ---- Major skeleton ----
    
    idx.maj = which(!c(F, x$V1[-1] == x$V1[-nrow(x)]))
    x.major = cbind(x[idx.maj,setdiff(colnames(x), c('V8', 'V9'))], idx.maj)
    
    # Sort to set the id
    x.major = x.major[order(x.major$V2),]  # not needed
    x.major$id = 1:nrow(x.major)
    
    # Sort by the position in the reference
    x.major = x.major[order(-x.major$V7),]
    x.major = x.major[order(x.major$p.beg),]
    
    # If complete overlap - remove shortest
    idx.overlap = which(x.major$p.end[-1] <= x.major$p.end[-nrow(x.major)]) + 1
    if(length(idx.overlap) > 0){
      x.major = x.major[-idx.overlap,]  
    }
    
    # Remove those, which are not in the correct place
    x.major$id = rank(x.major$id)
    idx = which(abs(diff(x.major$id)) > 1)
    idx = intersect(idx, idx + 1)
    
    if(length(idx) > 0){
      x.major = x.major[-idx,]
    }
    
    # ----  Define blocks  in the skeleton ----
    
    # re-arrange IDs
    x.major$id = rank(x.major$id)
    # [1 or 0] - [begin or end] of block
    x.major$block = c(1, abs(x.major$id[-1] - x.major$id[-nrow(x.major)]) != 1) 
    # block ID
    x.major$block.id = cumsum(x.major$block)
    
    # Analyse only beginnings of blocks
    x.block = x.major[x.major$block == 1,]
    x.block$id = rank(x.block$id)
    bl.len = abs(tapply(x.major$p.end, x.major$block.id, max) - tapply(x.major$p.beg, x.major$block.id, min))
    bl.len = bl.len[order(as.numeric(names(bl.len)))]
    x.block$len = bl.len
    
    # Remain blocks, only if they are in a correct place and long enough
    while(T){
      x.block$id = rank(x.block$id)
      idx = which(abs(diff(x.block$id)) > 1)
      idx = unique(c(idx, idx + 1))
      idx = idx[idx < nrow(x.block)]
      # idx = intersect(idx, idx +1)
      if(length(idx) == 0) break
      if(min(x.block$len[idx]) > 20000) break
      idx.remove = idx[x.block$len[idx] == min(x.block$len[idx])][1]
      pokaz('- Remove block', idx.remove, ';length:', x.block$len[idx.remove])
      x.block = x.block[-idx.remove,]
      # rownames(x.block) = NULL
    }
    
    # Remove wrong blocks completely
    remain.block = x.block$block.id
    if(length(remain.block) > 0){
      x.major = x.major[x.major$block.id %in% remain.block,]
    }
    
    file.aln.pre <- paste(path.aln, paste0(pref.comb, '_maj.rds', collapse = ''), sep = '')
    file.maj.idx <- paste(path.aln, paste0(pref.comb, '_maj_idx.rds', collapse = ''), sep = '')
    file.raw.idx <- paste(path.aln, paste0(pref.comb, '_raw_idx.rds', collapse = ''), sep = '')
    
    
    x.sk = x[x.major$idx.maj,]
    saveRDS(x.sk, file.aln.pre, compress = F)
    saveRDS(x.major, file.maj.idx, compress = F)
    saveRDS(x[, !(colnames(x) %in% c('V8', 'V9'))], file.raw.idx, compress = F)
    

  } else {
    print('reading')
    t = readRDS(file.aln.pre)
  } # if blast alignment exists
  
  # if((nrow(t) <= 1) || (is.null(t))) {
  #   message('no gaps')
  #   
  #   if(for.flag) next
  #   return(NULL)
  #   
  # }
  
  # 
  # # ---- Get gaps ----
  # print('Get gaps')
  # t = t[t[,'V7'] > 10000,,drop=F]
  # if((nrow(t) <= 1) || (is.null(t))) {
  #   message('no gaps')
  #   if(for.flag) next
  #   return(NULL)
  # }
  # 
  # t = cleanJumps2(t, base.len, min.len = 30000)
  # saveRDS(object = t, file = file.aln.pre) 
  # 
  # query.len <- length(query.fas.chr)
  # 
  # pos.q.free = rep(0, query.len)  # free positions if query
  # pos.b.free = rep(0, base.len)  # free positions if base
  # 
  # t = getBase(t, base.len)
  # for(irow in 1:nrow(t)){
  #   pos.q.free[t$V2[irow]:t$V3[irow]] <- 1
  #   pos.b.free[t$V4[irow]:t$V5[irow]] <- 1
  #   # print(sum(pos.q.free))
  #   # print(sum(pos.b.free))
  # }
  # 
  # print(sum(pos.q.free))
  # print(sum(pos.b.free))
  # 
  # idx.gaps.ok = c()
  # for(irow in 1:(nrow(t) - 1)){
  #   pos.gap.q = t$V3[irow]:t$V2[irow + 1] 
  #   pos = sort(c(t$V4[irow], t$V5[irow], t$V4[irow+1], t$V5[irow+1]))
  #   pos.gap.b = pos[2]:pos[3]
  #   
  #   if(length(pos.gap.q) == 2) next
  #   if(length(pos.gap.b) == 2) next
  #   
  #   pos.gap.q = pos.gap.q[-c(1, length(pos.gap.q))]
  #   pos.gap.b = pos.gap.b[-c(1, length(pos.gap.b))]
  #   
  #   n.q.occ = sum(pos.q.free[pos.gap.q])
  #   n.b.occ = sum(pos.b.free[pos.gap.b])
  #   if(!((n.q.occ == 0) && (n.b.occ == 0))) next
  # 
  #   i.gap = irow
  #   file.gap.query = paste0(path.gaps, query.name[i.query], '_', query.chr, '_', base.chr, '_query_',i.gap,'.txt', collapse = '')
  #   file.gap.base = paste0(path.gaps, query.name[i.query], '_', query.chr, '_', base.chr, '_base_',i.gap,'.txt', collapse = '')
  #   if(file.exists(file.gap.query)) next
  #     
  #   
  #   p.beg = min(pos.gap.q)
  #   p.end = max(pos.gap.q)
  #   if(p.end - p.beg > max.len) next
  #   
  #   # Split region into pieces of 500nt for more effective BLAST
  #   if(p.end-p.beg > 500){      
  #     p.tmp <- seq(p.beg, p.end, 500)
  #     
  #     # print(c(p.beg, p.end))
  #     # print(p.tmp)
  #     
  #     p.beg = p.tmp
  #     p.end = p.tmp-1
  #     p.beg = p.beg[-length(p.beg)]
  #     p.end = p.end[-1]
  #     p.end[length(p.end)] = max(pos.gap.q)
  #   }
  #   
  #   s.query = c()
  #   s.query.names = c()
  #   for(j.gap in 1:length(p.beg)){
  #     # print(c(p.beg[j.gap], p.end[j.gap], length(query.fas.chr)))
  #     s.query = c(s.query, paste0(query.fas.chr[p.beg[j.gap]:p.end[j.gap]], 
  #                                 collapse = ''))
  #     s.query.names = c(s.query.names, 
  #                       paste0('gap', 
  #                              '_', i.gap, '_', i.gap + 1,
  #                              '|',query.name[i.query] , '|', query.chr, '|',
  #                              'query', '|', p.beg[j.gap], '|', p.end[j.gap], collapse = ''))
  #   }
  #   
  #   names(s.query) <- s.query.names      
  #   # Write query
  #   # print(file.gap.query)
  #   writeFasta(s = s.query, file.gap.query)
  #   
  #   s.base= c()
  #   s.base.names = c()
  # 
  #   p.beg = min(pos.gap.b)
  #   p.end = max(pos.gap.b)
  #   if(p.end - p.beg > max.len) next
  #   
  #   s.base = paste0(base.fas.fw[p.beg:p.end], collapse = '')
  #   s.base.names = paste0('gap', 
  #                         '_', i.gap, '_', i.gap + 1,
  #                           '|',query.name[i.query] , '|', query.chr, '|',
  #                           'base','|', p.beg, '|', p.end, collapse = '')
  # 
  #   names(s.base) <- s.base.names
  #   # Write base
  #   writeFasta(s = s.base, file.gap.base)
  #   
  # }
}  # combinations
