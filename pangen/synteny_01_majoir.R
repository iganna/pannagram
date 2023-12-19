#library(ggplot2)

suppressMessages({
  library(Biostrings)
  library(seqinr)
  #library(spgs)  # reverseComplement("actg")
  library(foreach)
  library(doParallel)
  library(stringr)
  library(optparse)
})


source("utils/utils.R")
source("pangen/synteny_funcs.R")

pokazStage('Step 4. Alignment-1. Remaining syntenic (major) matches')

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-q", "--path.query"), type="character", default=NULL, 
              help="path to query chomosome fasta files", metavar="character"),  
  make_option(c("-i", "--path.blast"), type="character", default=NULL, 
              help="path to blast results", metavar="character"),
  make_option(c("-o", "--path.aln"), type="character", default=NULL, 
              help="path to the output directory with alignments", metavar="character"),
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
              help="alignment of all chromosomes vs all or not: T/F", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# Set the number of cores for parallel processing
num_cores <- ifelse(!is.null(opt$cores), opt$cores, 30)
myCluster <- makeCluster(num_cores, type = "PSOCK")
registerDoParallel(myCluster)

path.query    <- ifelse(!is.null(opt$path.query), opt$path.query, path.query)
path.blast.res <- ifelse(!is.null(opt$path.blast), opt$path.blast, path.blast.res)
path.aln      <- ifelse(!is.null(opt$path.aln), opt$path.aln, path.aln)
path.base     <- ifelse(!is.null(opt$path.ref), opt$path.ref, path.base)
base.acc      <- ifelse(!is.null(opt$pref), opt$pref, base.acc)
path.gaps     <- ifelse(!is.null(opt$path.gaps), opt$path.gaps, path.gaps)

n.chr.ref     <- ifelse(!is.null(opt$n.chr.ref), opt$n.chr.ref, n.chr.ref)
n.chr.acc     <- ifelse(!is.null(opt$n.chr.acc), opt$n.chr.acc, n.chr.acc)
all.vs.all    <- ifelse(!is.null(opt$all.vs.all), as.logical(opt$all.vs.all), all.vs.all)


# Create folders for the alignment results
if(!dir.exists(path.aln)) dir.create(path.aln)
if(!dir.exists(path.gaps)) dir.create(path.gaps)


#' ============================================================================

max.len = 10^6
len.blast = 50

files.blast <- list.files(path.blast.res, pattern = "\\.txt$")


pattern <- "(.+)_(\\d+)_(\\d+)\\.txt"

# Применение регулярного выражения
extracted_data <- lapply(files.blast, function(file) {
  matches <- regmatches(file, regexpr(pattern, file))
  if (length(matches) > 0) {
    parts <- strsplit(matches, "_")[[1]]
    return(list(NAME = parts[1], X = as.numeric(parts[2]), Y = as.numeric(parts[3])))
  } else {
    return(NULL)
  }
})



# for.flag = F
# tmp = foreach(f.blast = files.blast, .packages=c('crayon','stringr','Biostrings', 'seqinr'), .verbose = F)  %dopar% {  # which accession to use
for.flag = T
for(f.blast in files.blast){
  
  # Remove the '.txt' extension
  pref.comb <- sub("\\.txt$", "", f.blast)
  
  # Parser for BLAST-resulr file
  parts <- strsplit(pref.comb, "_")[[1]]
  
  base.chr <- parts[length(parts) - 1]
  query.chr <- parts[length(parts)]
  parts = parts[-c(length(parts) - 1, length(parts))]
  acc <- paste0(parts, collapse = '_')
  
  pokaz(acc, query.chr, base.chr)

  # file.aln.postgap3 <- paste(path.aln, paste0(pref.comb,  '_postgap3.rds', collapse = ''), sep = '')
  # if(file.exists(file.aln.postgap3)) {
  #   if(for.flag) next
  #   return(NULL)
  # }  
  
  # Read reference sequences
  base.file = paste0(base.acc, '_chr', base.chr , '.', 'fasta', collapse = '')
  pokaz('Base:', base.file)
  base.fas.fw = readFastaMy(paste(path.base, base.file, sep = ''))
  base.fas.fw = seq2nt(base.fas.fw)
  base.fas.bw = revCompl(base.fas.fw)
  base.len = length(base.fas.bw)

  # Read query sequences
  query.file = paste(acc, '_chr',query.chr, '.fasta', sep = '')
  pokaz('Query:', query.file)

  query.fas.chr = readFastaMy(paste(path.query, query.file, sep = ''))
  query.fas.chr = seq2nt(query.fas.chr)
  query.len = length(query.fas.chr)
  
  # Output files
  file.aln.pre <- paste(path.aln, paste0(pref.comb, '_maj.rds', collapse = ''), sep = '')
  # file.maj.idx <- paste(path.aln, paste0(pref.comb, '_maj_idx.rds', collapse = ''), sep = '')
  # file.raw.idx <- paste(path.aln, paste0(pref.comb, '_raw_idx.rds', collapse = ''), sep = '')
  
  # if(T){   
  if(!file.exists(file.aln.pre)){
    
    pokaz('Alignment:', acc, query.chr, base.chr)
    
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
    
    
    # x.dir = setDir(x, base.len = base.len)
    # checkCorrespToGenome(x.dir, query.fas = query.fas.chr,
    #                      base.fas.fw = base.fas.fw,
    #                      base.fas.bw = base.fas.bw)
    
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
    if(!isSorted(x.major$p.beg)) pokazAttention('1!!')
    
    # Remove wrong blocks completely
    remain.block = x.block$block.id
    if(length(remain.block) > 0){
      x.major = x.major[x.major$block.id %in% remain.block,]
    }
    
    # Check remained blocks
    x.major$id = rank(x.major$id)
    x.major$block = c(1, abs(x.major$id[-1] - x.major$id[-nrow(x.major)]) != 1) 
    x.major$block = x.major$block + c(1, abs(x.major$dir[-1] != x.major$dir[-nrow(x.major)])) 
    x.major$block = (x.major$block > 0) * 1
    
    x.major$block.id = cumsum(x.major$block)  # block ID
    
    if(length(unique((x.major$dir))) > 1){  # different dirertions exist
      cnt = table(x.major$block.id, x.major$dir)
      print(cnt)
      if(sum(cnt[,1] * cnt[,2]) != 0) stop('Blocks in x.major are wrongly defined')
    }
    
    # Save - Raw without sequences
    # saveRDS(x[, !(colnames(x) %in% c('V8', 'V9'))], file.raw.idx, compress = F)
    
    # ---- Filtration ----
    x = x[x.major$idx.maj,]
    x$block.id = x.major$block.id
    # saveRDS(x, file.aln.pre, compress = F)
    
    
    # ---- Remove short overlaps: twice, because from "both sides" ----
    for(i.tmp in 1:2){
      x = cleanBigOverlaps(x)
      x = cutSmallOverlaps(x)
    }
    
    # ---- Check sequences after cuts ----
    x.dir = setDir(x, base.len = base.len)
    checkCorrespToGenome(x.dir, query.fas = query.fas.chr,
                         base.fas.fw = base.fas.fw,
                         base.fas.bw = base.fas.bw)
    
    if(!isSorted(x$p.beg)) pokazAttention('6!!')
    
    # ---- Check uniquness of occupancy ----
    pos.q.occup = rep(0, base.len)
    for(irow in 1:nrow(x)){
      # pos.q.occup[x$V2[irow]:x$V3[irow]] = pos.q.occup[x$V2[irow]:x$V3[irow]] + 1
      # if(sum(pos.q.occup[x$V4[irow]:x$V5[irow]]) > 0) stop()
      pos.q.occup[x$V4[irow]:x$V5[irow]] = pos.q.occup[x$V4[irow]:x$V5[irow]] + 1
    }
    if(sum(pos.q.occup > 1) > 0) stop('Overlaps in base are remained')
    pokaz('Occupancy of base', sum(pos.q.occup))
    
    # Save
    saveRDS(x, file.aln.pre, compress = F)
    # saveRDS(x.major, file.maj.idx, compress = F)
    
  } else {
    pokaz('Reading instead of analysis', file.aln.pre)
    x = readRDS(file.aln.pre)
  } # if blast alignment exists

  if((nrow(x) <= 1) || (is.null(x))) {
    pokaz('No gaps')
    
    rmSafe(x)
    rmSafe(x.major)
    rmSafe(base.fas.bw)
    rmSafe(base.fas.fw)
    rmSafe(query.fas.chr)
    gc()
    
    if(for.flag) next
    return(NULL)
  }
  
  
  # ---- Get gaps ----
  pokaz('Get gaps')

  x.dir = setDir(x, base.len = base.len)
  checkCorrespToGenome(x.dir, query.fas = query.fas.chr, 
                       base.fas.fw = base.fas.fw, 
                       base.fas.bw = base.fas.bw)
  
  
  # Find occupied positions
  pos.q.free = rep(0, query.len)  # free positions in query
  pos.b.free = rep(0, base.len)  # free positions in base
  for(irow in 1:nrow(x)){
    pos.q.free[x$V2[irow]:x$V3[irow]] <- pos.q.free[x$V2[irow]:x$V3[irow]] + 1
    pos.b.free[x$V4[irow]:x$V5[irow]] <- pos.b.free[x$V4[irow]:x$V5[irow]] + 1
  }
  if((sum(pos.q.free > 1) != 0)) stop('Coverage query is wrong')
  if((sum(pos.b.free > 1) != 0)) stop('Coverage base is wrong')
  
  pos.q.free = -pos.q.free
  pos.b.free = -pos.b.free
  
  # ---- Write gaps ----
  # Within non-occupied positions find those, which can be

  pref.comarisson = paste('acc_', acc, '_qchr_', query.chr, '_bchr_', base.chr, '_', sep = '')
  # Query-file
  file.gap.query = paste0(path.gaps, pref.comarisson, 'query.fasta', collapse = '')
  # Base file
  file.gap.base = paste0(path.gaps, pref.comarisson, 'base.fasta', collapse = '')
  pokaz('Create gaps for', file.gap.query)
  pokaz('Create gaps for', file.gap.base)
  
  for(irow in 1:(nrow(x)-1)){
    
    # Common file name
    pref.gap = paste('gap_', irow, '_', irow + 1, '_', sep = '')
    
    if(x$bl[irow] != x$bl[irow+1]) next
    d1 = x$V2[irow+1] - x$V3[irow]
    
    pos = sort(c(x$V4[irow], x$V5[irow], x$V4[irow+1], x$V5[irow+1]))
    d2 = pos[3] - pos[2]
    if(d1 * d2 == 0) next  # Glue Zero, but not necessary
    
    # Short gaps will be covered later!!!
    # if((d1 <= len.blast) & (d2 <= len.blast)){
    #   x.tmp = x[c(irow,irow+1),]
    #   x.tmp = setDir(x.tmp, base.len = base.len)
    #   x.tmp = glueByThreshold(x.tmp, 100, query.fas = query.fas.chr,
    #                       base.fas.fw = base.fas.fw,
    #                       base.fas.bw = base.fas.bw, file.log=NULL)
    #   x.tmp = getBase(x.tmp, base.len = base.len)
    #   x[irow+1,] = x.tmp
    #   idx.remove = c(idx.remove, irow)
    # }
    
    # Don't consider for blast
    if(d1 <= len.blast){
      pos.q.free[(x$V3[irow]+1):(x$V2[irow+1]-1)] = -1
    }
    if(d2 <= len.blast){
      pos.b.free[(x$V5[irow]+1):(x$V4[irow+1]-1)] = -1
    }
    
    if(!((d1 >= len.blast) & (d2 >= len.blast))) next
    
    # Create files for BLAST
    pos.gap.q = x$V3[irow]:x$V2[irow + 1]
    pos = sort(c(x$V4[irow], x$V5[irow], x$V4[irow+1], x$V5[irow+1]))
    pos.gap.b = pos[2]:pos[3]
    pos.gap.q = pos.gap.q[-c(1, length(pos.gap.q))]
    pos.gap.b = pos.gap.b[-c(1, length(pos.gap.b))]
    
    # Save occupancy
    pos.q.free[pos.gap.q] = irow
    pos.b.free[pos.gap.b] = irow
    
    # if(file.exists(file.gap.query)) next
    
    if(abs(pos.gap.q[1] - pos.gap.q[length(pos.gap.q)]) > max.len) next
    if(abs(pos.gap.b[1] - pos.gap.b[length(pos.gap.b)]) > max.len) next
    
    # ---- Write query ----
    # Define Chunks
    s.q = query.fas.chr[pos.gap.q]
    s.q = nt2seq(s.q)
    n.bl = 500
    len.s.q = nchar(s.q)
    if(len.s.q > n.bl){
      p.beg = seq(1, len.s.q, n.bl)
      p.end = seq(n.bl, len.s.q, n.bl)
    } else {
      p.beg = 1
      p.end = len.s.q
    }
    s.q = splitSeq(s.q, n = n.bl)
    
    # Standsrd naming (as before)
    pref.q = paste(pref.comarisson, pref.gap,
                    'query', '|', pos.gap.q[p.beg], '|', pos.gap.q[p.end], sep = '')
    
    if(length(s.q) != length(pref.q)) stop('Chunk lengths do not much')
    names(s.q) = pref.q
    writeFastaMy(s.q, file.gap.query, append = T)
    
    # ---- Write base ----
    s.b = base.fas.fw[pos.gap.b]
    s.b = nt2seq(s.b)
    
    s.base.names = paste(pref.comarisson, pref.gap,
                         'base', '|', pos.gap.b[1], '|', pos.gap.b[length(pos.gap.b)], sep = '')
    
    names(s.b) = s.base.names
      
    writeFastaMy(s.b, file.gap.base, append = T)
    
  }  # irow search for gaps
  
  # ---- Write remained blocks ----
  pokaz('Create fasta for the remained sequences')
  file.gap.query = paste0(path.gaps, pref.comarisson, 'residual_query.fasta', collapse = '')
  # Base file
  file.gap.base = paste0(path.gaps, pref.comarisson, 'residual_base.fasta', collapse = '')
  pokaz('Create gaps for', file.gap.query)
  pokaz('Create gaps for', file.gap.base)
  
  
  ## ---- Write query ----
  # Query: Zero-coverage blocks
  diffs = diff(c(1, (pos.q.free != 0) * 1, 1))
  begs = which(diffs == -1)
  ends = which(diffs == 1) - 1
  for(irow in 1:length(begs)){
    if((ends[irow] - begs[irow]) < len.blast) next
    if((ends[irow] - begs[irow]) > max.len) next
    pos.gap.q = begs[irow]:ends[irow]
    
    # Define Chunks
    s.q = query.fas.chr[pos.gap.q]
    s.q = nt2seq(s.q)
    n.bl = 500
    len.s.q = nchar(s.q)
    if(len.s.q > n.bl){
      p.beg = seq(1, len.s.q, n.bl)
      p.end = seq(n.bl, len.s.q, n.bl)
    } else {
      p.beg = 1
      p.end = len.s.q
    }
    s.q = splitSeq(s.q, n = n.bl)
    
    # Standard naming (as before)
    pref.q = paste(pref.comarisson,
                   'resid_query', '|', pos.gap.q[p.beg], '|', pos.gap.q[p.end], sep = '')
    
    if(length(s.q) != length(pref.q)) stop('Chunk lengths do not much')
    names(s.q) = pref.q
    writeFastaMy(s.q, file.gap.query, append = T)
  }
  
  
  ## ---- Write base ----
  # Query: Zero-coverage blocks
  diffs = diff(c(1, (pos.b.free != 0) * 1, 1))
  begs = which(diffs == -1)
  ends = which(diffs == 1) - 1
  
  if(length(begs) > 0){
    for(irow in 1:length(begs)){
      if((ends[irow] - begs[irow]) < len.blast) next
      if((ends[irow] - begs[irow]) > max.len) next
      
      pos.gap.b = begs[irow]:ends[irow]
      
      s.b = base.fas.fw[pos.gap.b]
      s.b = nt2seq(s.b)
      
      s.base.names = paste(pref.comarisson,
                           'resid_base', '|', pos.gap.b[1], '|', pos.gap.b[length(pos.gap.b)], sep = '')
      
      names(s.b) = s.base.names
      
      writeFastaMy(s.b, file.gap.base, append = T)
    }
  }
  
  
  
  rmSafe(x)
  rmSafe(x.major)
  rmSafe(base.fas.bw)
  rmSafe(base.fas.fw)
  rmSafe(query.fas.chr)
  gc()
  
}  # combinations




