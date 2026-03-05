# Get SV positions, GFF files, dencity files and consensys sequences
# Find SVs and create GFF file

suppressMessages({ library(Biostrings)
  library(rhdf5)
  library(foreach)
  library(doParallel)
  library(optparse)
  library(pannagram)
  library(crayon)
})

# ***********************************************************************
# ---- Alignment types ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) 

# ***********************************************************************
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--ref"),  type = "character", default = NULL, help = "prefix of the reference file"),
  make_option("--path.features.msa", type = "character", default = NULL, help = "Path to msa dir (features)"),
  make_option("--path.seq", type = "character", default = NULL, help = "Path to seq dir"),
  make_option("--path.sv", type = "character", default = NULL, help = "Path to sv dir"),
  make_option("--path.gff", type = "character", default = NULL, help = "Path to gff dir"),
  make_option(c("--aln.type"),  type = "character", default = aln.type.msa, help = "type of alignment ('pan', 'ref', etc)"),
  make_option(c("--acc.anal"),  type = "character", default = NULL, help = "files with accessions to analyze"),
  make_option(c("--stat.only"), type = "character", default = NULL, help = "files with accessions to analyze"),
  make_option("--cutoff", type = "numeric", default = 0.9, help = "Frequency cutoff"),
  make_option("--min.len", type = "integer", default = 15, help = ""),
  make_option("--big.len", type = "integer", default = 50, help = ""),
  make_option("--max.len", type = "integer", default = 30000, help = ""),
  make_option("--cores",             type = "integer",   default = 1,            help = "number of cores to use for parallel processing"),
  make_option("--path.log",         type = "character", default = NULL, help = "Path for log files")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);


# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

path.log <- opt$path.log
if (!is.null(path.log) & !is.null(log.level)) {
  if (!dir.exists(path.log)){
    dir.create(path.log)
  }
}

# ***********************************************************************
# ---- Modes ----

# If only the statistics is needed
if (!is.null(opt$stat.only)) {
  flag.stat.only = T
} else {
  flag.stat.only = F
}

# Accessions to analyse
acc.anal <- opt$acc.anal
if(acc.anal == 'NULL') acc.anal = NULL
if(!is.null(acc.anal)){
  if (!file.exists(acc.anal)) {
    acc.anal = NULL
    pokazAttention('File', acc.anal, 'does NOT exists, so no accession filtration is applied.')
  } else {
    tmp = read.table(acc.anal, stringsAsFactors = F)
    acc.anal = tmp[,1]
  }
}

accessions = acc.anal
# pokaz('Accessions', acc.anal)

# ***********************************************************************
# ---- Paths ----
path.features.msa <- opt$path.features.msa
if(!dir.exists(path.features.msa)) stop(paste0('No Consensus directory found!', path.features.msa))

path.seq <- opt$path.seq
if(!dir.exists(path.seq)) stop(paste0('No Consensus directory found!', path.seq))

path.sv <- opt$path.sv
if(!dir.exists(path.sv)) stop(paste0('No SV directory found!', path.features.msa))

path.gff <- opt$path.gff
if(!dir.exists(path.gff)) stop(paste0('No GFF directory found!', path.features.msa))

# ***********************************************************************
# ---- Variables ----

num.cores <- opt$cores 
pokaz('Number of cores', num.cores)


cutoff <- opt$cutoff
min.len <- opt$min.len
big.len <- opt$big.len
max.len <- opt$max.len

# ***********************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

# Alignment prefix
if (!is.null(opt$aln.type)) {
  aln.type = opt$aln.type
} else {
  aln.type = aln.type.msa
}

# Reference genome
ref.name <- opt$ref
if(ref.name == "NULL" || is.null(ref.name)) ref.name <- ''

# Common code for aln.pref, ref.suffix and s.combinations
source(system.file("utils/chunk_combinations.R", package = "pannagram")) 

# ***********************************************************************
# ---- Positions of SVs ----

# sv.pos.all = c()
# sv.beg.all = c()
# sv.end.all = c()

sv.pos.list <- list()
sv.beg.list <- list()
sv.end.list <- list()

for(s.comb in s.combinations){
  
  
  file.sv.pos.log = paste0(path.log, 'sv.pos_', s.comb, '.log')
  file.sv.pos.rds = paste0(path.log, 'sv.pos_', s.comb, '.rds')
  file.sv.beg.log = paste0(path.log, 'sv.beg_', s.comb, '.log')
  file.sv.beg.rds = paste0(path.log, 'sv.beg_', s.comb, '.rds')
  file.sv.end.log = paste0(path.log, 'sv.end_', s.comb, '.log')
  file.sv.end.rds = paste0(path.log, 'sv.end_', s.comb, '.rds')
  
  
  if(!file.exists(file.sv.pos.log)) invisible(file.create(file.sv.pos.log))
  if(!file.exists(file.sv.beg.log)) invisible(file.create(file.sv.beg.log))
  if(!file.exists(file.sv.end.log)) invisible(file.create(file.sv.end.log))
  
  if(checkDone(file.sv.pos.log) &&
     checkDone(file.sv.beg.log) &&
     checkDone(file.sv.end.log)){
    
    pokaz('Reading beg-end', s.comb)
    sv.pos = readRDS(file.sv.pos.rds)
    sv.beg = readRDS(file.sv.beg.rds)
    sv.end = readRDS(file.sv.end.rds)
    
    k <- length(sv.pos.list) + 1
    sv.pos.list[[k]] <- sv.pos
    sv.beg.list[[k]] <- sv.beg
    sv.end.list[[k]] <- sv.end
    next
  } 
  
  
  q.chr = strsplit(s.comb, '_')[[1]][1]
  r.chr = strsplit(s.comb, '_')[[1]][2]
  if(q.chr != r.chr) {
    pokazAttention('Alignment of query chromosome', q.chr, 
                   'to the reference chromosome', r.chr, 'will be skipped.',
                   'For SV calling, only matching chromosome IDs are allowed.')
    next
  }
  
  pokaz('Run for combination', s.comb, "...")
  print(Sys.time())
  
  # Get file for the combination
  file.comb = paste0(path.features.msa, aln.pref, s.comb, ref.suff,'.h5')
  if(!file.exists(file.comb)) stop('Alignment file does not exist')
  # pokaz('Alignment file', file.comb)
  
  # Get accessions
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  
  # If the set of accessions to analyze has been provided, then focus only on them.
  if(!is.null(acc.anal)){
    accessions = intersect(accessions, acc.anal)
    if(length(setdiff(acc.anal, accessions)) > 0) {
      pokazAttention('Not all the accessions of interest are in the Alignment.') 
    }
  }
  
  n.acc = length(accessions)
  
  pokaz('Combine SV info from accessions...')
  print(Sys.time())
  
  # sv.cover = 0
  # for(acc in accessions){
  
  file.sv.cover = paste0(path.log, 'sv_cover_',s.comb,'.rds')
  file.sv.cover.log = paste0(path.log, 'sv_cover_',s.comb,'.log')
  if(!file.exists(file.sv.cover.log)) invisible(file.create(file.sv.cover.log))
  
  if(checkDone(file.sv.cover.log)){
    sv.cover = readRDS(file.sv.cover)
  } else {
    
    myCluster <- parallel::makeCluster(num.cores, type = "PSOCK")
    doParallel::registerDoParallel(myCluster)
    
    foreach::foreach(
      acc = accessions,
      .inorder = FALSE,
      .packages = c("rhdf5", "pannagram")
    ) %dopar% {
      
      
      pokaz('Positions of accession', acc)
      print(Sys.time())
      
      file.log.loop = paste0(path.log, 'log_', acc, '_', s.comb, '.log')
      if(!file.exists(file.log.loop)) invisible(file.create(file.log.loop))
      
      file.loop.save = paste0(path.log, acc, '_', s.comb, '.rds')
      
      # ---- Check log Done ----
      if(checkDone(file.log.loop)){
        return(NULL)
      }
      
      v <- h5read(file.comb, paste0(gr.accs.e, acc))
      v[is.na(v)] <- 0
      
      pokaz('Find gaps..')
      print(Sys.time())
      
      sv.acc <- findOnes(v == 0)
      if (nrow(sv.acc) == 0) {
        out <- 0
        rhdf5::H5close()
        return(out)
      }
      
      if (sv.acc$beg[1] == 1) sv.acc <- sv.acc[-1, , drop = FALSE]
      if (sv.acc$end[nrow(sv.acc)] == length(v)) sv.acc <- sv.acc[-nrow(sv.acc), , drop = FALSE]
      
      v.r <- rank(abs(v)) * sign(v)
      
      pokaz('Exclude gaps around inversions..')
      print(Sys.time())
      
      # sv.acc$beg.r = v.r[sv.acc$beg - 1]
      # sv.acc$end.r = v.r[sv.acc$end + 1]
      # sv.acc.na = sv.acc[(sv.acc$end.r - sv.acc$beg.r) != 1,]
      # if(nrow(sv.acc.na) > 0){
      #   for(irow in 1:nrow(sv.acc.na)){
      #     v[sv.acc.na$beg[irow]:sv.acc.na$end[irow]] = Inf
      #   }  
      # }
      
      idx.bad <- which((v.r[sv.acc$end + 1] - v.r[sv.acc$beg - 1]) != 1)
      # if (length(idx.bad) > 0) {
      #   for (i in idx.bad) {
      #     v[sv.acc$beg[i] : sv.acc$end[i]] <- Inf
      #   }
      # }
      if (length(idx.bad) > 0) {
        b <- sv.acc$beg[idx.bad]
        e <- sv.acc$end[idx.bad]
        
        # Difference array trick: mark interval coverage in O(L + k)
        d <- integer(length(v) + 1L)
        d[b] <- d[b] + 1L
        d[e + 1L] <- d[e + 1L] - 1L
        
        v[cumsum(d[-(length(v) + 1L)]) > 0L] <- Inf
      }
      
      out <- as.integer(v == 0)
      
      # Close HDF5 handles inside worker
      rhdf5::H5close()
      
      saveRDS(out, file.loop.save)
      pokaz('Done.', file=file.log.loop, echo=T)
      
      rm(v, v.r, out, sv.acc, idx.bad, b, e, d)
      gc()
      
      return(NULL)
    }
    
    parallel::stopCluster(myCluster)
    
    pokaz('Combine sv.cover')
    sv.cover <- 0
    for (acc in accessions) {
      pokaz('Acc', acc)
      file.loop.save = paste0(path.log, acc, '_', s.comb, '.rds')
      
      if (!file.exists(file.loop.save)) stop('Problem')
      
      x <- readRDS(file.loop.save)
      sv.cover <- sv.cover + x
    }
    
    
    saveRDS(sv.cover, file.sv.cover)
    pokaz('Done.', file=file.sv.cover.log, echo=T)
  }
  
  
  # SV groups
  pokaz('Create SV groups...')
  print(Sys.time())
  
  sv.pos = findOnes((sv.cover != 0)*1)
  if(nrow(sv.pos) == 0){
    pokazAttention('SVs were not generaed, and IT IS OK!')
    next
  } 
  sv.pos$len = abs(sv.pos$beg - sv.pos$end) + 1 # do not change
  sv.pos$beg = sv.pos$beg - 1
  sv.pos$end = sv.pos$end + 1
  sv.pos$beg[sv.pos$beg == 0] = 1
  sv.pos$end[sv.pos$end > length(sv.cover)] = length(sv.cover)
  
  # sv.beg = c()
  # sv.end = c()
  # for(acc in accessions){
  #   pokaz('Positions of accession', acc)
  #   print(Sys.time())
  #   
  #   v = h5read(file.comb, paste0(gr.accs.e, acc))
  #   v[is.na(v)] = 0
  #   sv.beg = cbind(sv.beg, v[sv.pos$beg])
  #   sv.end = cbind(sv.end, v[sv.pos$end])
  # }
  # colnames(sv.beg) = accessions
  # colnames(sv.end) = accessions
  
  n.sv = nrow(sv.pos)
  n.acc = length(accessions)
  
  sv.beg <- matrix(0, nrow = n.sv, ncol = n.acc,
    dimnames = list(NULL, accessions)
  )
  
  sv.end <- matrix(0, nrow = n.sv, ncol = n.acc,
    dimnames = list(NULL, accessions)
  )
  
  # for (acc in accessions) {
  #   pokaz('Positions of accession', acc)
  #   print(Sys.time())
  #   
  #   v <- h5read(file.comb, paste0(gr.accs.e, acc))
  #   v[is.na(v)] <- 0
  #   
  #   sv.beg[, acc] <- v[sv.pos$beg]
  #   sv.end[, acc] <- v[sv.pos$end]
  # }
  
  pokaz('Parallel starts')
  # The same code as above but paralleled
  myCluster <- parallel::makeCluster(num.cores, type = "PSOCK")
  doParallel::registerDoParallel(myCluster)
  
  res <- foreach(acc = accessions, 
                 .packages = "rhdf5"
                 ) %dopar% {
    
    v <- h5read(file.comb, paste0(gr.accs.e, acc))
    v[is.na(v)] <- 0
    
    list(
      acc = acc,
      beg = v[sv.pos$beg],
      end = v[sv.pos$end]
    )
  }
  parallel::stopCluster(myCluster)
  
  # Combine results together
  pokaz('Combine results together after parallel')
  accs <- vapply(res, `[[`, "", "acc")
  sv.beg <- do.call(cbind, lapply(res, `[[`, "beg"))
  sv.end <- do.call(cbind, lapply(res, `[[`, "end"))
  colnames(sv.beg) <- accs
  colnames(sv.end) <- accs
  
  # save(list = ls(), file = "tmp_workspace_sv.RData")
  
  # Check ranks
  pokaz('Check ranks...')
  print(Sys.time())
  
  for(acc in accessions){
    acc.r = c(sv.beg[,acc] + 0.1, sv.end[,acc] - 0.1)
    acc.r = rank(abs(acc.r)) * sign(acc.r)
    acc.r = matrix(acc.r, ncol = 2)
    
    d = acc.r[,2] - acc.r[,1]
    idx = which(d != 1)
    if(length(idx) > 0){
      sv.beg[idx, acc] = 0
      sv.end[idx, acc] = 0
    }
  }
  
  
  # Clean up mismatches
  for(acc in accessions){
    sv.beg[sv.end[,acc] == 0, acc] = 0
    sv.end[sv.beg[,acc] == 0, acc] = 0
    
    idx = which(sv.beg[, acc] * sv.end[,acc] < 0)
    if(length(idx) > 0){
      sv.beg[idx, acc] = 0
      sv.end[idx, acc] = 0
    }
  }
  
  # Clean up empty
  pokaz('Clean up empty...')
  print(Sys.time())
  
  sv.na = (rowSums(sv.beg) == 0) | (rowSums(sv.end) == 0)
  sv.pos = sv.pos[!sv.na,]
  sv.beg = sv.beg[!sv.na,]
  sv.end = sv.end[!sv.na,]
  
  # Calculate lengths
  sv.len.acc = abs(sv.beg - sv.end) - 1  # do not change
  colnames(sv.len.acc) = accessions
  
  
  # Remove those, which length in more that the length of SV
  idx = rowMax(sv.len.acc) <= sv.pos$len
  sv.pos = sv.pos[idx,]
  sv.beg = sv.beg[idx,]
  sv.end = sv.end[idx,]
  sv.len.acc = sv.len.acc[idx,]
  
  if(nrow(sv.pos) == 0){
    pokazAttention('SVs were not generaed, and IT IS OK!')
    next
  }
  
  # Calculate frequencies
  pokaz('Calculate frequencies')
  sv.pos$freq.min = 0
  sv.pos$freq.max = 0
  for(i.col in 1:n.acc){
    sv.pos$freq.min = sv.pos$freq.min + 1*((sv.len.acc[,i.col] <= (1-cutoff) * sv.pos$len) & (!is.na(sv.len.acc[,i.col])))
    sv.pos$freq.max = sv.pos$freq.max + 1*((sv.len.acc[,i.col] >= cutoff * sv.pos$len) & (!is.na(sv.len.acc[,i.col])))
  }
  
  sv.pos$freq.sum = sv.pos$freq.min + sv.pos$freq.max 
  sv.pos$single = (sv.pos$freq.sum == n.acc) * 1
  sv.pos = cbind(sv.pos, sv.len.acc[,1:n.acc])
  
  # save(list = ls(), file = "tmp_workspace_sv.RData")
  
  # Clean up
  pokaz('Additional Clean up..')
  idx = !((sv.pos$freq.min == 0) & (sv.pos$freq.sum == n.acc))
  sv.pos = sv.pos[idx,]
  sv.beg = sv.beg[idx,]
  sv.end = sv.end[idx,]
  sv.len.acc = sv.len.acc[idx,]
  
  idx = !((sv.pos$freq.max == 0) & (sv.pos$freq.sum == n.acc))
  sv.pos = sv.pos[idx,]
  sv.beg = sv.beg[idx,]
  sv.end = sv.end[idx,]
  sv.len.acc = sv.len.acc[idx,]
  
  if(sum((sv.pos$freq.min == 0) & (sv.pos$freq.sum == n.acc)) > 0) stop('WRONG1')
  if(sum((sv.pos$freq.max == 0) & (sv.pos$freq.sum == n.acc)) > 0) stop('WRONG2')
  
  s.gr.len <- nchar(as.character(nrow(sv.pos)))
  i.chr = strsplit(s.comb, '_')[[1]][1]
  gr = paste('SVgr', i.chr, 'id', sprintf("%0*d", s.gr.len, 1:nrow(sv.pos)), sep = '_')
  sv.pos = cbind(gr, sv.pos)
  sv.beg = cbind(gr, as.data.frame(sv.beg))
  sv.end = cbind(gr, as.data.frame(sv.end))
  sv.pos$chr = i.chr
  sv.beg$chr = i.chr
  sv.end$chr = i.chr
  
  # sv.pos.all = rbind(sv.pos.all, sv.pos)
  # sv.beg.all = rbind(sv.beg.all, sv.beg)
  # sv.end.all = rbind(sv.end.all, sv.end)
  
  if (nrow(sv.pos) == 0) next
  pokaz('Save posisiotns')
  saveRDS(sv.pos, file.sv.pos.rds)
  saveRDS(sv.beg, file.sv.beg.rds)
  saveRDS(sv.end, file.sv.end.rds)
  
  pokaz('Done.', file=file.sv.pos.log, echo=T)
  pokaz('Done.', file=file.sv.beg.log, echo=T)
  pokaz('Done.', file=file.sv.end.log, echo=T)
  
  k <- length(sv.pos.list) + 1
  sv.pos.list[[k]] <- sv.pos
  sv.beg.list[[k]] <- sv.beg
  sv.end.list[[k]] <- sv.end
  
  H5close()
  gc()
  
}

sv.pos.all <- do.call(rbind, sv.pos.list)
sv.beg.all <- do.call(rbind, sv.beg.list)
sv.end.all <- do.call(rbind, sv.end.list)

sv.pos.all$name = paste0(sv.pos.all$gr, '|', sv.pos.all$len)
rownames(sv.pos.all) = sv.pos.all$name

pokaz('Save pos-file...')
file.sv.pos = paste0(path.sv, 'sv_pangen_pos.rds')
saveRDS(sv.pos.all, file.sv.pos)

if(flag.stat.only){
  pokaz('Stat was generated')
  quit(save="no")
} 

## ---- Stop for Stat ----
pokaz('Save beg-end files...')
file.sv.pos.beg = paste0(path.sv, 'sv_pangen_beg.rds')
file.sv.pos.end = paste0(path.sv, 'sv_pangen_end.rds')
saveRDS(sv.beg.all, file.sv.pos.beg)
saveRDS(sv.end.all, file.sv.pos.end)


# ---- FASTA of seSVs ----
pokaz('Generate Fasta...')

file.sv.small =  paste0(path.sv, 'seq_sv_short.fasta')
file.sv.big =  paste0(path.sv, 'seq_sv_large.fasta')

seqs.small = c()
seqs.big = c()
for(s.comb in s.combinations){
  i.chr = comb2ref(s.comb)
  pokaz('Chromosome', i.chr)
  print(Sys.time())
  
  file.chr = paste0(path.seq, 'seq_cons_', s.comb, ref.suff ,'.fasta')
  if(!file.exists(file.chr)) stop(paste0('File with the consensus sequence does not exist:', file.chr))
  s.chr = readFasta(file.chr)
  s.chr = seq2nt(s.chr)
  
  # Small sequences  
  idx.small = which((sv.pos.all$single == 1) & 
                      (sv.pos.all$len >= min.len) & 
                      (sv.pos.all$len < big.len) & (sv.pos.all$chr == i.chr))
  for(irow in idx.small){
    s.tmp = s.chr[(sv.pos.all$beg[irow] + 1):(sv.pos.all$end[irow] - 1) ]
    if(sum(s.tmp == 'N') > (0.5 * length(s.tmp))) next
    seqs.small[sv.pos.all$name[irow]] = nt2seq(s.tmp)
  }
  
  # Big sequence
  idx.big = which((sv.pos.all$single == 1) & 
                    (sv.pos.all$len >= big.len) &
                    (sv.pos.all$len < max.len) & (sv.pos.all$chr == i.chr))
  # print(head(idx.big))
  for(irow in idx.big){
    s.tmp = s.chr[(sv.pos.all$beg[irow] + 1):(sv.pos.all$end[irow] - 1) ]
    if(sum(s.tmp == 'N') > (0.5 * length(s.tmp))) next
    seqs.big[sv.pos.all$name[irow]] = nt2seq(s.tmp)
  }
  pokaz('Number of large sequences', length(seqs.big))
  pokaz('Number of short sequences', length(seqs.small))
}

pokaz('Save Fasta...')
writeFasta(seqs.small, file.sv.small)
writeFasta(seqs.big, file.sv.big)


# ---- GFF files ----
pokaz('Gff files for accessions..')
print(Sys.time())

sv.version = 6
file.sv.gff = paste(path.gff, 'svs_pangen_v',sprintf("%02d", sv.version),'.gff', sep = '')

sv.pos.all$V10 = 1:nrow(sv.pos.all)
rownames(sv.pos.all) = sv.pos.all$gr

# save(list = ls(), file = 'tmx_workspace_sv.RData')
# stop('Enough')

## ---- Single-event ----
pokaz('Single-event..')
print(Sys.time())

sv.se = sv.pos.all[sv.pos.all$single == 1,]
sv.se.type = rep('indel', nrow(sv.se))

# percent.phasing = 0.11 # Was used for 27-genome paper
percent.phasing = 0.10
n.acc = length(accessions)
val.insert = c(1:n.acc)[(1:n.acc) <= (percent.phasing * n.acc)]
val.delet = c(1:n.acc)[(1:n.acc) >= (n.acc - percent.phasing * n.acc)]

sv.se.type[sv.se$freq.max %in% val.insert] = 'insertion'
sv.se.type[sv.se$freq.max %in% val.delet] = 'deletion'

sv.annot = paste('ID=', sv.se$gr,  
                 # ';te=', sv.se$te, 
                 ';presence=',sv.se$freq.max,sep = '',
                 ';len_init=', sv.se$len)

# save(list = ls(), file = "tmp_workspace_sv.RData")

sv.se.gff = data.frame(V1 = paste0('PanGen_Chr', sv.se$chr),
                       V2 = 'pannagram',
                       V3 = sv.se.type, 
                       V4 = sv.se$beg + 1, 
                       V5 = sv.se$end - 1,
                       V6 = '.', V7 = '+', V8 = '.', 
                       V9 = sv.annot, 
                       V10 = sv.pos.all[sv.se$gr, 'V10'])

## ---- Multiple-event ----
pokaz('Multiple-event..')
print(Sys.time())

sv.me = sv.pos.all[sv.pos.all$single != 1,, drop=F]

if(nrow(sv.me) > 0){
  s.multi = 'complex'
  sv.me.gff = data.frame(V1 = paste0('PanGen_Chr', sv.me$chr),
                         V2 = 'pannagram',
                         V3 = s.multi, V4 = sv.me$beg, V5 = sv.me$end,
                         V6 = '.', V7 = '+', V8 = '.', 
                         V9 = paste0('ID=', sv.me$gr, ';len_init=', sv.me$len), 
                         V10 = sv.pos.all[sv.me$gr, 'V10'])
  
  sv.gff = rbind(sv.se.gff, sv.me.gff)  
} else {
  sv.gff = sv.se.gff
}

sv.gff = sv.gff[order(sv.gff$V10),]

options(scipen = 999)
# write.table(sv.gff[,1:9], file.sv.gff, quote = F, row.names = F, col.names = F, sep = '\t')
writeGFF(sv.gff[,1:9], file.sv.gff)
options(scipen = 0)

# ---- GFF In accessions ----
pokaz('Gff files for accessions comb..')
print(Sys.time())

for(i.acc in 1:length(accessions)){
  acc = accessions[i.acc]
  pokaz('Generate GFF for accession', acc)
  print(Sys.time())
  
  df = sv.gff
  df$V1 = paste0(acc, '_Chr', sv.pos.all$chr)
  df$V4 = sv.beg.all[,acc] + 1
  df$V5 = sv.end.all[,acc] - 1
  
  df$V9 = paste('ID=', sv.pos.all$gr, '.', acc, 
                ';len_init=', sv.pos.all$len,
                ';len_acc=', abs(sv.end.all[,acc]-sv.beg.all[,acc])-1, sep = '')
  
  df = df[df$V5 > df$V4,]
  
  df = df[sv.pos.all$len > 0,]
  
  df = df[!is.na(df$V4),]
  df = df[!is.na(df$V5),]
  
  df = df[df$V4 != 0,]
  df = df[df$V5 != 0,]
  
  if(sum(df$V5 < df$V4) > 0){
    # save(list = ls(), file = 'tmx_workspace_sv_acc.RData')
    pokazAttention("sum(df$V5 < df$V4) = ", sum(df$V5 < df$V4))
    stop()
  }
  
  # Strand
  idx.strand = df$V4 < 0
  df$V7[idx.strand] = '-'
  tmp = abs(df$V4[idx.strand])
  df$V4[idx.strand] = abs(df$V5[idx.strand])
  df$V5[idx.strand] = tmp
  
  file.sv.gff = paste(path.gff, 'svs_acc_', acc, '_v',sprintf("%02d", sv.version),'.gff', sep = '')
  
  options(scipen = 999)
  # write.table(df[,1:9], file.sv.gff, quote = F, row.names = F, col.names = F, sep = '\t')
  writeGFF(df[,1:9], file.sv.gff)
  options(scipen = 0)
}



file.sv.small =  paste0(path.sv, 'seq_sv_short.fasta')
file.sv.big =  paste0(path.sv, 'seq_sv_large.fasta')


seqs.small = c()
seqs.big = c()
for(s.comb in s.combinations){
  i.chr = comb2ref(s.comb)
  pokaz('Chromosome', i.chr)
  print(Sys.time())
  
  file.chr = paste0(path.seq, 'seq_cons_', s.comb, ref.suff ,'.fasta')
  if(!file.exists(file.chr)) stop(paste0('File with the consensus sequence does not exist:', file.chr))
  s.chr = readFasta(file.chr)
  s.chr = seq2nt(s.chr)
  
  # Small sequences  
  idx.small = which((sv.pos.all$single == 1) & 
                      (sv.pos.all$len >= min.len) & 
                      (sv.pos.all$len < big.len) & (sv.pos.all$chr == i.chr))
  for(irow in idx.small){
    s.tmp = s.chr[(sv.pos.all$beg[irow] + 1):(sv.pos.all$end[irow] - 1) ]
    if(sum(s.tmp == 'N') > (0.5 * length(s.tmp))) next
    seqs.small[paste(sv.pos.all$gr[irow],sv.pos.all$len[irow], sep = '|')] = paste0(s.tmp, collapse = '')
  }
  
  # Big sequence
  idx.big = which((sv.pos.all$single == 1) & 
                    (sv.pos.all$len >= big.len) &
                    (sv.pos.all$len < max.len) & (sv.pos.all$chr == i.chr))
  # print(head(idx.big))
  for(irow in idx.big){
    s.tmp = s.chr[(sv.pos.all$beg[irow] + 1):(sv.pos.all$end[irow] - 1) ]
    # print(sum(s.tmp == 'N'))
    # print((0.5 * length(s.tmp)))
    if(sum(s.tmp == 'N') > (0.5 * length(s.tmp))) next
    seqs.big[paste(sv.pos.all$gr[irow],sv.pos.all$len[irow], sep = '|')] = paste0(s.tmp, collapse = '')
  }
  pokaz('Number of large sequences', length(seqs.big))
  pokaz('Number of short sequences', length(seqs.small))
}

writeFasta(seqs.small, file.sv.small)
writeFasta(seqs.big, file.sv.big)



# ---- GFF files ----
pokaz('Gff files for accessions..')
print(Sys.time())

sv.version = 6
file.sv.gff = paste(path.gff, 'svs_pangen_v',sprintf("%02d", sv.version),'.gff', sep = '')

sv.pos.all$V10 = 1:nrow(sv.pos.all)
rownames(sv.pos.all) = sv.pos.all$gr

# save(list = ls(), file = 'tmx_workspace_sv.RData')
# stop('Enough')

## ---- Single-event ----
pokaz('Single-event..')
print(Sys.time())

sv.se = sv.pos.all[sv.pos.all$single == 1,]
sv.se.type = rep('indel', nrow(sv.se))

# percent.phasing = 0.11 # Was used for 27-genome paper
percent.phasing = 0.10
n.acc = length(accessions)
val.insert = c(1:n.acc)[(1:n.acc) <= (percent.phasing * n.acc)]
val.delet = c(1:n.acc)[(1:n.acc) >= (n.acc - percent.phasing * n.acc)]

sv.se.type[sv.se$freq.max %in% val.insert] = 'insertion'
sv.se.type[sv.se$freq.max %in% val.delet] = 'deletion'

sv.annot = paste('ID=', sv.se$gr,  
                 # ';te=', sv.se$te, 
                 ';presence=',sv.se$freq.max,sep = '',
                 ';len_init=', sv.se$len)

# save(list = ls(), file = "tmp_workspace_sv.RData")

sv.se.gff = data.frame(V1 = paste0('PanGen_Chr', sv.se$chr),
                       V2 = 'pannagram',
                       V3 = sv.se.type, 
                       V4 = sv.se$beg + 1, 
                       V5 = sv.se$end - 1,
                       V6 = '.', V7 = '+', V8 = '.', 
                       V9 = sv.annot, 
                       V10 = sv.pos.all[sv.se$gr, 'V10'])

## ---- Multiple-event ----
pokaz('Multiple-event..')
print(Sys.time())

sv.me = sv.pos.all[sv.pos.all$single != 1,, drop=F]

if(nrow(sv.me) > 0){
  s.multi = 'complex'
  sv.me.gff = data.frame(V1 = paste0('PanGen_Chr', sv.me$chr),
                         V2 = 'pannagram',
                         V3 = s.multi, V4 = sv.me$beg, V5 = sv.me$end,
                         V6 = '.', V7 = '+', V8 = '.', 
                         V9 = paste0('ID=', sv.me$gr, ';len_init=', sv.me$len), 
                         V10 = sv.pos.all[sv.me$gr, 'V10'])
  
  sv.gff = rbind(sv.se.gff, sv.me.gff)  
} else {
  sv.gff = sv.se.gff
}

sv.gff = sv.gff[order(sv.gff$V10),]

options(scipen = 999)
# write.table(sv.gff[,1:9], file.sv.gff, quote = F, row.names = F, col.names = F, sep = '\t')
writeGFF(sv.gff[,1:9], file.sv.gff)
options(scipen = 0)

# ---- GFF In accessions ----
pokaz('Gff files for accessions comb..')
print(Sys.time())

for(i.acc in 1:length(accessions)){
  acc = accessions[i.acc]
  pokaz('Generate GFF for accession', acc)
  print(Sys.time())
  
  df = sv.gff
  df$V1 = paste0(acc, '_Chr', sv.pos.all$chr)
  df$V4 = sv.beg.all[,acc] + 1
  df$V5 = sv.end.all[,acc] - 1
  
  df$V9 = paste('ID=', sv.pos.all$gr, '.', acc, 
                ';len_init=', sv.pos.all$len,
                ';len_acc=', abs(sv.end.all[,acc]-sv.beg.all[,acc])-1, sep = '')
  
  df = df[df$V5 > df$V4,]
  
  df = df[sv.pos.all$len > 0,]
  
  df = df[!is.na(df$V4),]
  df = df[!is.na(df$V5),]
  
  df = df[df$V4 != 0,]
  df = df[df$V5 != 0,]
  
  if(sum(df$V5 < df$V4) > 0){
    # save(list = ls(), file = 'tmx_workspace_sv_acc.RData')
    pokazAttention("sum(df$V5 < df$V4) = ", sum(df$V5 < df$V4))
    stop()
  }
  
  # Strand
  idx.strand = df$V4 < 0
  df$V7[idx.strand] = '-'
  tmp = abs(df$V4[idx.strand])
  df$V4[idx.strand] = abs(df$V5[idx.strand])
  df$V5[idx.strand] = tmp
  
  file.sv.gff = paste(path.gff, 'svs_acc_', acc, '_v',sprintf("%02d", sv.version),'.gff', sep = '')
  
  options(scipen = 999)
  # write.table(df[,1:9], file.sv.gff, quote = F, row.names = F, col.names = F, sep = '\t')
  writeGFF(df[,1:9], file.sv.gff)
  options(scipen = 0)
}

