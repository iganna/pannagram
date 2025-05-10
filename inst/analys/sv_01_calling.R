# Get SV positions, GFF files, dencity files and consensys sequences
# Find SVs and create GFF file

suppressMessages({ library(Biostrings)
  library(rhdf5)
  library('foreach')
  library(doParallel)
  library("optparse")
  library(pannagram)
  library(crayon)
})

pokazAttention('Make sure that the consensus sequences for the pangenome chromosomes have been generated.')

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--ref.pref"),  type = "character", default = NULL, help = "prefix of the reference file"),
  make_option(c("--path.cons"), type = "character", default = NULL, help = "path to directory with the consensus"),
  make_option(c("--cores"),     type = "integer",   default = 1, help = "number of cores to use for parallel processing"),
  make_option(c("--aln.type"),  type = "character", default = "default", help = "type of alignment ('msa_', 'comb_', 'v_', etc)"),
  make_option(c("--acc.anal"),  type = "character", default = NULL, help = "files with accessions to analyze"),
  make_option(c("--stat.only"), type = "character", default = NULL, help = "files with accessions to analyze")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);


# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

# ***********************************************************************

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

# Alignment prefix
if (!is.null(opt$aln.type)) {
  aln.type = opt$aln.type
} else {
  aln.type = aln.type.msa
}

# Reference genome
if (is.null(opt$ref.pref) || (opt$ref.pref == "NULL")) {
  ref.pref <- ""
  ref.suff = ""
} else {
  ref.pref <- opt$ref.pref
  ref.suff <- paste0('_ref_', ref.pref)
}

if (!is.null(opt$path.cons)) path.cons <- opt$path.cons

if(!dir.exists(path.cons)) stop(paste0('Consensus folder does nto exist', path.cons))

path.sv = paste0(path.cons, 'sv/')
if (!dir.exists(path.sv)) dir.create(path.sv)
if(!dir.exists(path.sv)) stop(paste0('SV folder does nto exist', path.cons))

path.gff = paste0(path.sv, 'gff/')
if (!dir.exists(path.gff)) dir.create(path.gff)
if(!dir.exists(path.gff)) stop(paste0('GFF folder does nto exist', path.cons))

cutoff = 0.90


# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type, ".*", ref.suff, "\\.h5")
# pokaz(s.pattern)
pokaz(path.cons)
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
# pokaz(files)

pref.combinations = gsub(aln.type, "", files)
pref.combinations <- sub(ref.suff, "", pref.combinations)
pref.combinations <- sub(".h5", "", pref.combinations)

if(length(pref.combinations) == 0){
  stop('No Combinations found.')
} else {
  pokaz('Combinations', pref.combinations)  
}

# ---- Positions of SVs ----

sv.pos.all = c()
sv.beg.all = c()
sv.end.all = c()

for(s.comb in pref.combinations){
  pokaz('* Combination', s.comb)
  
  # Get file for the combination
  file.comb = paste0(path.cons, aln.type, s.comb, ref.suff,'.h5')
  if(!file.exists(file.comb)) stop('Alignment file does not exist')
  pokaz('Alignment file', file.comb)
  
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
  sv.cover = 0
  for(acc in accessions){
    pokaz('Positions of accession', acc)
    v = h5read(file.comb, paste0(gr.accs.e, acc))
    v[is.na(v)] = 0
    
    pokaz('Find gaps..')
    sv.acc = findOnes(v == 0)
    if(nrow(sv.acc) == 0) next
    
    if(sv.acc$beg[1] == 1) sv.acc = sv.acc[-1,]
    if(sv.acc$end[nrow(sv.acc)] == length(v)) sv.acc = sv.acc[-nrow(sv.acc),]
    v.r = rank(abs(v)) * sign(v)
    
    pokaz('Exclude gaps around inversions..')
    sv.acc$beg.r = v.r[sv.acc$beg - 1]
    sv.acc$end.r = v.r[sv.acc$end + 1]
    sv.acc.na = sv.acc[(sv.acc$end.r - sv.acc$beg.r) != 1,]
    if(nrow(sv.acc.na) > 0){
      for(irow in 1:nrow(sv.acc.na)){
        v[sv.acc.na$beg[irow]:sv.acc.na$end[irow]] = Inf
      }  
    }
    sv.cover = sv.cover + (v == 0) * 1
  }
  
  # save(list = ls(), file = "tmp_workspace_sv.RData")
  # stop()
  
  # SV groups
  pokaz('Create SV groups...')
  
  sv.pos = findOnes((sv.cover != 0)*1)
  sv.pos$len = abs(sv.pos$beg - sv.pos$end) + 1 # do not change
  sv.pos$beg = sv.pos$beg - 1
  sv.pos$end = sv.pos$end + 1
  sv.pos$beg[sv.pos$beg == 0] = 1
  sv.pos$end[sv.pos$end > length(sv.cover)] = length(sv.cover)
  
  sv.beg = c()
  sv.end = c()
  for(acc in accessions){
    pokaz('Positions of accession', acc)
    v = h5read(file.comb, paste0(gr.accs.e, acc))
    v[is.na(v)] = 0
    sv.beg = cbind(sv.beg, v[sv.pos$beg])
    sv.end = cbind(sv.end, v[sv.pos$end])
  }
  colnames(sv.beg) = accessions
  colnames(sv.end) = accessions
  
  # save(list = ls(), file = "tmp_workspace_sv.RData")
  
  # Check ranks
  pokaz('Check ranks...')
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
  
  # Clean up empty
  pokaz('Clean up empty...')
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
  
  sv.pos.all = rbind(sv.pos.all, sv.pos)
  sv.beg.all = rbind(sv.beg.all, sv.beg)
  sv.end.all = rbind(sv.end.all, sv.end)
  
  H5close()
  gc()
  
}

pokaz('Saving....')

file.sv.pos = paste0(path.sv, 'sv_pangen_pos.rds')
saveRDS(sv.pos.all, file.sv.pos)

if(flag.stat.only){
  pokaz('Stat was generated')
  quit(save="no")
} 
## ---- Stop for Stat ----
pokaz('Stat....')
file.sv.pos.beg = paste0(path.sv, 'sv_pangen_beg.rds')
file.sv.pos.end = paste0(path.sv, 'sv_pangen_end.rds')
saveRDS(sv.beg.all, file.sv.pos.beg)
saveRDS(sv.end.all, file.sv.pos.end)



sv.mismatch = (sv.beg.all[, accessions] * sv.end.all[,accessions])  < 0
for(acc in accessions){
  sv.beg.all[sv.end.all[,acc] == 0, acc] = 0
  sv.end.all[sv.beg.all[,acc] == 0, acc] = 0
  
  idx = which(sv.beg.all[, acc] * sv.end.all[,acc] < 0)
  if(length(idx) > 0){
    sv.beg.all[idx, acc] = 0
    sv.end.all[idx, acc] = 0
  }
}

# ---- GFF files ----
pokaz('Gff files for accessions..')
sv.version = 6
file.sv.gff = paste(path.gff, 'svs_pangen_v',sprintf("%02d", sv.version),'.gff', sep = '')

save(list = ls(), file = "tmp_workspace_sv.RData")

sv.pos.all$V10 = 1:nrow(sv.pos.all)
rownames(sv.pos.all) = sv.pos.all$gr

# save(list = ls(), file = 'tmx_workspace_sv.RData')
# stop('Enough')

## ---- Single-event ----
pokaz('Single-event..')
sv.se = sv.pos.all[sv.pos.all$single == 1,]
sv.se.type = rep('indel', nrow(sv.se))

n.acc = length(accessions)
val.insert = c(1:n.acc)[(1:n.acc) <= (0.11 * n.acc)]
val.delet = c(1:n.acc)[(1:n.acc) >= (n.acc - 0.11 * n.acc)]

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
sv.me = sv.pos.all[sv.pos.all$single != 1,, drop=F]

if(nrow(sv.me) > 0){
  s.multi = 'multi'
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
write.table(sv.gff[,1:9], file.sv.gff, quote = F, row.names = F, col.names = F, sep = '\t')
options(scipen = 0)

# ---- GFF In accessions ----
pokaz('Gff files for accessions comb..')
for(i.acc in 1:length(accessions)){
  acc = accessions[i.acc]
  pokaz('Generate GFF for accession', acc)
  
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
    save(list = ls(), file = 'tmx_workspace_sv_acc.RData')
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
  write.table(df[,1:9], file.sv.gff, quote = F, row.names = F, col.names = F, sep = '\t')
  options(scipen = 0)
}

# ---- FASTA of seSVs ----
path.seq = paste0(path.cons, 'seq/')
if (!dir.exists(path.seq)) {
  pokazAttention('Consensus sequence doesn’t exist')
  stop('Please generate the pangenome consensus sequence first')
}

min.len = 15
big.len = 50
max.len = 30000


file.sv.small =  paste0(path.sv, 'seq_sv_small.fasta')
file.sv.big =  paste0(path.sv, 'seq_sv_big.fasta')

pokaz(file.sv.small)

seqs.small = c()
seqs.big = c()
for(s.comb in pref.combinations){
  i.chr = comb2ref(s.comb)
  pokaz('Chromosome', i.chr)
  file.chr = paste0(path.seq, 'seq_cons_', i.chr, '.fasta')
  if(!file.exists(file.chr)) stop(paste0('File with the consensus sequence does not exist:', file.chr))
  s.chr = readFasta(file.chr)
  s.chr = seq2nt(s.chr)
  
  # Small sequences  
  idx.small = which((sv.pos.all$single == 1) & 
                      (sv.pos.all$len >= min.len) & 
                      (sv.pos.all$len < big.len) & (sv.pos.all$chr == i.chr))
  # print(head(idx.small))
  for(irow in idx.small){
    s.tmp = s.chr[(sv.pos.all$beg[irow] + 1):(sv.pos.all$end[irow] - 1) ]
    # print('---')
    # print(s.tmp)
    # print((sv.pos.all$beg[irow] + 1))
    # print(sv.pos.all$end[irow] - 1)
    # print(sum(s.tmp == 'N'))
    # print((0.5 * length(s.tmp)))
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
  pokaz('Number of big sequences', length(seqs.big))
  pokaz('Number of small sequences', length(seqs.small))
}

writeFasta(seqs.small, file.sv.small)
writeFasta(seqs.big, file.sv.big)


# ---- GFF densities ----

# ***********************************************************************
# ---- Manual testing ----

if(F){
  library(rhdf5)
source(system.file("pannagram/utils/utils.R", package = "pannagram"))
  path.cons = './'
  ref.pref = '0'  
} 


