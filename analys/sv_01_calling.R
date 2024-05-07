# Find SVs and create GFF file

suppressMessages({ library(Biostrings)
  library(rhdf5)
  library('seqinr')
  library('foreach')
  library(doParallel)
  library("optparse")
})

source("utils/utils.R")

pokazStage('Get SV positions, GFF files, dencity files and consensys sequences')
pokazAttention('Be sure, that consensus sequence for the pangenome chromosomes have been generated')

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--ref.pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to directory with the consensus", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("--aln.type"), type="character", default="default", 
              help="type of alignment ('msa_', 'comb_', 'v_', etc)", metavar="character"),
  make_option(c("--acc.anal"), type = "character", default = NULL,
              help = "files with accessions to analyze", metavar = "character"),
  make_option(c("--stat.only"), type = "character", default = NULL,
              help = "files with accessions to analyze", metavar = "character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);


# If only the statistics is needed
if (!is.null(opt$stat.only)) {
  flag.stat = F
} else {
  flag.stat = T
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
  aln.type = 'msa_'
}

# Reference genome
if (is.null(opt$ref.pref) || (opt$ref.pref == 'NULL')) {
  ref.pref = NULL
  # stop("ref.pref is NULL")
} else {
  ref.pref <- opt$ref.pref
}

if (!is.null(opt$path.cons)) path.cons <- opt$path.cons

path.sv = paste(path.cons, 'sv/', sep = '')
if (!dir.exists(path.sv)) dir.create(path.sv)
path.gff = paste(path.sv, 'gff/', sep = '')
if (!dir.exists(path.gff)) dir.create(path.gff)

gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'

cutoff = 0.90


# ---- Combinations of chromosomes query-base to create the alignments ----

if(is.null(ref.pref)){
  
  pokaz('Reference genome:', ref.pref)
  
  s.pattern <- paste("^",aln.type,"\\d+_\\d+[^.]*\\.h5$", sep = '')
  files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
  
  # Extract reference names
  files.suff <- sapply(files, function(filename) {
    matches <- regmatches(filename, regexec(paste(aln.type, "\\d+_\\d+([^.]*)\\.h5", sep =""), filename))
    return(matches[[1]][2])  
  })
  if(length(unique(files.suff)) != 1){
    stop('Specify the genome, which was used for sorting')
  }

  pref.combinations <- sapply(files, function(filename) {
    matches <- regmatches(filename, regexec(paste(aln.type, "(\\d+)_(\\d+)[^.]*\\.h5", sep = ''), filename))
    return(paste(matches[[1]][2], matches[[1]][3], sep = "_"))
  })
  names(pref.combinations) = NULL
  
} else {
  pokaz('Genome for sorting:', ref.pref)
  # Old working version
  s.pattern <- paste("^", aln.type, ".*", '_ref_', ref.pref, sep = '')
  files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
  pref.combinations = gsub(aln.type, "", files)
  pref.combinations <- sub("_ref.*$", "", pref.combinations)
  pref.combinations <- pref.combinations[grep("^[0-9]+_[0-9]+$", pref.combinations)]
  
}


if(length(pref.combinations == 0)){
  pokazAttention('No Combinations found.')
} else {
  pokaz('Combinations', pref.combinations)  
}


# ---- Positions of SVs ----

sv.pos.all = c()
sv.beg.all = c()
sv.end.all = c()

flag.for = T
for(s.comb in pref.combinations){
  
  pokaz('* Combination', s.comb)
  
  # Get file for the combination
  if(!is.null(ref.pref)){
    file.comb = paste(path.cons, aln.type, s.comb,'_ref_',ref.pref,'.h5', sep = '')
  } else {
    s.pattern.comb <- paste("^",aln.type,s.comb,"[^.]*\\.h5$", sep = '')
    file.comb <- list.files(path = path.cons, pattern = s.pattern.comb, full.names = FALSE)
    file.comb = paste(path.cons, file.comb, sep = '')
    pokaz(file.comb)
  }
  
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
    pokaz('Sequence of accession', acc)
    v = h5read(file.comb, paste(gr.accs.e, acc, sep = ''))
    sv.cover = sv.cover + ((v == 0) & (!is.na(v))) * 1
  }
  
  # SV groups
  pokaz('Create SV groups...')
  pos.non.sv = (sv.cover == 0)
  idx.beg = which((pos.non.sv[-length(pos.non.sv)] == 1) & (pos.non.sv[-1] == 0))+ 1
  idx.end = which((pos.non.sv[-length(pos.non.sv)] == 0) & (pos.non.sv[-1] == 1)) 
  if(pos.non.sv[1] == 0) idx.end = idx.end[-1]
  idx.beg = idx.beg[1:length(idx.end)]
  len.sv = idx.end - idx.beg + 1
  
  sv.pos = data.frame(beg = idx.beg-1, end = idx.end+1)
  sv.pos$len = abs(sv.pos$beg - sv.pos$end) - 1 # do not change
  
  sv.beg = c()
  sv.end = c()
  for(acc in accessions){
    pokaz('Sequence of accession', acc)
    v = h5read(file.comb, paste(gr.accs.e, acc, sep = ''))
    sv.beg = cbind(sv.beg, v[sv.pos$beg])
    sv.end = cbind(sv.end, v[sv.pos$end])
  }
  colnames(sv.beg) = accessions
  colnames(sv.end) = accessions
  
  sv.len.acc = abs(sv.beg - sv.end) - 1  # do not change
  colnames(sv.len.acc) = accessions
  sv.pos$freq.min = 0
  sv.pos$freq.max = 0
  for(i.col in 1:n.acc){
    sv.pos$freq.min = sv.pos$freq.min + 1*((sv.len.acc[,i.col] <= (1-cutoff) * sv.pos$len) & (!is.na(sv.len.acc[,i.col])))
    sv.pos$freq.max = sv.pos$freq.max + 1*((sv.len.acc[,i.col] >= cutoff * sv.pos$len) & (!is.na(sv.len.acc[,i.col])))
  }
  
  sv.pos$freq.na = rowSums(is.na(sv.len.acc))
  
  sv.pos$freq.sum = sv.pos$freq.min + sv.pos$freq.max + sv.pos$freq.na
  sv.pos$single = (sv.pos$freq.sum == n.acc) * 1
  sv.pos = cbind(sv.pos, sv.len.acc[,1:n.acc])
  
  # Clean up NA
  sv.na = (sv.pos$freq.na == n.acc) | (sv.pos$freq.min == 0) | (sv.pos$freq.max == 0)
  sv.pos = sv.pos[!sv.na,]
  sv.beg = sv.beg[!sv.na,]
  sv.end = sv.end[!sv.na,]
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


file.sv.pos = paste(path.sv, 'sv_pangen_pos.rds', sep='')
saveRDS(sv.pos.all, file.sv.pos)

if(flag.stat){
  pokaz('Stat was generated')
  quit(save="no")
} 
## ---- Stop for Stat ----

file.sv.pos.beg = paste(path.sv, 'sv_pangen_beg.rds', sep='')
file.sv.pos.end = paste(path.sv, 'sv_pangen_end.rds', sep='')
saveRDS(sv.beg.all, file.sv.pos.beg)
saveRDS(sv.end.all, file.sv.pos.end)


# ---- GFF files ----

sv.version = 6
file.sv.gff = paste(path.gff, 'svs_pangen_v',sprintf("%02d", sv.version),'.gff', sep = '')

sv.pos.all$V10 = 1:nrow(sv.pos.all)
rownames(sv.pos.all) = sv.pos.all$gr

## ---- Single-event ----
sv.se = sv.pos.all[sv.pos.all$single == 1,]
sv.se.type = rep('indel', nrow(sv.se))

n.acc = length(accessions)
val.insert = c(1:n.acc)[(1:n.acc) <= (0.11 * n.acc)]
val.delet = c(1:n.acc)[(1:n.acc) >= (1 - 0.11 * n.acc)]

sv.se.type[sv.se$freq.max %in% val.insert] = 'insertion'
sv.se.type[sv.se$freq.max %in% val.delet] = 'deletion'

sv.annot = paste('ID=', sv.se$gr,  
                 # ';te=', sv.se$te, 
                 ';presence=',sv.se$freq.max,sep = '',
                 ';len_init=', sv.se$len)

sv.se.gff = data.frame(V1 = paste('PanGen_Chr', sv.se$chr, sep = ''),
                       V2 = 'pannagram',
                       V3 = sv.se.type, 
                       V4 = sv.se$beg + 1, 
                       V5 = sv.se$end - 1,
                       V6 = '.', V7 = '+', V8 = '.', 
                       V9 = sv.annot, 
                       V10 = sv.pos.all[sv.se$gr, 'V10'])

## ---- Multiple-event ----
sv.me = sv.pos.all[sv.pos.all$single != 1,]

s.multi = 'multi'
sv.me.gff = data.frame(V1 = paste('PanGen_Chr', sv.me$chr, sep = ''),
                       V2 = 'pannagram',
                       V3 = s.multi, V4 = sv.me$beg, V5 = sv.me$end,
                       V6 = '.', V7 = '+', V8 = '.', 
                       V9 = paste('ID=', sv.me$gr, ';len_init=', sv.me$len, sep = ''), 
                       V10 = sv.pos.all[sv.me$gr, 'V10'])

sv.gff = rbind(sv.se.gff, sv.me.gff)
sv.gff = sv.gff[order(sv.gff$V10),]

options(scipen = 999)
write.table(sv.gff[,1:9], file.sv.gff, quote = F, row.names = F, col.names = F, sep = '\t')
options(scipen = 0)

# ---- GFF In accessions ----


for(i.acc in 1:length(accessions)){
  acc = accessions[i.acc]
  pokaz('Generate GFF for accession', acc)
  
  df = sv.gff
  df$V1 = paste(acc, '_Chr', sv.pos.all$chr, sep = '')
  df$V4 = sv.beg.all[,acc] + 1
  df$V5 = sv.end.all[,acc] - 1
  df$V9 = paste('ID=', sv.me$gr, '.', acc, 
                ';len_init=', sv.pos.all$len,
                ';len_acc=', abs(sv.end.all[,acc]-sv.beg.all[,acc])-1, sep = '')
  
  df = df[sv.pos.all$len > 0,]
  
  df = df[!is.na(df$V4),]
  df = df[!is.na(df$V5),]
  
  df = df[df$V4 != 0,]
  df = df[df$V5 != 0,]
  
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
path.seq = paste(path.cons, 'seq/', sep = '')
if (!dir.exists(path.seq)) {
  pokazAttention('Consensus sequence doesn’t exist')
  stop('Please generate the pangenome consensus sequence first')
}

min.len = 15
big.len = 50
max.len = 30000


file.sv.small =  paste(path.sv, 'seq_sv_small.fasta', sep = '')
file.sv.big =  paste(path.sv, 'seq_sv_big.fasta', sep = '')

pokaz(file.sv.small)

seqs.small = c()
seqs.big = c()
for(s.comb in pref.combinations){
  i.chr = comb2ref(s.comb)
  pokaz('Chromosome', i.chr)
  file.chr = paste(path.seq, 'seq_cons_', i.chr, '.fasta', sep = '')
  s.chr = readFastaMy(file.chr)
  s.chr = seq2nt(s.chr)
  
  # Small sequences  
  idx.small = which((sv.pos.all$single == 1) & 
                      (sv.pos.all$len >= min.len) & 
                      (sv.pos.all$len < big.len) & (sv.pos.all$chr == i.chr))
  print(head(idx.small))
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
  print(head(idx.big))
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

writeFastaMy(seqs.small, file.sv.small)
writeFastaMy(seqs.big, file.sv.big)


# ---- GFF densities ----

# ***********************************************************************
# ---- Manual testing ----

if(F){
  library(rhdf5)
  source('../../../pannagram/utils/utils.R')
  path.cons = './'
  ref.pref = '0'  
} 


