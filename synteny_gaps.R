suppressMessages({
  library(Biostrings)
  library(seqinr)
  library(spgs)  # reverseComplement("actg")
  library(foreach)
  library(doParallel)
  library(stringr)
  library(optparse)
})


source("utils.R")
source("synteny_funcs.R")

pokazStage('Alignment2: fill the gaps between synteny blocks')


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


#' ============================================================================

files.query = list.files(path = path.query, pattern = "\\.fasta$")
query.name = gsub("*.fasta","",files.query)

query.name = unique(sapply(query.name, function(s){ strsplit(s, '_')[[1]][1] }))
names(query.name) <- NULL

#query.name = c('0')

pokaz('Accessions:', query.name)
pokaz('Base accession:', base.acc)


# Combinations of chromosomes query-base to chreate the alignments
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



for.flag = F
tmp = foreach(i.chr.pair = 1:nrow(chromosome.pairs), .packages=c('stringr','Biostrings', 'seqinr', 'spgs'))  %dopar% {  # which accession to use

#for.flag = T
#for(i.chr.pair in 1:nrow(chromosome.pairs)){
  
  pokaz('A pair number', i.chr.pair, ':')
  
  i.query = chromosome.pairs[i.chr.pair, 1]
  
  if(query.name[i.query] == base.acc){
    if(for.flag) next
    return(NULL)
  }
  
  query.chr = chromosome.pairs[i.chr.pair, 2]
  base.chr = chromosome.pairs[i.chr.pair, 3]
  
  pref.comb = paste0(query.name[i.query], '_', query.chr, '_', base.chr, collapse = '')
  
  file.aln.pre <- paste(path.aln, paste0(pref.comb, '_maj.rds', collapse = ''), sep = '')
  if(!file.exists(file.aln.pre)) {
    if(for.flag) next
    return(NULL)
  }
  
  base.file = paste0(base.acc, '_chr', base.chr , '.', base.suff, collapse = '')
  pokaz('Base:', base.file, sep=' ')
  
  file.aln.full <- paste(path.aln, paste0(pref.comb,  '_full.rds', collapse = ''), sep = '')
  if(file.exists(file.aln.full)) {
    if(for.flag) next
    return(NULL)
  }
  
  pokaz('Alignment:', query.name[i.query], query.chr, base.chr)
  
  # # Read reference sequences
  # base.file = paste0(base.acc, '_chr', base.chr , '.', base.suff, collapse = '')
  # pokaz('Base:', base.file)
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
  
  
  x = readRDS(file = file.aln.pre) 
  
  file.gaps.out = paste0(path.gaps,
                         'acc_', query.name[i.query], 
                         '_qchr', query.chr, '_bchr', base.chr, '_out.txt', collapse = '')
  
  if(!file.exists(file.gaps.out)){
    saveRDS(object = x, file = file.aln.full) 
    if(for.flag) next
    return(NULL)
  }
  
  all_lines_start_with_hash <- system(sprintf("grep -v '^#' %s | wc -l", file.gaps.out), intern = TRUE) == 0
  if(all_lines_start_with_hash){
    if(for.flag) next
    return(NULL)
  }
  
  t.gap <- read.table(file.gaps.out, stringsAsFactors = F, header = F)
  
  
  start.pos.query = as.numeric(sapply(strsplit(t.gap$V1, "\\|"), "[", 5)) - 1
  start.pos.base = as.numeric(sapply(strsplit(t.gap$V10, "\\|"), "[", 5)) - 1
  t.gap[,2:3] = t.gap[,2:3] + start.pos.query
  t.gap[,4:5] = t.gap[,4:5] + start.pos.base
  # t.gap$base = t.gap$V10
  t.gap$V10 = x$V10[1]
  t.gap <- setDir(t.gap, base.len)
  t.gap <- t.gap[order(t.gap$V2),]
  
  checkCorrespToGenome(t.gap[1:min(1000, nrow(t.gap)),], query.fas = query.fas.chr,
                       base.fas.fw = base.fas.fw,
                       base.fas.bw = base.fas.bw)
  
  
  # wnd.q = sapply(strsplit(t.gap$V1, "\\|"), "[", 1)
  # wnd.q = as.numeric(sapply(strsplit(wnd.q, "_"), "[", 2))
  # idx.q = sort(unique(wnd.q))
  
  x = x[order(x$V2),]
  irow = 1
  idx.gap.used = c()
  while(irow < nrow(x)){
    p1 = x$V3[irow]
    p2 = x$V2[irow+1]
    idx.gap = which((t.gap$V3 >= p1) & (t.gap$V3 <= p2))
    idx.gap = setdiff(idx.gap, idx.gap.used)
    idx.gap.used = c(idx.gap.used, idx.gap)
    if(length(idx.gap) == 0){
      irow = irow + 1
      next
    }
  
    t.gap.tmp = t.gap[idx.gap, ]
  
    idx = c(irow, irow + 1)
    t.gap.tmp =  rbind(x[idx,], t.gap.tmp)
    t.gap.tmp = t.gap.tmp[order(t.gap.tmp$V2),]
    
    t.gap.tmp = removeOverlapdAfterBlastGap(t.gap.tmp)
    t.gap.tmp = t.gap.tmp[order(t.gap.tmp$V2),]
    t.gap.tmp = glueZero(t.gap.tmp)
    t.gap.tmp = t.gap.tmp[order(t.gap.tmp$V2),]
    t.gap.tmp = cleanJumps2(t.gap.tmp, min(t.gap.tmp$V7[c(1, nrow(t.gap.tmp))]))
    
    
    t.gap.tmp = glueZero(t.gap.tmp)
    t.gap.tmp <- glueAlmostZero(t.gap.tmp, threshold = 2000,
                                query.fas.chr = query.fas.chr,
                                base.fas.fw = base.fas.fw,
                                base.fas.bw = base.fas.bw)
    if(nrow(t.gap.tmp) != 1){
      thresholds = c(100, 500, 1000, 1500, 2000)
      # print('glueByThreshold...')
      t.gap.tmp = glueByThreshold(t.gap.tmp, thresholds, query.fas = query.fas.chr,
                                  base.fas.fw = base.fas.fw,
                                  base.fas.bw = base.fas.bw, file.log=file.log)
      
      t.gap.tmp = removeCompleteOverlaps2(t.gap.tmp, base.len)
      t.gap.tmp = removeShortOverlaps2(t.gap.tmp, base.len)
    }
    
    x = x[-idx,]
    x = rbind(x, t.gap.tmp)
    x = x[order(x$V2),]
  }
  
  
  # t <- hardCleaning(t, base.len)
  
  checkCorrespToGenome(x, query.fas = query.fas.chr,
                       base.fas.fw = base.fas.fw,
                       base.fas.bw = base.fas.bw)
  
  saveRDS(object = x, file = file.aln.full) 

}  # accessions




