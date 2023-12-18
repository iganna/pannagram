#' How the output files look like:
#'     group           name       otype   dclass      dim
#'         /           accs   H5I_GROUP                  
#'     /accs              0 H5I_DATASET    FLOAT 28940631
#'     /accs          10002 H5I_DATASET    FLOAT 28940631
#'     /accs          10015 H5I_DATASET    FLOAT 28940631

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
})

source("utils/utils.R")
source("pangen/synteny_funcs.R")

pokazStage('Step 7. Combine reference-based alignments by chromosomes')

args = commandArgs(trailingOnly=TRUE)


option_list = list(
  make_option(c("--path.chr.len"), type="character", default=NULL, 
              help="file with lengths of chromosomes", metavar="character"),
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to consensus directory", metavar="character"),
  make_option(c("-o", "--path.aln"), type="character", default=NULL, 
              help="path to the output directory with alignments", metavar="character"),
  make_option(c("-p", "--pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("-n", "--n.chr.ref"), type="character", default=NULL, 
              help="number of chromosomes in the reference genome", metavar="character"),
  make_option(c("-m", "--n.chr.acc"), type="character", default=NULL, 
              help="number of chromosomes in the accessions", metavar="character"),
  make_option(c("-a", "--all.vs.all"), type="character", default=NULL, 
              help="alignment of all chromosomes vs all or not: T/F", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("-r", "--path.ref"), type="character", default=NULL, 
              help="path to the reference file", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="type of fasta files", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# Set the number of cores for parallel processing
num.cores.max = 10
num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))
myCluster <- makeCluster(num.cores, type = "PSOCK")
registerDoParallel(myCluster)

# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) system(paste('mkdir ', path.cons, sep = ''))

# File with chromosomal lengths, which should be in the consensus dir
path.chr.len = ifelse(!is.null(opt$path.chr.len), opt$path.chr.len, 'chr_len/')
path.chr.len = paste(path.cons, path.chr.len, sep = '')
if(!dir.exists(path.chr.len)) system(paste('mkdir ', path.chr.len, sep = ''))

if (!is.null(opt$path.ref)) path.base <- opt$path.ref  # to know the chromosomal lengths
if (!is.null(opt$pref)) base.acc.ref <- opt$pref

# Path with alignments
if (!is.null(opt$path.aln)) path.aln <- opt$path.aln

# To create combinations
if (!is.null(opt$n.chr.ref)) n.chr.ref <- as.numeric(opt$n.chr.ref)
if (!is.null(opt$n.chr.acc)) n.chr.acc <- as.numeric(opt$n.chr.acc)
if (!is.null(opt$all.vs.all)) all.vs.all <- as.logical(opt$all.vs.all)


# ---- Get accession names ----

aln.suff <- "_full.rds"
aln.files <- list.files(path.aln)

pokaz('Files', aln.files)

aln.files <- aln.files[grep(paste0(aln.suff, "$"), aln.files)]

accessions <- sapply(aln.files, function(filename){
  parts <- unlist(strsplit(filename, "_", fixed = TRUE))
  name <- paste(parts[1:(length(parts) - 3)], collapse = "_")
  return(name)})
names(accessions) = NULL



accessions <- sort(unique(accessions))

pokaz('Accessions:', accessions)

# ---- Combinations of chromosomes query-base to create the alignments ----

chromosome.pairs <- unique(do.call(rbind, lapply(aln.files, function(filename){
  parts <- unlist(strsplit(filename, "_", fixed = TRUE))
  s.comb <- c(as.numeric(parts[length(parts) - 2]),
              as.numeric(parts[length(parts) - 1]))
  return(s.comb)})))

pokaz('Combinations:', paste(chromosome.pairs[,1], chromosome.pairs[,2], sep = '_'))

# ---- Length of reference chromosomes ----

file.chr.len = paste(path.chr.len, base.acc.ref, '_len.rds', sep = '')
# pokaz('File with chromosomal lengths', file.chr.len)
if(!file.exists(file.chr.len)){
  chr.len = c()
  for(i.chr in 1:n.chr.ref){
    acc = base.acc.ref
    # print(c(i.chr, acc))
    
    # Read base chromosome
    base.file = paste0(base.acc.ref, '_chr', i.chr , '.fasta', collapse = '')
    # pokaz('Base:', base.file)
    base.fas.fw = readFastaMy(paste(path.base, base.file, sep = ''))
    base.fas.fw = seq2nt(base.fas.fw)
    chr.len = c(chr.len, length(base.fas.fw))
  }
  saveRDS(chr.len, file.chr.len)
} else {
  chr.len = readRDS(file.chr.len)
}

# pokaz('Chromosomal lengths', chr.len)

# ----  Combine correspondence  ----

pokaz('Reference:', base.acc.ref)


max.len.gap = 20000

flag.for = F
tmp = foreach(i.chr.pair = 1:nrow(chromosome.pairs), .packages=c('rhdf5', 'crayon'))  %dopar% {  # which accession to use
 # flag.for = T
 # for(i.chr.pair in 1:nrow(chromosome.pairs)){
  
  query.chr = chromosome.pairs[i.chr.pair, 1]
  base.chr = chromosome.pairs[i.chr.pair, 2]
  
  base.len = chr.len[base.chr]
  
  file.comb = paste(path.cons, 'comb_', query.chr, '_', base.chr,'_ref_',base.acc.ref,'.h5', sep = '')
  if (file.exists(file.comb)) file.remove(file.comb)
  h5createFile(file.comb)
  
  # Path to accessions chunks
  gr.accs <- "accs/"
  # TODO: Check the availability of the group before creating it
  h5createGroup(file.comb, gr.accs)
  
  
  # gr.break = 'break/'
  # h5createGroup(file.comb, gr.break)
  
  idx.break = 0
  # idx.gaps = rep(0, base.len)
  
  for(acc in accessions){
    
    # pokaz('Accession', acc, 'qchr', query.chr, 'bchr', base.chr)
    
    pref.comb = paste0(acc, '_', query.chr, '_', base.chr, collapse = '')
    file.aln.full <- paste(path.aln, paste0(pref.comb,  '_full.rds', collapse = ''), sep = '')
    if(!file.exists(file.aln.full)) next
    
    # Reading the alignment
    x = readRDS(file.aln.full)
    
    # Get query coordinates in base order
    x.corr = getCorresp2BaseSign(x, base.len)
    
    if(sum(duplicated(x.corr[x.corr != 0])) > 0) stop('DUPLICSTIONS')
    
    # Write into file
    suppressMessages({
      h5write(x.corr, file.comb, paste(gr.accs, '', acc, sep = ''))
    })
    
    # # ----  Find gaps  ----
    # 
    # idx.gaps[x.corr == 0] = idx.gaps[x.corr == 0] + 1
    # 
    # # ----  Find breaks  ----
    # v = x.corr
    # 
    # # Find blocks of additional breaks
    # v = cbind(v, 1:length(v))                       # 2 - in ref-based coordinates
    # v = v[v[,1] != 0,]                                   # 1 - existing coordinates of accessions
    # v = cbind(v, 1:nrow(v))                       # 3 - ranked order in ref-based coordinates
    # v = cbind(v, rank(abs(v[,1])) * sign(v[,1]))  # 4 - signed-ranked-order in accessions coordinates 
    # 
    # # v = v[order(v[,1]),]  # not necessary
    # 
    # # with the absence, but neighbouring
    # idx.tmp = which( (abs(diff(v[,4])) == 1) &  # Neighbouring in accession-based order
    #                    (abs(diff(abs(v[,3])) == 1)) &  # Neighbouring in ref-based order
    #                    (abs(diff(v[,1])) <= max.len.gap) &  # Filtering by length in accession coordinates
    #                    (abs(diff(v[,2])) <= max.len.gap) &  # Filtering by length in reference coordinates
    #                    (abs(diff(v[,1])) > 1))  # NOT neighbouring in accession-specific coordinates
    # 
    # # Fix (beg < end) order
    # idx.tmp.acc = data.frame(beg = v[idx.tmp,2], end = v[idx.tmp+1,2], acc = acc)
    # idx.ord = which(idx.tmp.acc$beg > idx.tmp.acc$end)
    # if(length(idx.ord) > 0){
    #   tmp = idx.tmp.acc$beg[idx.ord]
    #   idx.tmp.acc$beg[idx.ord] = idx.tmp.acc$end[idx.ord]
    #   idx.tmp.acc$end[idx.ord] = tmp
    # }
    # # idx.tmp.acc = idx.tmp.acc[order(idx.tmp.acc$beg),]  # order ONLY if ordered before
    # 
    # # Remove overlaps
    # idx.overlap = which( (idx.tmp.acc$beg[-1] - idx.tmp.acc$end[-nrow(idx.tmp.acc)]) <= 3)
    # 
    # i.cnt = 0
    # if(length(idx.overlap) > 0){
    #   j.ov = 0
    #   for(i.ov in idx.overlap){
    #     if(i.ov <= j.ov) next
    #     j.ov = i.ov + 1
    #     while(j.ov %in% idx.overlap){
    #       j.ov = j.ov + 1
    #     }
    #     # print(c(i.ov, j.ov))
    #     i.cnt = i.cnt + 1
    #     idx.tmp.acc$end[i.ov] = idx.tmp.acc$end[j.ov]
    #   }
    #   idx.tmp.acc = idx.tmp.acc[-(idx.overlap+1),]
    # }
    # 
    # 
    # # Write into file
    # suppressMessages({
    #   h5write(idx.tmp.acc, file.comb, paste(gr.break, '', acc, sep = ''))
    # })
    # 
    # 
    # # Fill up positions with breaks
    # idx.break.acc = rep(0, base.len)
    # idx.break.acc[idx.tmp.acc$beg] = 1
    # idx.break.acc[idx.tmp.acc$end] = -1
    # idx.break.acc = cumsum(idx.break.acc)
    # idx.break.acc[idx.tmp.acc$end] = 1
    # 
    # # Save breaks
    # idx.break = idx.break + idx.break.acc
    
    rmSafe(x.corr)
    rmSafe(x)
    rmSafe(v)
    rmSafe(idx.tmp.acc)
    # rmSafe(idx.break.acc)
    
  }
  
  suppressMessages({
    # h5write(idx.break, file.comb, 'breaks_all')
    # h5write(idx.gaps, file.comb, 'gaps_all')
    h5write(base.acc.ref, file.comb, 'ref')
    
    h5write(1:base.len, file.comb, paste(gr.accs, '', base.acc.ref, sep = ''))
    # h5write(NULL, file.comb, paste(gr.break, base.acc.ref, sep = ''))
  })
  
  # rmSafe(idx.break)
  rmSafe(idx.gaps)
  
  H5close()
  gc()
  
}





