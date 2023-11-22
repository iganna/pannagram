# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
})

source("utils.R")
source("synteny_funcs.R")

pokazStage('Combine: alignments by chromosomes')

args = commandArgs(trailingOnly=TRUE)


option_list = list(
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to consensus directory", metavar="character"),
  make_option(c("-p", "--ref.pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("-n", "--n.chr.ref"), type="character", default=NULL, 
              help="number of chromosomes in the reference genome", metavar="character"),
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# Set the number of cores for parallel processing
num.cores.max = 10
num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))
myCluster <- makeCluster(num.cores, type = "PSOCK")
registerDoParallel(myCluster)


if (is.null(opt$ref.pref)) {
  stop("ref.pref is NULL")
} else {
  ref.pref <- opt$ref.pref
}

# ---- Combinations of chromosomes query-base to create the alignments ----

# path.cons = './'
# ref.pref = '0'

s.pattern <- paste("^", 'res_', ".*", '_ref_', ref.pref, sep = '')
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub("res_", "", files)
pref.combinations <- sub("_ref.*$", "", pref.combinations)

pokaz('Reference:', base.acc.ref)
pokaz('Combinations', pref.combinations)

# ----  Combine correspondence  ----

gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'
max.len.gap = 20000

#flag.for = F
#tmp = foreach(s.comb = pref.combinations, .packages=c('rhdf5', 'crayon'))  %dopar% {  # which accession to use
flag.for = T
for(s.comb in pref.combinations){
  
  query.chr = chromosome.pairs[i.chr.pair, 1]
  base.chr = chromosome.pairs[i.chr.pair, 2]
  
  base.len = chr.len[base.chr]
  
  file.comb = paste(path.cons, 'res_', s.comb,'_ref_',ref.pref,'.h5', sep = '')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b], 
                         groups1$name[groups1$group == gr.accs.b])  # full name of accessions
  # accessions <- sub("^.*_", "", accessions)
  
  
  
  if (file.exists(file.comb)) file.remove(file.comb)
  h5createFile(file.comb)
  
  # Path to accessions chunks
  gr.accs <- "accs/"
  # TODO: Check the availability of the group before creating it
  h5createGroup(file.comb, gr.accs)
  
  
  gr.break = 'break/'
  h5createGroup(file.comb, gr.break)
  
  idx.break = 0
  idx.gaps = rep(0, base.len)
  
  for(acc in accessions){
    
    pokaz('Accession', acc, 'qchr', query.chr, 'bchr', base.chr)
    
    pref.comb = paste0(acc, '_', query.chr, '_', base.chr, collapse = '')
    file.aln.full <- paste(path.aln, paste0(pref.comb,  '_full.rds', collapse = ''), sep = '')
    if(!file.exists(file.aln.full)) next
    
    # Reading the alignment
    x = readRDS(file.aln.full)
    
    # Get query coordinates in base order
    x.corr = getCorresp2BaseSign(x, base.len)
    
    # Write into file
    suppressMessages({
      h5write(x.corr, file.comb, paste(gr.accs, 'acc_', acc, sep = ''))
    })
    
    # ----  Find gaps  ----
    
    idx.gaps[x.corr == 0] = idx.gaps[x.corr == 0] + 1
    
    # ----  Find breaks  ----
    v = x.corr
    
    # Find blocks of additional breaks
    v = cbind(v, 1:length(v))                       # 2 - in ref-based coordinates
    v = v[v[,1] != 0,]                                   # 1 - existing coordinates of accessions
    v = cbind(v, 1:nrow(v))                       # 3 - ranked order in ref-based coordinates
    v = cbind(v, rank(abs(v[,1])) * sign(v[,1]))  # 4 - signed-ranked-order in accessions coordinates 
    
    # v = v[order(v[,1]),]  # not necessary
    
    # with the absence, but neighbouring
    idx.tmp = which( (abs(diff(v[,4])) == 1) &  # Neighbouring in accession-based order
                       (abs(diff(abs(v[,3])) == 1)) &  # Neighbouring in ref-based order
                       (abs(diff(v[,1])) <= max.len.gap) &  # Filtering by length in accession coordinates
                       (abs(diff(v[,2])) <= max.len.gap) &  # Filtering by length in reference coordinates
                       (abs(diff(v[,1])) > 1))  # NOT neighbouring in accession-specific coordinates
    
    # Fix (beg < end) order
    idx.tmp.acc = data.frame(beg = v[idx.tmp,2], end = v[idx.tmp+1,2], acc = acc)
    idx.ord = which(idx.tmp.acc$beg > idx.tmp.acc$end)
    if(length(idx.ord) > 0){
      tmp = idx.tmp.acc$beg[idx.ord]
      idx.tmp.acc$beg[idx.ord] = idx.tmp.acc$end[idx.ord]
      idx.tmp.acc$end[idx.ord] = tmp
    }
    # idx.tmp.acc = idx.tmp.acc[order(idx.tmp.acc$beg),]  # order ONLY if ordered before
    
    # Remove overlaps
    idx.overlap = which( (idx.tmp.acc$beg[-1] - idx.tmp.acc$end[-nrow(idx.tmp.acc)]) <= 3)
    
    i.cnt = 0
    if(length(idx.overlap) > 0){
      j.ov = 0
      for(i.ov in idx.overlap){
        if(i.ov <= j.ov) next
        j.ov = i.ov + 1
        while(j.ov %in% idx.overlap){
          j.ov = j.ov + 1
        }
        # print(c(i.ov, j.ov))
        i.cnt = i.cnt + 1
        idx.tmp.acc$end[i.ov] = idx.tmp.acc$end[j.ov]
      }
      idx.tmp.acc = idx.tmp.acc[-(idx.overlap+1),]
    }
    
    
    # Write into file
    suppressMessages({
      h5write(idx.tmp.acc, file.comb, paste(gr.break, 'acc_', acc, sep = ''))
    })
    
    
    # Fill up positions with breaks
    idx.break.acc = rep(0, base.len)
    idx.break.acc[idx.tmp.acc$beg] = 1
    idx.break.acc[idx.tmp.acc$end] = -1
    idx.break.acc = cumsum(idx.break.acc)
    idx.break.acc[idx.tmp.acc$end] = 1
    
    # Save breaks
    idx.break = idx.break + idx.break.acc
    
    rmSafe(x.corr)
    rmSafe(x)
    rmSafe(v)
    rmSafe(idx.tmp.acc)
    rmSafe(idx.break.acc)
    
  }
  
  suppressMessages({
    h5write(idx.break, file.comb, 'breaks_all')
    h5write(idx.gaps, file.comb, 'gaps_all')
    h5write(base.acc.ref, file.comb, 'ref')
    
    h5write(1:base.len, file.comb, paste(gr.accs, 'acc_', base.acc.ref, sep = ''))
    # h5write(NULL, file.comb, paste(gr.break, base.acc.ref, sep = ''))
  })
  
  rmSafe(idx.break)
  rmSafe(idx.gaps)
  
  H5close()
  gc()
  
}





