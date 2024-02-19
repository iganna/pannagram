suppressMessages({ library(Biostrings)
library(rhdf5)
library('seqinr')
library('foreach')
library(doParallel)
library("optparse")
source('utils/utils.R')
})

pokazStage('Step 12. Combine all alignments together into the final one')

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.mafft.in"), type="character", default=NULL, 
              help="path to directory, where to combine fasta files for mafft runs", metavar="character"),
  make_option(c("-p", "--ref.pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("--path.mafft.out"), type="character", default=NULL, 
              help="path to directory, where to mafft results are", metavar="character"),
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to directory with the consensus", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores.max = 10
num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))

# Reference genome
if (is.null(opt$ref.pref)) {
  stop("ref.pref is NULL")
} else {
  ref.pref <- opt$ref.pref
}

if (!is.null(opt$path.mafft.in)) path.mafft.in <- opt$path.mafft.in
if (!is.null(opt$path.mafft.out)) path.mafft.out <- opt$path.mafft.out
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons


# ***********************************************************************
# ---- Preparation ----

n.flank = 30

gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'

max.block.elemnt = 3 * 10^ 6


# ---- Combinations of chromosomes query-base to create the alignments ----


s.pattern <- paste("^", 'res_', ".*", '_ref_', ref.pref, sep = '')
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pokaz('Path', path.cons)
pokaz('Files', files)
pref.combinations = gsub("res_", "", files)
pref.combinations <- sub("_ref.*$", "", pref.combinations)

pokaz('Reference:', ref.pref)
pokaz('Combinations', pref.combinations)


# ***********************************************************************
# ---- MAIN program body ----

stat.comb <- data.frame(comb = character(),
                        coverage = numeric(),
                        stringsAsFactors = FALSE)

for(s.comb in pref.combinations){
  
  pokaz('* Combination', s.comb)
  
  # Get accessions
  file.comb = paste(path.cons, 'res_', s.comb,'_ref_',ref.pref,'.h5', sep = '')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)
  base.len = length(h5read(file.comb, paste(gr.accs.e, accessions[1], sep = '')))
  
  
  # n.rows.block = round(max.block.elemnt / n.acc)
  # seq.blocks = seq(1, base.len - n.rows.block, n.rows.block)
  
  # ---- All MAFFT results for the combination ----
  pref = paste('Gap', s.comb, sep = '_')
  mafft.res = data.frame(file = list.files(path = path.mafft.out, 
                                           pattern = paste('^', pref, '.*_flank_', n.flank, '_aligned.fasta$', sep='')))
  
  
  mafft.res$comb = lapply(mafft.res$file, function(s) strsplit(s, '_')[[1]])
  
  y = lapply(mafft.res$file, function(s) strsplit(s, '_')[[1]])
  z = t(matrix(unlist(y), ncol = length(y)))
  mafft.res$comb = paste(z[,2], z[,3], sep = '_')
  mafft.res$beg = as.numeric(z[,5])
  mafft.res$end = as.numeric(z[,6])
  mafft.res$id = as.numeric(z[,3])
  
  
  # ---- Get long alignment positions ----
  pokaz('Read Long alignments..')
  idx.skip = c()
  mafft.aln.pos = list()
  for(i in 1:nrow(mafft.res)){
    # if((i %% 100) == 0) pokaz('Aln', i)
    file.aln = paste(path.mafft.out, mafft.res$file[i], sep = '')
    
    if(!file.exists(file.aln)) {
      idx.skip = c(idx.skip, i)
      next
    }

    # command <- paste("wc -l", file.aln)
    # result <- system(command, intern = TRUE)
    # num.lines <- as.integer(strsplit(result, " ")[[1]][1])
    # if(num.lines < 2) {
    #   idx.skip = c(idx.skip, i)
    #   next
    # }

    
    # pokaz(file.aln)
    aln.seq = readFastaMy(file.aln)
    n.aln.seq = length(aln.seq)
    name.aln.seq = names(aln.seq)
    name.acc = sapply(name.aln.seq, function(s) strsplit(s, '\\|')[[1]][1])
    pos.aln = sapply(name.aln.seq, function(s) strsplit(s, '\\|')[[1]][3:4])
    
    
    len.aln.seq <- nchar(aln.seq[1])
    # aln.mx <- matrix(, nrow = n.aln.seq, ncol = len.aln.seq)
    pos.mx <- matrix(0, nrow = n.aln.seq, ncol = len.aln.seq)
    for (i.seq in 1:length(aln.seq)) {
      tmp = strsplit(aln.seq[i.seq], "")[[1]]
      tmp.nongap = which(tmp != '-')
      # tmp[tail(tmp.nongap, (n.flank) )] = '-'
      # tmp[tmp.nongap[1:(n.flank) ]] = '-'
      # aln.mx[i, ] <- tmp

      tmp.nongap = tmp.nongap[-(1:(n.flank))]
      tmp.nongap <- tmp.nongap[1:(length(tmp.nongap) - n.flank)]

      p1 = pos.aln[1, i.seq]
      p2 = pos.aln[2, i.seq]
      pos.tmp = p1:p2
      pos.mx[i.seq, tmp.nongap] = pos.tmp
    }
    # aln.mx = aln.mx[,colSums(aln.mx != '-') != 0]
    pos.mx = pos.mx[,colSums(pos.mx != 0) != 0]
    row.names(pos.mx) = name.acc

    mafft.aln.pos[[i]] = pos.mx
  }
  
  if(length(idx.skip) > 0){
    mafft.aln.pos = mafft.aln.pos[-idx.skip]
    mafft.res = mafft.res[-idx.skip,]
  }
  
  # warnings()
  mafft.res$len = unlist(lapply(mafft.aln.pos, ncol))
  mafft.res$extra = mafft.res$len - (mafft.res$end - mafft.res$beg - 1)
  # if(min(mafft.res$extra) < 0) stop('Long: Wrong lengths of alignment and gaps')
  mafft.res$extra[mafft.res$extra < 0] = 0
  
  # ---- Short alignments ----
  pokaz('Read Short alignments..')
  msa.res = readRDS(paste(path.cons, 'aln_short_', s.comb, '.rds', sep = ''))
  msa.res$len = unlist(lapply(msa.res$aln, nrow))
  msa.res$extra = msa.res$len - (msa.res$ref.pos$end - msa.res$ref.pos$beg - 1)
  # if(min(msa.res$extra) < 0) stop('Short: Wrong lengths of alignment and gaps')
  msa.res$extra[msa.res$extra < 0] = 0

  # ---- Singletons alignments ----
  pokaz('Read Singletons..')
  single.res = readRDS(paste(path.cons, 'singletons_', s.comb, '.rds', sep = ''))
  single.res$len = rowSums(single.res$pos.end) - rowSums(single.res$pos.beg)  + 1
  single.res$extra = single.res$len - (single.res$ref.pos$end - single.res$ref.pos$beg - 1)
  # if(min(single.res$extra) < 0) stop('Wrong lengths of alignment and gaps')
  single.res$extra[single.res$extra < 0] = 0
  
  # ---- Analysis of positions ----
  # Here I wouls like fo find function of positions corresponcences between 4 things: 
  # old coordinates, long, short and singleton coordinates
  
  n.shift = rep(0, base.len)
  
  n.shift[mafft.res$end] = mafft.res$extra  # Long extra
  n.shift[msa.res$ref.pos$end] = msa.res$extra  # Short extra
  n.shift[single.res$ref.pos$end] = single.res$extra # Singletons extra
  n.shift = cumsum(n.shift)
  
  fp.main = (1:base.len) + n.shift
  
  # -- 
  # Singletons
  fp.single = list()
  for(i in 1:length(single.res$len)){
    n.pos = single.res$len[i] - 2
    fp.single[[i]] = fp.main[single.res$ref.pos$beg[i]] + (1:n.pos)
  }
  
  # Short
  fp.short = list()
  for(i in 1:length(msa.res$len)){
    n.pos = msa.res$len[i]
    fp.short[[i]] = fp.main[msa.res$ref.pos$beg[i]] + (1:n.pos)
  }
  
  # Check short
  for(i in 1:length(msa.res$len)){
    if(is.null(msa.res$aln[[i]])) next
    if(length(fp.short[[i]]) != nrow(msa.res$aln[[i]])) stop(i)
  }
  
  
  # Long
  fp.long = list()
  for(i in 1:length(mafft.aln.pos)){
    n.pos = ncol(mafft.aln.pos[[i]])
    fp.long[[i]] = fp.main[mafft.res$beg[i]] + (1:n.pos)
  }
  
  
  pos.beg.all = list(single.res$ref.pos$beg, msa.res$ref.pos$beg, mafft.res$beg)
  pos.end.all = list(single.res$ref.pos$end, msa.res$ref.pos$end, mafft.res$end)
  
  pos.delete.all = 0
  for(i.pos in 1:3){
    pos.beg = pos.beg.all[[i.pos]]
    pos.end = pos.end.all[[i.pos]]
    pos.delete = rep(0, base.len)
    pos.delete[pos.beg] = 1
    # pokaz(  sum(pos.delete == 1), sum(pos.delete == -1))
    pos.delete[pos.end] = pos.delete[pos.end]-1
    # pokaz(  sum(pos.delete == 1), sum(pos.delete == -1))
    pos.delete = cumsum(pos.delete)
    pos.delete[pos.beg] = 0
    pos.delete[pos.end] = 0
    
    pos.delete.all = pos.delete.all + pos.delete
  }
  pos.delete = pos.delete.all

  if(sum(unlist(pos.end.all) - unlist(pos.beg.all) - 1) != sum(pos.delete)) stop('Wrong identification of positions to delete')
  
  pos.remain = pos.delete == 0
  
  fp.main[pos.delete != 0] = 0
  
  base.len.aln = max(fp.main)
  
  # Check-points
  fp.add = c(unlist(fp.single),unlist(fp.short),unlist(fp.long))
  if(sum(duplicated(c(fp.main[fp.main != 0], fp.add))) != 0) stop('Something if wrotng with positions; Point A')
  # if(length(unique(c(fp.main, fp.add))) != (max(fp.main) + 1)) stop('Something if wrotng with positions; Point B')  # it's not trow anymore
  
  
  # ---- Define blocks before the big alignments ----
  # pos.beg = fp.main[mafft.res$beg] + 1
  # pos.beg.bins <- cut(pos.beg, breaks = c(seq.blocks, Inf), labels = FALSE)
  # pos.block.end = tapply(pos.beg, pos.beg.bins, max)
  # pos.block.end[length(pos.block.end)] = base.len
  
  file.res = paste(path.cons, 'msa_', s.comb,'_ref_',ref.pref,'.h5', sep = '')
  if (file.exists(file.res)) file.remove(file.res)
  h5createFile(file.res)
  
  suppressMessages({
    h5createGroup(file.res, gr.accs.e)
  })
  
  pos.nonzero = fp.main != 0
  for(acc in accessions){
    pokaz('Accession', acc)
    v = h5read(file.comb, paste(gr.accs.e, acc, sep = ''))
    v.aln = rep(0, base.len.aln)
    v.aln[fp.main[pos.nonzero]] = v[pos.nonzero]
    # v.aln[fp.main[pos.remain]] = v[pos.remain]
    
    # Add singletons
    for(i in 1:length(single.res$len)){
      if(single.res$pos.beg[i, acc] != 0){
        # if(i == 2) stop('670')
        pos = single.res$pos.beg[i, acc]:single.res$pos.end[i, acc]
        pos = pos[-c(1, length(pos))]
          
        v.aln[fp.single[[i]]] = pos
        # if(length(unique(v.aln)) != (sum(v.aln != 0) + 1)) stop('1')
      } 
    }
    if(length(unique(v.aln)) != (sum(v.aln != 0) + 1)) stop('1: Duplicated positions in Singletons')
    
    # Add short
    for(i in 1:length(msa.res$len)){
      if(acc %in% colnames(msa.res$aln[[i]])){
        
        v.aln[fp.short[[i]]] = msa.res$aln[[i]][,acc]
      } 
    }
    if(length(unique(v.aln)) != (sum(v.aln != 0) + 1)) stop('2: Duplicated positions in short alignments')
    
    # add long
    for(i in 1:length(mafft.aln.pos)){
      if(acc %in% rownames(mafft.aln.pos[[i]])){
        v.aln[fp.long[[i]]] = mafft.aln.pos[[i]][acc,]
        
        # if(7916909 %in%  mafft.aln.pos[[i]][acc,]) stop(i)
        # if(length(unique(v.aln)) != (sum(v.aln != 0) + 1)) stop('3')
      } 
    }
    if(length(unique(v.aln)) != (sum(v.aln != 0) + 1)) stop('3: Duplicated positions in long alignments')
    
    # Maybe something was overlapped by accident
    
    suppressMessages({
      h5write(v.aln, file.res, paste(gr.accs.e, acc, sep = ''))
      
      stat.comb <- rbind(stat.comb, 
                         data.frame(comb = acc, coverage = sum(v.aln != 0)))
      
    })
    
  }
  H5close()
  gc()
  
}

warnings()


saveRDS(stat.comb, paste(path.cons, 'stat', s.comb, '_', ref.pref,'.rds', sep = ''))



# ***********************************************************************
# ---- Manual testing ----


if(F){

  # ********************
  # Help testing
  file.comb = 'res_1_1_ref_0.h5'
  file.comb = 'msa_1_1_ref_0.h5'
  
  library(rhdf5)
  gr.accs.b <- "/accs"
  gr.accs.e <- "accs/"
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  v = c()
  for(acc in accessions){
    v.acc = h5read(file.comb, paste(gr.accs.e, acc, sep = ''))
    v = cbind(v, v.acc)
    print(acc)
    print(sum((v.acc != 0) & (!is.na(v.acc))))
  }
  
  # ********************
  # Tom-Greg
  library(rhdf5)
  source('../../../pannagram/utils/utils.R')
  path.cons = './'
  path.chromosomes = '/home/anna/storage/arabidopsis/pacbio/pan_test/tom/chromosomes/'
  ref.pref = '0'
  path.mafft.out = '../mafft_out/'
  n.flank = 30
  
  gr.accs.e <- "accs/"
  gr.accs.b <- "/accs"
  gr.break.e = 'break/'
  gr.break.b = '/break'
  
  max.block.elemnt = 3 * 10^ 6
  
  # ********************
  # Rhizobia
  library(rhdf5)
  source('../../../arabidopsis/pacbio/pannagram/utils/utils.R')
  path.cons = './'
  path.chromosomes = '../chromosomes/'
  ref.pref = 'GR3013_prokka'
  path.mafft.out = '../mafft_out/'
  n.flank = 30
  
  gr.accs.e <- "accs/"
  gr.accs.b <- "/accs"
  gr.break.e = 'break/'
  gr.break.b = '/break'
  
  max.block.elemnt = 3 * 10^ 6
  
  
}