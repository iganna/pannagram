# Combine all alignments together into the final one

suppressMessages({
library(rhdf5)
library('foreach')
library(doParallel)
library("optparse")
source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func_mafft_refine.R", package = "pannagram"))
})

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.mafft.out"), type="character", default=NULL, 
              help="path to directory, where to mafft results are", metavar="character"),
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to directory with the consensus", metavar="character"),
  make_option(c("--path.out"), type="character", default=NULL, 
              help="path to directory with MSA output", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("--path.log"), type = "character", default = NULL,
              help = "Path for log files", metavar = "character"),
  make_option(c("--log.level"), type = "character", default = NULL,
              help = "Level of log to be shown on the screen", metavar = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# TODO: SHOULD BE PARAMETERS

n.flank = 30
max.block.elemnt = 3 * 10^ 6

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores.max = 10
num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))

if (!is.null(opt$path.mafft.out)) path.mafft.out <- opt$path.mafft.out
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if (!is.null(opt$path.out)) path.out <- opt$path.out


# ***********************************************************************

# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type.comb, ".*")
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.comb, "", files)
pref.combinations <- sub(".h5", "", pref.combinations)

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)


# ***********************************************************************
# ---- MAIN program body ----

stat.comb <- data.frame(comb = character(),
                        coverage = numeric(),
                        stringsAsFactors = FALSE)


for(s.comb in pref.combinations){
  
  pokaz('* Combination', s.comb, file=file.log.main, echo=echo.main)
  
  # Get accessions
  file.comb = paste0(path.cons, aln.type.comb, s.comb,'.h5')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)
  base.len = length(h5read(file.comb, paste0(gr.accs.e, accessions[1])))
  
  
  # n.rows.block = round(max.block.elemnt / n.acc)
  # seq.blocks = seq(1, base.len - n.rows.block, n.rows.block)
  
  # ---- All MAFFT results for the combination ----
  pref = paste('Gap', s.comb, sep = '_')
  mafft.res = data.frame(file = list.files(path = path.mafft.out, 
                                           pattern = paste0('^', pref, '.*_flank_', n.flank, '.*_aligned.*\\.fasta$')))
  
  
  y = lapply(mafft.res$file, function(s) strsplit(s, '_')[[1]][1:6])
  z = t(matrix(unlist(y), ncol = length(y)))
  mafft.res$comb = paste(z[,2], z[,3], sep = '_')
  mafft.res$beg = as.numeric(z[,5])
  mafft.res$end = as.numeric(z[,6])
  mafft.res$id = as.numeric(z[,3])
  
  
  # ---- Get long alignment positions ----
  pokaz('Read Long alignments..', file=file.log.main, echo=echo.main)
  idx.skip = c()
  mafft.aln.pos = list()
  
  pokaz('Number of mafft files', nrow(mafft.res), file=file.log.main, echo=echo.main)
  for(i in 1:nrow(mafft.res)){
    

    pokaz('Aln', i, mafft.res$file[i], file=file.log.main, echo=echo.main)
    
    file.aln = paste0(path.mafft.out, mafft.res$file[i])
    
    if(!file.exists(file.aln)) {
      idx.skip = c(idx.skip, i)
      next
    }

    if (grepl("_aligned2", mafft.res$file[i])) {
      remove.flank = F
    } else {
      remove.flank = T
    }
    
    aln.seq = readFasta(file.aln)
    
    # ---
    # WITHOUT REFINEMENT
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
      
      if(remove.flank){
        tmp.nongap = tmp.nongap[-(1:(n.flank))]
        tmp.nongap <- tmp.nongap[1:(length(tmp.nongap) - n.flank)]  
      }
      
      p1 = pos.aln[1, i.seq]
      p2 = pos.aln[2, i.seq]
      pos.tmp = p1:p2
      pos.mx[i.seq, tmp.nongap] = pos.tmp
    }
    # aln.mx = aln.mx[,colSums(aln.mx != '-') != 0]
    pos.mx = pos.mx[,colSums(pos.mx != 0) != 0, drop=F]
    row.names(pos.mx) = name.acc


    mafft.aln.pos[[mafft.res$file[i]]] = pos.mx
    # ---

    
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
  pokaz('Read Short alignments..', file=file.log.main, echo=echo.main)
  msa.res = readRDS(paste0(path.cons, 'aln_short_', s.comb, '.rds'))
  msa.res$len = unlist(lapply(msa.res$aln, nrow))
  msa.res$extra = msa.res$len - (msa.res$ref.pos$end - msa.res$ref.pos$beg - 1)
  # if(min(msa.res$extra) < 0) stop('Short: Wrong lengths of alignment and gaps')
  msa.res$extra[msa.res$extra < 0] = 0

  # ---- Singletons alignments ----
  pokaz('Read Singletons..', file=file.log.main, echo=echo.main)
  single.res = readRDS(paste0(path.cons, 'singletons_', s.comb, '.rds'))
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
    if(length(fp.short[[i]]) != nrow(msa.res$aln[[i]])){
      
      
      file.ws = "tmp_workspace.RData"
      all.local.objects <- ls()
      save(list = all.local.objects, file = file.ws)
      
      stop(paste0('Short', i)) 
    }
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
    # pokaz(  sum(pos.delete == 1), sum(pos.delete == -1), file=file.log.main, echo=echo.main)
    pos.delete[pos.end] = pos.delete[pos.end]-1
    # pokaz(  sum(pos.delete == 1), sum(pos.delete == -1), file=file.log.main, echo=echo.main)
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
  
  file.res = paste0(path.out, aln.type.msa, s.comb,'.h5')
  if (file.exists(file.res)) file.remove(file.res)
  h5createFile(file.res)
  
  suppressMessages({
    h5createGroup(file.res, gr.accs.e)
  })
  
  pos.nonzero = fp.main != 0
  for(acc in accessions){
    pokaz('Accession', acc, file=file.log.main, echo=echo.main)
    v = h5read(file.comb, paste0(gr.accs.e, acc))
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
        
        # if(1835154 %in%  mafft.aln.pos[[i]][acc,]) stop(i)
        # if(length(unique(v.aln)) != (sum(v.aln != 0) + 1)) stop('3')
      } 
    }
    if(length(unique(v.aln)) != (sum(v.aln != 0) + 1)) stop('3: Duplicated positions in long alignments')
    
    # Maybe something was overlapped by accident
    
    suppressMessages({
      h5write(v.aln, file.res, paste0(gr.accs.e, acc))
      
      stat.comb <- rbind(stat.comb, 
                         data.frame(comb = acc, coverage = sum(v.aln != 0)))
      
    })
    
  }
  H5close()
  gc()
  
}

warnings()

saveRDS(stat.comb, paste0(path.cons, 'stat', s.comb,'.rds'))


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
    v.acc = h5read(file.comb, paste0(gr.accs.e, acc))
    v = cbind(v, v.acc)
    print(acc)
    print(sum((v.acc != 0) & (!is.na(v.acc))))
  }
  
  # ********************
  # Tom-Greg
  library(rhdf5)
source(system.file("pannagram/utils/utils.R", package = "pannagram"))
  path.cons = './'
  path.chromosomes = '/home/anna/storage/arabidopsis/pacbio/pan_test/tom/chromosomes/'
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
source(system.file("utils/utils.R", package = "pannagram"))
  path.cons = './'
  path.chromosomes = '../chromosomes/'
  path.mafft.out = '../mafft_out/'
  n.flank = 30
  
  gr.accs.e <- "accs/"
  gr.accs.b <- "/accs"
  gr.break.e = 'break/'
  gr.break.b = '/break'
  
  max.block.elemnt = 3 * 10^ 6
  
  
}

pokaz('Done.', file=file.log.main, echo=echo.main)


