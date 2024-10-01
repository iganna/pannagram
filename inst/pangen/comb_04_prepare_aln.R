# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
  library(muscle)  #BiocManager::install("muscle")
  # library(Biostrings)
})

source(system.file("utils/utils.R", package = "pannagram"))
# source("synteny_funcs.R")

# pokazStage('Step 10. Prepare sequences for MAFFT')

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to consensus directory", metavar="character"),
  make_option(c("--path.chromosomes"), type="character", default=NULL, 
              help="path to directory with chromosomes", metavar="character"),
  make_option(c("--path.mafft.in"), type="character", default=NULL, 
              help="path to directory, where to combine fasta files for mafft runs", metavar="character"),
  make_option(c("-p", "--ref.pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
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

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores.max = 10
num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))
if(num.cores > 1){
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster)
}

# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop('Consensus folder doesn???t exist')

# Path to chromosomes
if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes

# Reference genome
if (is.null(opt$ref.pref)) {
  stop("ref.pref is NULL")
} else {
  ref.pref <- opt$ref.pref
}


# Path to mafft input
if (!is.null(opt$path.mafft.in)) path.mafft.in <- opt$path.mafft.in

# **************************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", 'res_', ".*", '_ref_', ref.pref)
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pokaz('Files', files, file=file.log.main, echo=echo.main)
pref.combinations = gsub("res_", "", files)
pref.combinations <- sub("_ref.*$", "", pref.combinations)

pokaz('Reference:', ref.pref, file=file.log.main, echo=echo.main)
pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ----  Combine correspondence  ----

gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'
max.len.gap = 100000
len.short = 50
n.flank = 30

gr.blocks = 'blocks/'

# ***********************************************************************
# ---- MAIN program body ----

echo = F
for(s.comb in pref.combinations){
  
  if(echo) pokaz('* Combination', s.comb)
  q.chr = strsplit(s.comb, '_')[[1]][1]
  
  file.comb = paste0(path.cons, 'res_', s.comb,'_ref_',ref.pref,'.h5')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)
  
  # # For testing
  # v = c()
  # for(acc in accessions){
  #   v.acc = h5read(file.comb, paste0(gr.accs.e, acc))
  #   v = cbind(v, v.acc)
  # }
  
  # ---- Merge coverages ----
  file.breaks = paste0(path.cons, 'breaks_', s.comb,'_ref_',ref.pref,'.rds')
  idx.break = readRDS(file.breaks)
  # idx.break = idx.break[idx.break$acc != paste0('', ref.pref),]  # Remove the reference correspondence
  
  
  
  # Merge full coverages 
  idx.break = idx.break[order(-idx.break$end),]
  idx.break = idx.break[order(idx.break$beg),]
  
  n = 0
  while(n != nrow(idx.break)){
    n = nrow(idx.break)
    # print(n)
    idx.full.cover = which(diff(idx.break$end) <= 0) + 1
    if(length(idx.full.cover) == 0) break
    idx.break = idx.break[-idx.full.cover,]
  }
  
  # Merge overlaps
  n = 0
  while(n != nrow(idx.break)){
    n = nrow(idx.break)
    # print(n)
    idx.over = which((idx.break$beg[-1] < idx.break$end[-n]))
    idx.over = setdiff(idx.over, idx.over+1)
    
    if(length(idx.over) == 0) break
    
    idx.break$beg[idx.over + 1] = idx.break$beg[idx.over]
    idx.break = idx.break[-idx.over,]
  }
  
  # Length of gaps
  idx.break$idx = paste('gap', 1:nrow(idx.break), sep = '|')
  idx.break$len = idx.break$end - idx.break$beg + 1
  
  
  ## ---- Get begin-end positions of gaps ----
  v.beg = c()
  v.end = c()
  for(acc in accessions){
    pokaz(acc)
    
    # if(sub('', '', acc) ==ref.pref){
    #   v.beg = cbind(v.beg, idx.break$beg)
    #   v.end = cbind(v.end, idx.break$end)
    #   next
    # }
    x.acc = h5read(file.comb, paste0(gr.accs.e, acc))
    blocks.acc = h5read(file.comb, paste0(gr.blocks, acc))
    
    # ### ---- Find prev and next ----
    # n = nrow(x.acc)
    # for(i.tmp in 1:2){
    #   if(echo) pokaz('Accession', acc, 'fill prev/next', i.tmp)
    #   v.acc = x.acc
    #   if(i.tmp == 2) {
    #     v.acc = rev(v.acc)
    #   }
    #   v.rank = rank(v.acc)
    #   v.rank[v.acc == 0] = 0
    #   
    #   v.zero.beg = which((v.acc[-1] == 0) & (v.acc[-n] != 0)) + 1
    #   v.zero.end = which((v.acc[-1] != 0) & (v.acc[-n] == 0)) 
    #   if(v.acc[1] == 0) v.zero.end = v.zero.end[-1]
    #   if(v.acc[n] == 0) v.zero.end = c(v.zero.end, n)
    #   
    #   # ..... WITHIN ONE STRATCH BLOCK .....
    #   idx = which(abs(v.rank[v.zero.beg-1] - v.rank[v.zero.end+1]) != 1)
    #   v.zero.beg = v.zero.beg[-idx]
    #   v.zero.end = v.zero.end[-idx]
    #   # .....
    #   
    #   tmp = rep(0, n)
    #   tmp[v.zero.beg] = v.acc[v.zero.beg-1]
    #   tmp[v.zero.end] = -v.acc[v.zero.beg-1]
    #   tmp[v.zero.end[v.zero.beg == v.zero.end]] = 0
    #   tmp = cumsum(tmp)
    #   tmp[v.zero.end] = v.acc[v.zero.beg-1]
    #   
    #   if(i.tmp == 1){
    #     v.prev = x.acc + tmp
    #   } else {
    #     tmp = rev(tmp)
    #     v.next = x.acc + tmp
    #   }
    # }
    # 
    # v.beg = cbind(v.beg, v.prev[idx.break$beg] * (v.prev[idx.break$beg+1] != 0))
    # v.end = cbind(v.end, v.next[idx.break$end] * (v.next[idx.break$end-1] != 0))
    # v = cbind(v, v.acc)
    
    b.acc.beg = blocks.acc[idx.break$beg]
    b.acc.end = blocks.acc[idx.break$end]
    idx.same.bl = (b.acc.beg == b.acc.end) * 1
    
    
    # v.beg = cbind(v.beg, x.acc[idx.break$beg] )
    # v.end = cbind(v.end, x.acc[idx.break$end] )
    
    v.beg = cbind(v.beg, fillPrev(x.acc)[idx.break$beg] * idx.same.bl )
    v.end = cbind(v.end, fillNext(x.acc)[idx.break$end] * idx.same.bl )
    # 
    
    # print(cbind(v.beg[!idx.same.bl,ncol(v.beg)], v.end[!idx.same.bl,ncol(v.beg)]))
    
  }
  colnames(v.beg) = accessions
  colnames(v.end) = accessions
  
  ## ---- Check lengths ----
  v.len = v.end - v.beg -1
  v.len[v.end == 0] = 0
  v.len[v.beg == 0] = 0
  
  v.len[sign(v.beg * v.end) < 0] = 0
  
  if(echo) pokaz('Amount of long fragments', sum(v.len > max.len.gap))
  v.len[v.len > max.len.gap] = 0
  for(icol in 1:ncol(v.len)){
    idx.dup = v.beg[duplicated(v.beg[,icol]),icol]
    if(echo) pokaz('Duplicated in column', icol, ', amount:', sum(v.len[v.beg[,icol] %in% idx.dup, icol] != 0))
    v.len[v.beg[,icol] %in% idx.dup, icol] = 0
    idx.dup = v.end[duplicated(v.end[,icol]),icol]
    if(echo) pokaz('Duplicated in column', icol, ', amount:', sum(v.len[v.end[,icol] %in% idx.dup, icol] != 0))
    v.len[v.end[,icol] %in% idx.dup, icol] = 0
  }
  v.beg[v.len == 0] = 0
  v.end[v.len == 0] = 0
  
  # ---- Check direction ----
  idx.wrong.dir = sign(v.end - v.beg) < 0
  v.end[idx.wrong.dir] = 0
  v.beg[idx.wrong.dir] = 0
  
  # ---- Get sequences for the alignment, but not singletons ----
  idx.singletons = which(rowSums(v.len != 0) == 1)  # don't have to align
  idx.aln = which(rowSums(v.len != 0) > 1)  # have to align
  
  # Distinguish long (MAFFT) and short(muscle) alignments
  max.len = apply(v.len, 1, max)
  idx.add.flank = max.len > len.short
  s.flank.beg = rep('A', n.flank)
  s.flank.end = rep('T', n.flank)
  
  # Names of breaks
  n.digits <- nchar(as.character(nrow(idx.break)))
  format.digits <- paste0("%0", n.digits, "d")
  s.break = paste('Gap', s.comb, 
                  sprintf(format.digits, 1:nrow(idx.break)), 
                  idx.break$beg, idx.break$end,
                  sep = '_')
  
  # Get sequences
  aln.seqs <- vector("list", length = nrow(idx.break))
  aln.pos <- vector("list", length = nrow(idx.break))
  for(acc in accessions){
    if(echo) pokaz('Accession', acc)
    # read the genome and chromosome
    file.chromosome = paste(path.chromosomes, 
                            acc, 
                            '_chr', q.chr, '.fasta', sep = '')
    genome = readFastaMy(file.chromosome)
    genome = seq2nt(genome)
    
    for(irow in idx.aln){
      if(v.beg[irow, acc] == 0) next
      
      if(v.beg[irow, acc] > 0){
        s.strand = '+'
        p1 = v.beg[irow, acc] + 1
        p2 = v.end[irow, acc] - 1
        if(p2 < p1) {
          print(irow)
          stop('wrong direction, +')
        }
        seq = genome[p1:p2]
        pos = p1:p2
      } else {
        s.strand = '-'
        p1 = -v.end[irow, acc] + 1
        p2 = -v.beg[irow, acc] - 1
        if(p2 < p1)  {
          print(irow)
          stop('wrong direction, -')
        }
        seq = genome[p1:p2]
        seq = revCompl(seq)
        pos = (-p2):(-p1)
      }
      
      # Add flanking sequence, if to MAFFT
      if(idx.add.flank[irow]){
        seq = c(s.flank.beg, seq, s.flank.end)
        seq.name = paste(acc, q.chr, pos[1], pos[length(pos)], s.strand, p2 - p1 + 1, sep = '|')
      } else {
        seq.name = acc
      }
      
      # seq.name = paste(s.break[irow], acc, q.chr, p1, p2, s.strand, p2 - p1 + 1, sep = '|')
      # seq.name = acc
      aln.seqs[[irow]][seq.name] = nt2seq(seq)
      aln.pos[[irow]][[seq.name]] = pos
    }
  }
  
  
  # ---- Align short ----
  idx.short = idx.aln[!idx.add.flank[idx.aln]]
  
  
  if(echo) pokaz('Align short seqs')
  
  
  # Core core for the short alignmgnets
  CODE_ALN_SHORT <- function(echo=F){
    
    set = DNAStringSet(seqs)
    aln = muscle(set, quiet = T)
    
    set = as.character(aln)
    n.pos = nchar(set[1])
    
    if(echo & (n.pos > 10)){
      if(echo) pokaz('Iteration', irow)
      print(aln)
    }
    
    val.acc.pos = matrix(0, nrow=n.pos, ncol=n.acc)
    colnames(val.acc.pos) = accessions
    for(s.acc in names(set)){
      s = strsplit(set[[s.acc]],'')[[1]]
      val.acc.pos[s != '-',s.acc] = pos.idx[[s.acc]]
    }
    
    rm(aln)
    rm(set)
    gc()
    
    return(val.acc.pos)
  }
  
  # Two possible loops depending on the number of cores
  if(num.cores == 1){
    pokaz('No parallel computing: short sequences')
    # One core
    res.msa = list()
    for(i.tmp in 1:length(idx.short)){
      irow = idx.short[i.tmp]
      seqs = aln.seqs[[irow]]
      pos.idx = aln.pos[[irow]]
      idx.gap.pos = idx.break$beg[irow]
      
      res.msa[[irow]] = CODE_ALN_SHORT() # --- COMMON CODE, one core ----
    }
    res.msa = res.msa[idx.short]
  } else {
    # Many cores
    pokaz('Parallel computing: short sequences')
    res.msa <- foreach(seqs = aln.seqs[idx.short],
                       pos.idx = aln.pos[idx.short],
                       idx.gap.pos = idx.break$beg[idx.short],
                       .packages=c('muscle', 'Biostrings', 'crayon'))  %dopar% {
                         
                         return(CODE_ALN_SHORT()) # --- COMMON CODE, many cores ----
                         
                       }
  }
  
  saveRDS(list(aln = res.msa,
               seqs = aln.seqs[idx.short], 
               pos.idx = aln.pos[idx.short],
               ref.pos = idx.break[idx.short, c('beg', 'end')]), paste0(path.cons, 'aln_short_',s.comb,'.rds'), compress = F)
  
  saveRDS(list(pos.beg = v.beg[idx.singletons,],
               pos.end = v.end[idx.singletons,],
               ref.pos = idx.break[idx.singletons, c('beg', 'end')]), paste0(path.cons, 'singletons_',s.comb,'.rds'), compress = F)
  
  # ---- Create files for mafft ----
  
  if(echo) pokaz('Prepare seqs for mafft')
  
  idx.long = idx.aln[idx.add.flank[idx.aln]]
  
  
  
  if(num.cores == 1){
    pokaz('No parallel computing: long sequences')
    for (i.tmp in idx.long) {
      seqs <- aln.seqs[[i.tmp]]
      pos.idx <- aln.pos[[i.tmp]]
      s.aln <- s.break[[i.tmp]]
      
      f.pref <- paste0(path.mafft.in, s.aln, '_flank_', n.flank, '.fasta')
      
      writeFastaMy(seqs, f.pref)
    }
    
  } else {
    pokaz('Parallel computing: long sequences')
    tmp <- foreach(seqs = aln.seqs[idx.long], 
                       pos.idx = aln.pos[idx.long], 
                       s.aln = s.break[idx.long],
                       .packages=c('muscle'))  %dopar% {
                         
                         f.pref = paste(path.mafft.in, s.aln,
                                        '_flank_', n.flank,'.fasta', sep = '')
                         writeFastaMy(seqs, f.pref)
                       }
  }
  H5close()
  gc()
}


if(num.cores > 1){
  stopCluster(myCluster)
}

warnings()


# ***********************************************************************
# ---- Manual testing ----

if(F){
  library(rhdf5)
source(system.file("utils/utils.R", package = "pannagram"))
  path.cons = './'
  path.chromosomes = '/home/anna/storage/arabidopsis/pacbio/pan_test/p27/chromosomes/'
  ref.pref = '0'
  
  
  library(rhdf5)
source(system.file("/Users/annaigolkina/Library/CloudStorage/OneDrive-Personal/vienn/pacbio/pannagram/utils/utils.R", package = "pannagram"))
  path.cons = '/Volumes/Samsung_T5/vienn/alignment/new/consensus/'
  path.chromosomes = '/Volumes/Samsung_T5/vienn/pb_chromosomes/'
  ref.pref = '0'
  
  
}




