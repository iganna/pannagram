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
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("--path.log"), type = "character", default = NULL,
              help = "Path for log files", metavar = "character"),
  make_option(c("--log.level"), type = "character", default = NULL,
              help = "Level of log to be shown on the screen", metavar = "character"),
  make_option(c("--max.len.gap"), type = "integer", default = NULL,
              help = "Max len of the gap", metavar = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);


#TODO: SHOULD BE PARAMATERS
len.short = 50
# len.large = 40000
n.flank = 30

# print(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

# ***********************************************************************
# ---- Values of parameters ----

# Max len gap
if (is.null(opt$max.len.gap)) {
  stop("Error: max.len.gap is NULL")
} else {
  len.large <- opt$max.len.gap
}

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

# Path to mafft input
if (!is.null(opt$path.mafft.in)) path.mafft.in <- opt$path.mafft.in

# **************************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type.comb, ".*")
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.comb, "", files)
pref.combinations <- sub(".h5", "", pref.combinations)

if(length(pref.combinations) == 0) {
  stop('No files with the ref-based alignments are found')
}

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

echo = F
for(s.comb in pref.combinations){
  
  if(echo) pokaz('* Combination', s.comb)
  q.chr = strsplit(s.comb, '_')[[1]][1]
  
  file.comb = paste0(path.cons, aln.type.comb, s.comb,'.h5')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)
  
  # ---- Read Breaks  ----
  file.breaks = paste0(path.cons, 'breaks_', s.comb,'.rds')
  idx.breaks = readRDS(file.breaks)
  
  # ---- Merge coverages ----
  n.init = nrow(idx.breaks)
  idx.breaks = idx.breaks[,c('idx.beg', 'idx.end')]
  
  idx.breaks = idx.breaks[order(-idx.breaks$idx.end),]
  idx.breaks = idx.breaks[order(idx.breaks$idx.beg),]
  
  idx.breaks$id = 1:nrow(idx.breaks)
  idx.breaks = idx.breaks[!duplicated(idx.breaks[, c('idx.beg', 'idx.end')]),]
  idx.breaks$cnt = c(idx.breaks$id[-1], n.init+1) - idx.breaks$id
  if(!(sum(idx.breaks$cnt) == n.init)) stop('Checkpoint1')
  if (any(idx.breaks$cnt < 0)) stop('Checkpoint2')
  
  # 
  # pos = rep(0, 40000000)
  # for(irow in 1:nrow(idx.breaks)){
  #   pos[idx.breaks$idx.beg[irow]:idx.breaks$idx.end[irow]] = irow
  # }
  # sum(pos == 0)
  
  idx.breaks <- idx.breaks[order(-idx.breaks$idx.end), ]
  idx.breaks <- idx.breaks[order(idx.breaks$idx.beg), ]
  
  # Merge overlaps
  n = 0
  while (n != nrow(idx.breaks)) {
    n = nrow(idx.breaks)
    
    idx_full_cover = which(idx.breaks$idx.beg[-1] <= idx.breaks$idx.end[-nrow(idx.breaks)])
    idx_full_cover = setdiff(idx_full_cover, idx_full_cover + 1)
    
    if (length(idx_full_cover) == 0) break
    
    idx.breaks$cnt[idx_full_cover] = idx.breaks$cnt[idx_full_cover] + idx.breaks$cnt[idx_full_cover + 1]
    idx.breaks$idx.end[idx_full_cover] = pmax(idx.breaks$idx.end[idx_full_cover], idx.breaks$idx.end[idx_full_cover + 1])
    
    idx.breaks = idx.breaks[-(idx_full_cover + 1), ]
  }
  
  # Length of gaps
  idx.breaks$idx = paste('gap', 1:nrow(idx.breaks), sep = '|')
  idx.breaks$len = idx.breaks$idx.end - idx.breaks$idx.beg + 1
  
  if(sum(idx.breaks$cnt) != n.init) stop('Checkpoint3')
  
  
  ## ---- Get begin-end positions of gaps ----
  v.beg = c()
  v.end = c()
  for(acc in accessions){
    pokaz(acc)
    
    x.acc = h5read(file.comb, paste0(gr.accs.e, acc))
    b.acc = h5read(file.comb, paste0(gr.blocks, acc))
    
    x.beg = fillPrev(x.acc)[idx.breaks$idx.beg]
    x.end = fillNext(x.acc)[idx.breaks$idx.end]
    
    idx.no.zero = (x.beg != 0) & (x.end != 0)
    idx.no.zero[idx.no.zero] = b.acc[abs(x.beg[idx.no.zero])] == b.acc[abs(x.end[idx.no.zero])]
  
    x.beg[!idx.no.zero] = 0
    x.end[!idx.no.zero] = 0
    
    v.beg = cbind(v.beg, x.beg)
    v.end = cbind(v.end, x.end)
    
  }
  colnames(v.beg) = accessions
  colnames(v.end) = accessions
  
  # Check inversions
  if (any(sign(v.beg * v.end) < 0)) stop('Checkpoint4')
  # To test: which(rowSums(sign(v.beg * v.end) < 0) > 0)
  
  # Check direction
  if (any(sign(v.end - v.beg) < 0)) stop('Checkpoint5')
  
  # ---- Zero-positions mask ----

  zero.mask = (v.end == 0) | (v.beg == 0)
  v.end[zero.mask] = 0
  v.beg[zero.mask] = 0
  
  # ---- Check lengths ----
  v.len = v.end - v.beg - 1
  v.len[zero.mask] = 0
  
  if (any(v.len < 0)) stop('Checkpoint6')
  v.len[zero.mask] = 0
  
  zero.len.mask = (v.len == 0)
  v.end[zero.len.mask] = 0
  v.beg[zero.len.mask] = 0
  
  if(echo) pokaz('Amount of long fragments', sum(v.len > len.large))
  
  # ---- Checkups for duplicates ----
  for(icol in 1:ncol(v.len)){
    idx.dup = unique(v.beg[duplicated(v.beg[,icol]),icol])
    if(length(setdiff(idx.dup, 0)) != 0) {
      stop(paste('Duplicated in column', icol, 'in v.beg, amount:', length(idx.dup) - 1))
    }
    # v.len[v.beg[,icol] %in% idx.dup, icol] = 0
    
    idx.dup = unique(v.end[duplicated(v.end[,icol]),icol])
    if(length(setdiff(idx.dup, 0)) != 0) {
      stop(paste('Duplicated in column', icol, 'in v.end, amount:', length(idx.dup) - 1))
    }
    # v.len[v.end[,icol] %in% idx.dup, icol] = 0
  }

  
  # ---- Subdivide into categories ----
  
  idx.breaks$single = rowSums(v.len != 0)
  idx.breaks$len.acc = rowMax(v.len)
  idx.singl = which(idx.breaks$single == 1)
  idx.short = which((idx.breaks$single != 1) & (idx.breaks$len.acc <= len.short))
  idx.large = which((idx.breaks$single != 1) & (idx.breaks$len.acc > len.short) & (idx.breaks$len.acc <= len.large))
  idx.extra = which((idx.breaks$single != 1) & (idx.breaks$len.acc > len.large))
  
  if(sum(length(idx.singl) + 
         length(idx.short) + 
         length(idx.large) + 
         length(idx.extra)) != nrow(idx.breaks)) stop('Chrckpoint7')
  
  # Names of files
  n.digits <- nchar(as.character(nrow(idx.breaks)))
  format.digits <- paste0("%0", n.digits, "d")
  idx.breaks$file = paste0('Gap_', s.comb,  '_',
                  sprintf(format.digits, 1:nrow(idx.breaks)), '_',
                  idx.breaks$idx.beg, '_',
                  idx.breaks$idx.end, '_flank_', n.flank, '.fasta')
  
  # Save breaks
  file.breaks.merged = paste0(path.cons, 'breaks_merged_', s.comb,'.rds')
  saveRDS(idx.breaks, file.breaks.merged)
  
  ## ---- Save singletons ----
  saveRDS(list(pos.beg = v.beg[idx.singl,],
               pos.end = v.end[idx.singl,],
               ref.pos = idx.breaks[idx.singl, c('idx.beg', 'idx.end')]), paste0(path.cons, 'singletons_',s.comb,'.rds'), compress = F)
  
  ## ---- Save Long and Keep short ----
  aln.seqs <- vector("list", length = nrow(idx.breaks))
  aln.pos <- vector("list", length = nrow(idx.breaks))
  for(acc in accessions){
    
    # Read the chromosome
    file.chromosome = paste(path.chromosomes, 
                            acc, 
                            '_chr', q.chr, '.fasta', sep = '')
    genome = readFasta(file.chromosome)
    genome = seq2nt(genome)
    
    ### ---- Get sequences ----
    s.flank.beg = rep('A', n.flank)
    s.flank.end = rep('T', n.flank)
    
    getSeq <- function(irow, for.mafft = F){
      
      if(v.beg[irow, acc] > 0){
        s.strand = '+'
        p1 = v.beg[irow, acc] + 1
        p2 = v.end[irow, acc] - 1
        if(p2 < p1) {
          stop(paste('Wrong direction in strand (+) in row', irow))
        }
        seq = genome[p1:p2]
        pos = p1:p2
      } else {
        s.strand = '-'
        p1 = -v.end[irow, acc] + 1
        p2 = -v.beg[irow, acc] - 1
        if(p2 < p1)  {
          stop(paste('Wrong direction in strand (-) in row', irow))
        }
        seq = genome[p1:p2]
        seq = revCompl(seq)
        pos = (-p2):(-p1)
      }
      
      if(for.mafft){
        seq = c(s.flank.beg, seq, s.flank.end)
      }
      
      seq = nt2seq(seq)
      
      seq.name = paste(acc, q.chr, pos[1], pos[length(pos)], s.strand, p2 - p1 + 1, sep = '|')  
      names(seq) = seq.name
      
      return(list(seq = seq, pos = pos))
    }
    
    ### ---- MAFFT ----
    for(irow in idx.large){
      if(v.beg[irow, acc] == 0) next
      
      res = getSeq(irow)
      seq = res$seq
      
      # Write to the fasta file
      file.out = paste0(path.mafft.in, idx.breaks$file[irow])
      if(file.exists(file.out)){
        writeFasta(seq, file.out, append = T)
      } else {
        writeFasta(seq, file.out)
      }
      
    }
    
    ### ---- Short sequences ----
    for(irow in idx.short){
      if(v.beg[irow, acc] == 0) next
      
      res = getSeq(irow)
      seq = res$seq
      pos = res$pos

      # Save to the further slignment
      aln.seqs[[irow]][acc] = seq
      aln.pos[[irow]][[acc]] = pos
    }
    
    # ### ---- Extra sequences ----
    # aln.seqs <- vector("list", length = nrow(idx.breaks))
    # aln.pos <- vector("list", length = nrow(idx.breaks))
    # for(irow in idx.large){
    #   seq = getSeq(irow)
    #   
    #   # Save to the further alignment
    #   aln.seqs[[irow]][seq.name] = nt2seq(seq)
    #   aln.pos[[irow]][[seq.name]] = pos
    # }
    
  }

  # ---- Align Short sequences ----
  if(echo) pokaz('Align short seqs')
  
  # Checkuip the numer of sequences in alignments
  tmp = unlist(lapply(aln.seqs, length))
  if(sum(tmp == 1) > 0) stop('Checkup-short1')
  
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
      idx.gap.pos = idx.breaks$idx.beg[irow]
      
      res.msa[[i.tmp]] = CODE_ALN_SHORT()
    }
  } else {
    # Many cores
    pokaz('Parallel computing: short sequences')
    res.msa <- foreach(seqs = aln.seqs[idx.short],
                       pos.idx = aln.pos[idx.short],
                       idx.gap.pos = idx.breaks$idx.beg[idx.short],
                       .packages=c('muscle', 'Biostrings', 'crayon'))  %dopar% {
                         return(CODE_ALN_SHORT()) 
                       }
  }
  
  saveRDS(list(aln = res.msa,
               seqs = aln.seqs[idx.short], 
               pos.idx = aln.pos[idx.short],
               ref.pos = idx.breaks[idx.short, c('idx.beg', 'idx.end')]), paste0(path.cons, 'aln_short_',s.comb,'.rds'), compress = F)
  

  
  # ---- Extra alignments ----
  
  
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
  path.cons = '/Volumes/Samsung_T5/vienn/test/symA_test_0/intermediate/consensus/'
  path.mafft.in = '/Volumes/Samsung_T5/vienn/test/symA_test_0/intermediate/mafft_in/'
  path.chromosomes = '/Volumes/Samsung_T5/vienn/test/symA_test_0/intermediate/chromosomes/'
  
}




