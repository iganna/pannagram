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
source(system.file("pangen/comb_func.R", package = "pannagram"))
# source("synteny_funcs.R")

# pokazStage('Step 10. Prepare sequences for MAFFT')

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.cons",        type = "character", default = NULL, help = "Path to consensus directory"),
  make_option("--path.chromosomes", type = "character", default = NULL, help = "Path to directory with chromosomes"),
  make_option("--path.mafft.in",    type = "character", default = NULL, help = "Path to directory where to combine fasta files for mafft runs"),
  
  make_option("--max.len.gap",      type = "integer",   default = NULL, help = "Max length of the gap"),
  
  make_option("--cores",            type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",         type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",        type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

#TODO: SHOULD BE PARAMATERS
len.short = 50
# len.large = 40000
n.flank = 30
len.large.mafft = 15000

s.flank.beg = rep('A', n.flank)
s.flank.end = rep('T', n.flank)

# print(opt)

# ***********************************************************************
# ---- Logging ----
source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----
source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

aln.type.in = aln.type.clean

# ***********************************************************************
# ---- Values of parameters ----

# Max len gap
if (is.null(opt$max.len.gap)) {
  stop("Error: max.len.gap is NULL")
} else {
  len.large <- opt$max.len.gap
}

# Number of cores for parallel processing
num.cores = opt$cores
if(is.null(num.cores)) stop('Whong number of cores: NULL')

pokaz('Number of cores', num.cores)
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

s.pattern <- paste0("^", aln.type.in, ".*")
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.in, "", files)
pref.combinations <- sub(".h5", "", pref.combinations)

if(length(pref.combinations) == 0) {
  stop('No files with the ref-based alignments are found')
}

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

echo = T
for(s.comb in pref.combinations){
  
  if(echo) pokaz('* Combination', s.comb)
  q.chr = strsplit(s.comb, '_')[[1]][1]
  
  file.ws = paste0(path.cons, 'breaks_ws_', s.comb, '.RData')
  load(file.ws)
  
  # Define indexes for short and singletons
  idx.singl = which((breaks$single == 1) & ((breaks$idx.end - breaks$idx.beg - 1) == 0))
  idx.short = which((breaks$single != 1) & (breaks$len.acc <= len.short))
  
  # Define indexes for long sequences
  idx.large = which((breaks$single != 1) & (breaks$len.acc > len.short) & (breaks$len.mean <= len.large.mafft))
  idx.extra = which((breaks$single != 1) & (breaks$len.mean > len.large.mafft))
  
  pokaz('Lengths singl/short/large/extra',
        length(idx.singl), 
        length(idx.short), 
        length(idx.large), 
        length(idx.extra))
  
  # ----
  
  # IMPORTANT: THERE ARE SOME BREAKS WITH LOOK LIKE SINGLETONS< BUT THEY ARE NOT
  # if(sum(length(idx.singl) +
  #        length(idx.short) +
  #        length(idx.large) +
  #        length(idx.extra)) != nrow(breaks)) {
  #   save(list = ls(), file = "tmp_wrong_Chrckpoint7.RData")
  #   stop('Chrckpoint7')
  # } 
  
  # Names of files
  n.digits <- nchar(as.character(nrow(breaks)))
  format.digits <- paste0("%0", n.digits, "d")
  breaks$file = paste0('Gap_', s.comb,  '_',
                       sprintf(format.digits, 1:nrow(breaks)), '_',
                       breaks$idx.beg, '_',
                       breaks$idx.end, '_flank_', n.flank, '.fasta')
  
  breaks$file[idx.extra] = sub('.fasta', '_extra.fasta', breaks$file[idx.extra])
  
  # Save breaks
  file.breaks.merged = paste0(path.cons, 'breaks_merged_', s.comb,'.rds')
  saveRDS(breaks, file.breaks.merged)
  
  ## ---- Save singletons ----
  saveRDS(list(pos.beg = v.beg[idx.singl,],
               pos.end = v.end[idx.singl,],
               ref.pos = data.frame(beg = breaks$idx.beg[idx.singl],
                                    end = breaks$idx.end[idx.singl]) ), 
          paste0(path.cons, 'singletons_',s.comb,'.rds'), compress = F)
  
  # if(s.comb == '10_10'){
  #   save(list = ls(), file = "tmp_workspace.RData")
  #   stop('Enough..')
  # }
  
  ## ---- Analyse by portions ----
  
  idx.remained = setdiff(1:nrow(breaks), idx.singl)
  
  aln.seqs <- vector("list", length = nrow(breaks))
  aln.seqs.names <- vector("list", length = nrow(breaks))
  
  k = 10
  order.acc = ceiling(1:length(accessions) / k)
  for(i.k in min(order.acc):max(order.acc)){
    
    accessions.tmp = accessions[which(order.acc == i.k)]
    for(acc in accessions.tmp){
      pokaz(acc)
      file.chromosome = paste(path.chromosomes, 
                              acc, 
                              '_chr', q.chr, '.fasta', sep = '')
      genome = readFasta(file.chromosome)
      genome = seq2nt(genome)
      
      getSeq <- function(p1, p2, for.mafft = F){
        
        if(p1 > 0){
          s.strand = '+'
          p1 = p1 + 1
          p2 = p2 - 1
          if(p2 < p1) {
            stop(paste('Wrong direction in strand (+) in row', irow))
          }
          seq = genome[p1:p2]
          pos = p1:p2
        } else {
          s.strand = '-'
          tmp = p1
          p1 = -p2 + 1
          p2 = -tmp - 1
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
        
        return(seq = seq)
      }
      
      # Get only those breaks where the accession is not empty
      p1 = v.beg[idx.remained, acc]
      p2 = v.end[idx.remained, acc]
      idx.acc = (p1 != 0) & (p2 != 0)
      idx.tmp.acc = idx.remained[idx.acc]
      p1 = p1[idx.acc]
      p2 = p2[idx.acc]
      
      # Get sequences
      subsets <- mapply(function(b, e, for.mafft) getSeq(b, e, for.mafft), p1, p2, idx.tmp.acc %in% c(idx.large, idx.extra))
      
      if(names(subsets) == 'GCA_040115645.1.GCA_040115645.1|10|1384970|1385000|+|31'){
        save(list = ls(), file = "tmp_workspace.RData")
        stop()
      }
      
      # Save sequences
      pokaz('---')
      pokaz(names(subsets))
      aln.seqs[idx.tmp.acc] <- mapply(function(x, y) c(x, y), aln.seqs[idx.tmp.acc], subsets, SIMPLIFY = FALSE)
      aln.seqs.names[idx.tmp.acc] <- mapply(function(x, y) c(x, y), aln.seqs.names[idx.tmp.acc], names(subsets), SIMPLIFY = FALSE)
      pokaz(aln.seqs.names[idx.tmp.acc])
      
      rm(genome)
    }
    idx.save = c(idx.large, idx.extra)
    
    if(length(idx.save) == 0) next
    
    idx.save = idx.save[sapply(idx.save, function(x) !is.null(aln.seqs[[x]]))]
    
    if(length(intersect(idx.save, idx.short)) > 0) stop('Wrong idx are saved')
    
    if(num.cores == 1){
      pokaz('Save sequences...')
      for(i in idx.save){
        writeFasta(aln.seqs[[i]], 
                   file = paste0(path.mafft.in,breaks$file[i]), 
                   seq.names = aln.seqs.names[[i]],
                   append = T)
      }
    } else { # Many cores
      pokaz('Save sequences with parallel...')
      foreach(i = idx.save,
              .packages=c('crayon'))  %dopar% {
                writeFasta(aln.seqs[[i]], 
                           file = paste0(path.mafft.in,breaks$file[i]), 
                           seq.names = aln.seqs.names[[i]],
                           append = T)
              }
    }
    if(echo) pokaz('.. done!')
    
    aln.seqs[idx.save] <- list(NULL)
    aln.seqs.names[idx.save] <- list(NULL)
    
    n.null <- sum(sapply(aln.seqs, Negate(is.null)))
    
  }
  
    
  # save(list = ls(), file = "tmp_workspace.RData")

  for(i in setdiff(1:length(aln.seqs), idx.short)){
    aln.seqs[i] <- list(NULL)
  }
  
  n.null <- sum(sapply(aln.seqs, Negate(is.null)))
  pokaz(n.null, length(idx.short))
  if(n.null != length(idx.short)) {
    # pokazAttention('fix short')
    # save(list = ls(), file = "tmp_wrong_number_of_short.RData")
    stop('Wrong number of short')
  }
  
  all.local.objects <- c("breaks", "aln.seqs", "aln.seqs.names", "idx.short", "accessions")
  file.ws = paste0(path.cons, 'small_ws_', s.comb, '.RData')
  save(list = all.local.objects, file = file.ws)
  
  rm(aln.seqs)
  rm(aln.seqs.names)

}


if(num.cores > 1){
  stopCluster(myCluster)
}

warnings()


# ***********************************************************************
# ---- Manual testing ----
