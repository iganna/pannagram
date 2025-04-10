suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(pryr)
})


source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/synteny_func.R", package = "pannagram"))
source(system.file("pangen/synteny_func_gap.R", package = "pannagram"))
# source("visualisation/visualisation.R")

# pokazStage('Step 6. Alignment-2. Fill the gaps between synteny blocks')

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("--path.chr"),     type = "character", default = NULL, help = "Path to query chromosome fasta files"),
  make_option(c("--path.aln"),     type = "character", default = NULL, help = "Path to the output directory with alignments"),
  make_option(c("--path.gaps"),    type = "character", default = NULL, help = "Path to the directory with blast of gaps"),
  
  make_option(c("--ref"),          type = "character", default = NULL, help = "Name of the reference genome"),
  make_option(c("--accessions"),   type = "character", default = NULL, help = "File containing accessions to analyze"),
  make_option(c("--combinations"), type = "character", default = NULL, help = "File containing combinations to analyze"),
  
  make_option(c("--cores"),        type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option(c("--path.log"),     type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"),    type = "character", default = NULL, help = "Level of log to be shown on the screen")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

pos.beg.info = 2  # Position in sequence's name, where the genome's begin is started

# print(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- Values of parameters ----

# Number of cores
num.cores <- opt$cores

# Reference genome
if (!is.null(opt$ref)) base.acc <- opt$ref

# Paths
if (!is.null(opt$path.chr)) path.chr <- opt$path.chr
if (!is.null(opt$path.aln)) path.aln <- opt$path.aln
if (!is.null(opt$path.gaps)) path.gaps <- opt$path.gaps

if(!dir.exists(path.aln)) dir.create(path.aln)
if(!dir.exists(path.gaps)) dir.create(path.gaps)

# Accessions
file.acc <- ifelse(!is.null(opt$accessions), opt$accessions, stop("File with accessions are not specified"))
tmp <- read.table(file.acc, stringsAsFactors = F)
accessions <- as.character(tmp[,1])
pokaz('Names of genomes for the analysis:', accessions, 
      file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- Combinations ----

file.combinations <- ifelse(!is.null(opt$combinations), opt$combinations, stop("File with combinations are not specified"))
files.maj <- list.files(path.aln, pattern = "\\maj.rds$")

# Filter files that start with one of the values in accessions
files.maj <- files.maj[sapply(files.maj, function(x) any(sapply(accessions, function(a) startsWith(x, a))))]

# Additional filtration when combinations are provided
if (length(readLines(file.combinations)) != 0) {
  combinations = read.table(file.combinations)
  files.maj.comb = c()
  for(acc in accessions){
    for(i.comb in 1:nrow(combinations)){
      file.tmp = paste0(acc, '_', 
                        combinations[i.comb,1], '_',
                        combinations[i.comb,2], '_maj.rds')
      if(file.exists(paste0(path.aln, file.tmp))){
        files.maj.comb = c(files.maj.comb, file.tmp)  
      }
    }
  }
  files.maj = intersect(files.maj, files.maj.comb)
}

if(length(files.maj) == 0) stop('No alignments provided')
pokaz('Number of alignments:', length(files.maj), file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

check.genomes = F

loop.function <- function(f.maj,
                          echo.loop=T){
  
  initial.vars <- ls()
  
  # Remove extensions
  pref.comb <- sub("\\_maj.rds$", "", f.maj)
  
  # Log files
  file.log.loop = paste0(path.log, 'loop_file_', 
                         pref.comb,
                         '.log')
  if(!file.exists(file.log.loop)){
    invisible(file.create(file.log.loop))
  }
  
  # ---- Check log Done ----
  if(checkDone(file.log.loop)){
    return()
  }
  
  pokaz(file.log.loop)
  
  # --- --- --- --- --- --- --- --- --- --- ---
  # Parser prefix of the result file
  info <- pref2info(pref.comb)
  query.chr <- info$query.chr
  base.chr <- info$base.chr
  acc <- info$acc
  
  # Files: in and out
  file.aln.pre <- paste(path.aln, paste0(pref.comb, '_maj.rds', collapse = ''), sep = '')
  file.aln.full <- paste(path.aln, paste0(pref.comb,  '_full.rds', collapse = ''), sep = '')
  
  # If the previous file was not created - next
  if(!file.exists(file.aln.pre)) {
    pokaz('No file', file.aln.pre, file=file.log.loop, echo=echo.loop)
    pokaz('Done.', file=file.log.loop, echo=echo.loop)
    
    final.vars <- ls()
    new.vars <- setdiff(final.vars, initial.vars)
    rm(list = new.vars)
    gc()
    
    return(NULL)
  }
  
  pokaz('Alignment:', acc, query.chr, base.chr, file=file.log.loop, echo=echo.loop)
  
  # ---- Read genomes ----

  if(check.genomes){
    # Read reference sequences
    base.file = paste0(base.acc, '_chr', base.chr , '.fasta', collapse = '')
    pokaz('Base:', base.file, file=file.log.loop, echo=echo.loop)
    base.fas.fw = readFastaMy(paste0(path.chr, base.file))
    base.fas.fw = seq2nt(base.fas.fw)
    base.fas.bw = revCompl(base.fas.fw)
    base.len = length(base.fas.bw)
    
    # Read query sequences
    query.file = paste0(acc, '_chr',query.chr, '.fasta')
    pokaz('Query:', query.file, file=file.log.loop, echo=echo.loop)
    
    query.fas.chr = readFastaMy(paste0(path.chr, query.file))
    query.fas.chr = seq2nt(query.fas.chr)
    query.len = length(query.fas.chr)
  }
  
  # ---- Read Major ----
  
  pokaz('Read the skeleton alignment..', file=file.log.loop, echo=echo.loop)
  x.sk = readRDS(file.aln.pre)
  x.sk = cleanOverlaps(x.sk)
  
  if(check.genomes){
    x.dir = setDir(x.sk, base.len = base.len)
    checkCorrespToGenome(x.dir, query.fas = query.fas.chr,
                         base.fas.fw = base.fas.fw,
                         base.fas.bw = base.fas.bw)
  }
  
  pokaz('after skeleton')
  
  # Length of the base chromosome
  # OLD APPROXIMATION: max.chr.len = max(max(x.sk$V4), max(x.sk$V5)) + 10^6  # length of the reference genome
  
  file.acc.len = paste0(path.chr, base.acc, '_chr_len.txt', collapse = '')
  if(!file.exists(file.acc.len)){
    stop('Chromosole length file is nor provided!')
  }
  chr.len = read.table(file.acc.len, stringsAsFactors = F, header = 1)
  max.chr.len = chr.len$len[as.numeric(base.chr)]
  pokaz('Length of the reference chromosome is ', max.chr.len)
  
  # ---- Read Gaps between normal blocks ----
  
  complexity.threshold = 200  # Max number of blast hits between two synteny blocks
  
  file.gaps.out = paste0(path.gaps,
                         'acc_', acc, 
                         '_qchr_', query.chr, '_bchr_', base.chr, '_out.txt', collapse = '')
  
  pokaz('gap file', file.gaps.out, file=file.log.loop, echo=echo.loop)

  if(file.exists(file.gaps.out)){
    pokaz('Read blast of good gaps..', file=file.log.loop, echo=echo.loop)
    x.gap = readBlast(file.gaps.out)
    
    # save(list = ls(), file = "tmp_workspace_good.RData")
    
  } else {
    x.gap = NULL
  }
  
  # To catch possible bugs
  # if(file.gaps.out == 'your_filename.txt'){
  #   file.ws = "tmp_workspace.RData"
  #   
  #   all.local.objects <- ls()
  #   save(list = all.local.objects, file = file.ws)
  #   
  #   pokaz('Workspace is saved in', file.ws, file=file.log.loop, echo=echo.loop)
  #   stop('Enough..')
  # }
  
  
  if(!is.null(x.gap)){
    pokaz('Number of gaps', nrow(x.gap), file=file.log.loop, echo=echo.loop)
    
    x.gap$pref1 = sapply(x.gap$V1, function(s) strsplit(s, '_query')[[1]][1])
    x.gap$pref2 = sapply(x.gap$V10, function(s) strsplit(s, '_base')[[1]][1])
    x.gap = x.gap[x.gap$pref1 == x.gap$pref2,]
    x.gap = unique(x.gap)  # UNIQUE
    pokaz('nrow', nrow(x.gap), file=file.log.loop, echo=echo.loop)
    if(nrow(x.gap) == 0) x.gap = NULL
  }  # KOSTYL
  
  if(!is.null(x.gap)){  # KOSTYL
    
    x.gap$q.beg = as.numeric(sapply(x.gap$V1, function(s) strsplit(s, '\\|')[[1]][pos.beg.info])) - 1
    x.gap$b.beg = as.numeric(sapply(x.gap$V10, function(s) strsplit(s, '\\|')[[1]][pos.beg.info])) - 1
    
    # x.gap$q.end = as.numeric(sapply(x.gap$V1, function(s) strsplit(s, '\\|')[[1]][pos.beg.info+1]))
    # x.gap$b.end = as.numeric(sapply(x.gap$V10, function(s) strsplit(s, '\\|')[[1]][pos.beg.info+1]))
    
    x.gap$V2 = x.gap$V2 + x.gap$q.beg
    x.gap$V3 = x.gap$V3 + x.gap$q.beg
    x.gap$V4 = x.gap$V4 + x.gap$b.beg
    x.gap$V5 = x.gap$V5 + x.gap$b.beg
    x.gap$dir = (x.gap$V4 > x.gap$V5) * 1
    
    # if(sum(x.gap$V3 > x.gap$q.end) > 0) stop('query')
    # if(sum(x.gap$V5 > x.gap$b.end) > 0) stop('base')
    
    pokaz('before1')
    if(check.genomes){
      checkCorrespToGenome(setDir(x.gap, base.len = base.len), 
                           query.fas = query.fas.chr,
                           base.fas.fw = base.fas.fw,
                           base.fas.bw = base.fas.bw)
      pokaz('after1')
    }
    
    # save(list = ls(), file = "tmp_workspace.RData")
    
    x.gap = glueZero(x.gap)
    x.gap$idx = 1:nrow(x.gap)
    
    pokaz('before11')
    if(check.genomes){
      checkCorrespToGenome(setDir(x.gap, base.len = base.len), 
                           query.fas = query.fas.chr,
                           base.fas.fw = base.fas.fw,
                           base.fas.bw = base.fas.bw)
      pokaz('after11')
    }
    
    cnt = table(x.gap$pref1)
    
    idx.good = c()
    for(i in 1:length(cnt)){
      if(i %% 100 == 0) pokaz('Pgress: Number of analysed gaps', i, file=file.log.loop, echo=echo.loop)
      
      # name of the node
      s = names(cnt)[i]
      x.tmp = x.gap[x.gap$pref1 == s,]
      
      # If only one BLAST-hit - get it as it is.
      if(nrow(x.tmp) == 1){
        idx.good = c(idx.good, x.tmp$idx)
        next
      }

      # Clean overlaps
      x.tmp = cleanOverlaps(x.tmp)
      if(nrow(x.tmp) == 1){
        idx.good = c(idx.good, x.tmp$idx)
        next
      } 
      
      # Transform positions to positive
      df.gap = transform_positions(x.tmp)
      
      ## ---- Greedy loop ----
      
      overlap.cutoff = 0.2
      idx.added = greedy_loop(df.gap, overlap.cutoff)
      
      # ---- Remaining ----
      idx.good = c(idx.good, df.gap$idx[idx.added])
      
    }
    
    if(length(idx.good) > 0){
      x.res = x.gap[idx.good,]  
      # Clean overlaps from both (base and query) sides
      x.res = cleanOverlaps(x.res)
      
      pokaz('before13')
      if(check.genomes){
        checkCorrespToGenome(setDir(x.tmp, base.len = base.len), 
                             query.fas = query.fas.chr,
                             base.fas.fw = base.fas.fw,
                             base.fas.bw = base.fas.bw)
        pokaz('after13')
      }
      
    } else {
      # Create empty
      x.res <- data.frame(matrix(NA, nrow = 0, ncol = length(colnames(x.sk)), 
                                 dimnames = list(NULL, colnames(x.sk))))
    }
    rmSafe(x.gap)
    
  } else {
    # Create empty
    x.res <- data.frame(matrix(NA, nrow = 0, ncol = length(colnames(x.sk)), 
                               dimnames = list(NULL, colnames(x.sk))))
  }
  
  
  pokaz('before2')
  if(nrow(x.res) > 0){
    if(check.genomes){
      checkCorrespToGenome(setDir(x.res, base.len = base.len), 
                           query.fas = query.fas.chr,
                           base.fas.fw = base.fas.fw,
                           base.fas.bw = base.fas.bw)  
      pokaz('after2')
    }
  }
  
  
  
  # ---- Read additional alignments ----
  
  file.gaps.out = paste0(path.gaps,
                         'acc_', acc, 
                         '_qchr_', query.chr, '_bchr_', base.chr, '_residual_out.txt', collapse = '')
  x.gap = NULL
  if(file.exists(file.gaps.out)){
    pokaz('Read blast of "bad" gaps..', file.gaps.out, file=file.log.loop, echo=echo.loop)
    x.gap = readBlast(file.gaps.out)
    x.gap = unique(x.gap)
  }

  # save(list = ls(), file = "tmp_workspace1.RData")

  if(!is.null(x.gap)) {
    
    ## ---- Change positions ----
    
    pos.shift = posShift(x.gap)
    
    x.gap[,2:3] = x.gap[,2:3] + pos.shift$q[x.gap$V1,]$shift
    x.gap[,4:5] = x.gap[,4:5] + pos.shift$r[x.gap$V10,]$shift
    x.gap$dir = (x.gap$V4 > x.gap$V5) * 1
    
    # pokaz('Saving..')
    # save(list = ls(), file = "tmp_workspace.RData")
    
    x.tmp = glueZero(x.gap)
    x.tmp$idx = 1:nrow(x.tmp)
    # plotSynDot(x.tmp)
    
    # ---- Remain only those, that have the intersection in numbers ----
    getConnectingBlocks <- function(s){
      s = strsplit(s, '_connect_')[[1]][2] 
      s = strsplit(s, '_')[[1]][1:2]
      return(s)
    }
    num.q = sapply(x.tmp$V1, getConnectingBlocks)
    num.r = sapply(x.tmp$V10, getConnectingBlocks)
    
    idx.remain = ((num.q[1,] == num.r[1,]) | 
              (num.q[1,] == num.r[2,]) | 
              (num.q[2,] == num.r[1,]) | 
              (num.q[2,] == num.r[2,]))
    
    x.tmp = x.tmp[idx.remain,]
    rownames(x.tmp) = NULL
    
    # plotSynDot(x.tmp)
  } else {
    x.tmp = data.frame(matrix(NA, nrow = 0, ncol = length(colnames(x.sk)), 
                            dimnames = list(NULL, colnames(x.sk))))
  }
    
  if(nrow(x.tmp) > 0){
    # Clean the overlap
    
    x.tmp = cleanOverlaps(x.tmp)
   
    id.corresp = c()
    for(irow in 1:nrow(x.tmp)){
      if(x.tmp$dir[irow] == 0){
        tmp = which((x.tmp$V2[irow] <= x.gap$V3) & (x.gap$V2 <= x.tmp$V3[irow]) & 
                      (x.tmp$V4[irow] <= x.gap$V5) & (x.gap$V4 <= x.tmp$V5[irow]) & (x.gap$dir == 0))  
      } else {
        tmp = which((x.tmp$V2[irow] <= x.gap$V3) & (x.gap$V2 <= x.tmp$V3[irow]) & 
                      (x.tmp$V5[irow] <= x.gap$V4) & (x.gap$V5 <= x.tmp$V4[irow]) & (x.gap$dir != 0))  
      }
      if(length(tmp) == 0) stop('Wrong, no correspondence')
      id.corresp = rbind(id.corresp, cbind(tmp, irow))
    }
    id.corresp <- setNames(as.data.frame(id.corresp), c('init', 'new'))
    
    
    # IDX
    if(nrow(x.tmp) == 1){
      idx.bw = unique(id.corresp$init[id.corresp$new == 1])  # UNIQUE
    } else {
      
      # ---- The same greedy as before ----
      # Transform positions to positive
      
      df.gap = transform_positions(x.tmp)
      
      ## ---- Greedy loop ----
      
      overlap.cutoff = 0.2
      idx.remain = greedy_loop(df.gap, overlap.cutoff)
      
      idx.bw = unique(id.corresp$init[id.corresp$new %in% idx.remain])  # UNIQUE
    }
    
    # ---- Remaining ----
    
    idx.bw = unique(idx.bw)   # UNIQUE
    x.bw = x.gap[idx.bw,]
    
    # Change positions back
    x.bw[,2:3] = x.bw[,2:3] - pos.shift$q[x.bw$V1,]$shift
    x.bw[,4:5] = x.bw[,4:5] - pos.shift$r[x.bw$V10,]$shift
    
    # Shift positions to the initial
    x.bw$q.beg = as.numeric(sapply(x.bw$V1, function(s) strsplit(s, '\\|')[[1]][pos.beg.info])) - 1
    x.bw$b.beg = as.numeric(sapply(x.bw$V10, function(s) strsplit(s, '\\|')[[1]][pos.beg.info])) - 1
    
    # x.bw$q.end = as.numeric(sapply(x.bw$V1, function(s) strsplit(s, '\\|')[[1]][pos.beg.info+1]))
    # x.bw$b.end = as.numeric(sapply(x.bw$V10, function(s) strsplit(s, '\\|')[[1]][pos.beg.info+1]))
    
    x.bw$V2 = x.bw$V2 + x.bw$q.beg
    x.bw$V3 = x.bw$V3 + x.bw$q.beg
    x.bw$V4 = x.bw$V4 + x.bw$b.beg
    x.bw$V5 = x.bw$V5 + x.bw$b.beg
    x.bw$dir = (x.bw$V4 > x.bw$V5) * 1
    
    x.bw = glueZero(x.bw)
    x.bw = cleanOverlaps(x.bw)
    
    pokaz('before3')
    if(check.genomes){
      checkCorrespToGenome(setDir(x.bw, base.len = base.len),
                           query.fas = query.fas.chr,
                           base.fas.fw = base.fas.fw,
                           base.fas.bw = base.fas.bw)
      pokaz('after3')
    }
    
  } else {
    
    # Create empty
    x.bw <- data.frame(matrix(NA, nrow = 0, ncol = length(colnames(x.sk)), 
                              dimnames = list(NULL, colnames(x.sk))))
  }
  
  # ---- Combine all together ----
  
  comb.names = intersect(intersect(colnames(x.sk), colnames(x.bw)), colnames(x.res))
  x.comb = rbind(x.res[,comb.names], x.sk[,comb.names])
  x.comb = rbind(x.comb, x.bw[,comb.names])
  x.comb = x.comb[order(x.comb$V2),]
  rownames(x.comb) = NULL
  
  x.comb = glueZero(x.comb)
  
  # # To catch possible bugs
  # if(T){
  #   file.ws = "tmp_workspace.RData"
  # 
  #   all.local.objects <- ls()
  #   save(list = all.local.objects, file = file.ws)
  # 
  #   pokaz('Workspace is saved in', file.ws, file=file.log.loop, echo=echo.loop)
  #   stop('Enough..')
  # }
  
  # ---- Check uniqueness ---- 
  
  x.sk1 = x.comb
  # x.sk1 = x.res
  # x.sk1 = x.bw
  
  # save(list = ls(), file = "tmp_workspace.RData")
  
  rownames(x.sk1) = NULL
  pos.q.occup = rep(0, max.chr.len)
  for(irow in 1:nrow(x.sk1)){
    pp = x.sk1$V4[irow]:x.sk1$V5[irow]
    if(sum(pos.q.occup[pp]) > 0) {
      pokaz('Non-unique', file=file.log.loop, echo=echo.loop)
      stop('non-unique') 
    }
    # pos.q.occup[pp] = pos.q.occup[pp] + 1
    pos.q.occup[pp] = irow
  }
  sum(pos.q.occup > 1)
  sum(pos.q.occup)
  
  # ---- Check genomes ----
  
  # save(list = ls(), file = "tmp_workspace.RData")
  # stop('Enough..')
  
  if(check.genomes){
    x.dir = setDir(x.comb, base.len = base.len)
    checkCorrespToGenome(x.dir, query.fas = query.fas.chr,
                         base.fas.fw = base.fas.fw,
                         base.fas.bw = base.fas.bw)
  }
  
  saveRDS(object = x.comb, file = file.aln.full)
  

  # Done
  pokaz('Done.', file=file.log.loop, echo=echo.loop)
  
  # Cleanup variables
  final.vars <- ls()
  new.vars <- setdiff(final.vars, initial.vars)
  rm(list = new.vars)
  gc()
  return(NULL)
}
  

# ***********************************************************************
# ---- Loop  ----


if(num.cores == 1){
  for(f.maj in files.maj){
    loop.function(f.maj, echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  
  batch.size <- 2 * num.cores
  
  # Loop through files.maj in batches
  tmp <- list()
  for (start in seq(1, length(files.maj), by = batch.size)) {
    
    end <- min(start + batch.size - 1, length(files.maj))  # Define the end of the current batch
    batch.files <- files.maj[start:end]
    
    # Create and register a new cluster for the current batch
    myCluster <- makeCluster(num.cores, type = "PSOCK")
    registerDoParallel(myCluster)
    
    # Run parallel loop for files in the current batch
    batch.results <- foreach(f.maj = batch.files, .packages = c('crayon')) %dopar% {
      loop.function(f.maj, echo.loop = echo.loop)
    }
    tmp <- c(tmp, batch.results)
    
    stopCluster(myCluster)  # Stop the cluster after completing the batch
    
    # pokaz('Memory', mem_used())
    # gc()
    # pokaz('Memory after gc', mem_used())
  }
}

warnings()

pokaz('Done.', file=file.log.main, echo=echo.main)

# ---- Manual testing  ----

# Bash:
# for file in *; do
#  if [ -f "$file" ]; then
#    last_line=$(tail -n 1 "$file")
#    if [[ "$last_line" != *Done* ]]; then
#      echo "$file: $last_line"
#    fi
#  fi
# done


# To catch possible bugs
# if(file.gaps.out == 'your_filename.txt'){
#   file.ws = "tmp_workspace.RData"
#   
#   all.local.objects <- ls()
#   save(list = all.local.objects, file = file.ws)
#   
#   pokaz('Workspace is saved in', file.ws, file=file.log.loop, echo=echo.loop)
#   stop('Enough..')
# }


# file.ws = "tmp_workspace.RData"
# all.local.objects <- ls()
# save(list = all.local.objects, file = file.ws)
# pokaz('Workspace is saved in', file.ws, file=file.log.loop, echo=echo.loop)
# stop('Enough..')



# save(list = ls(), file = "tmp_workspace.RData")
# stop('Enough..')

# 
# checkCorrespToGenome(setDir(x.comb, base.len = base.len), query.fas = query.fas.chr,
#                      base.fas.fw = base.fas.fw,
#                      base.fas.bw = base.fas.bw)
