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

# files.maj = c(
#   "100326_1_1_maj.rds",
#   "100324_1_1_maj.rds",
#   "100325_1_1_maj.rds",
#   "100323_3_3_maj.rds",
#   "100324_3_3_maj.rds",
#   "100323_4_4_maj.rds",
#   "100326_4_4_maj.rds",
#   "100325_3_3_maj.rds")
# 
# pokaz(files.maj)

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
                          done.set = character(0),
                          echo.loop=T){

  initial.vars <- ls()

  # Remove extensions
  pref.comb <- sub("\\_maj.rds$", "", f.maj)

  # ---- Checkpoint: skip already completed items ----
  item.id <- pref.comb
  if(item.id %in% done.set){
    return(NULL)
  }

  # One log file per worker (bounded number of files); also the checkpoint ledger
  file.log.loop = initLoopLog(path.log)

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
    markDone(item.id, file=file.log.loop, echo=echo.loop)

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
  
  # Normal gaps are now blasted in BATCHES (synteny_03/04): read & rbind all
  # <pref>b<N>_out.txt for this comparison (excludes *_residual_out.txt).
  gap.pref = paste0('acc_', acc, '_qchr_', query.chr, '_bchr_', base.chr, '_')
  all.out.files = list.files(path.gaps, full.names = TRUE)
  batch.out.files = all.out.files[startsWith(basename(all.out.files), gap.pref) &
                                  grepl('_b[0-9]+_out\\.txt$', basename(all.out.files))]

  pokaz('gap batch files', length(batch.out.files), file=file.log.loop, echo=echo.loop)

  if(length(batch.out.files) > 0){
    pokaz('Read blast of good gaps (batched)..', file=file.log.loop, echo=echo.loop)
    x.gap = do.call(rbind, lapply(batch.out.files, readBlast))
    if(!is.null(x.gap) && nrow(x.gap) == 0) x.gap = NULL
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
    
    x.gap$pref1 = sub('_query.*', '', x.gap$V1)   # vectorised (was sapply+strsplit per row)
    x.gap$pref2 = sub('_base.*', '', x.gap$V10)
    x.gap = x.gap[x.gap$pref1 == x.gap$pref2,]
    x.gap = unique(x.gap)  # UNIQUE
    pokaz('nrow', nrow(x.gap), file=file.log.loop, echo=echo.loop)
    if(nrow(x.gap) == 0) x.gap = NULL
  }  # KOSTYL
  
  if(!is.null(x.gap)){  # KOSTYL
    
    x.gap$q.beg = as.numeric(data.table::tstrsplit(x.gap$V1, '|', fixed = TRUE)[[pos.beg.info]]) - 1   # vectorised split
    x.gap$b.beg = as.numeric(data.table::tstrsplit(x.gap$V10, '|', fixed = TRUE)[[pos.beg.info]]) - 1
    
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
    
    # Split gaps by their connecting-block id ONCE (was x.gap[x.gap$pref1 == s,] per
    # group -> O(n_groups * n_hits)). split() groups by sorted factor levels == names(table),
    # preserving within-group row order, so the result is identical.
    gap.groups = split(x.gap, x.gap$pref1)

    idx.good = c()
    for(i in seq_along(gap.groups)){
      if(i %% 100 == 0) pokaz('Pgress: Number of analysed gaps', i, file=file.log.loop, echo=echo.loop)

      x.tmp = gap.groups[[i]]
      
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
    rm(x.gap, x.tmp, df.gap, cnt, idx.good)
    
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
    # Vectorised extraction of the two connecting-block numbers from '..._connect_<a>_<b>_...'
    # (was sapply + two strsplit per row). Rows: [1,]=<a>, [2,]=<b> (character, as before).
    cre = '.*_connect_([0-9]+)_([0-9]+)_.*'
    num.q = rbind(sub(cre, '\\1', x.tmp$V1),  sub(cre, '\\2', x.tmp$V1))
    num.r = rbind(sub(cre, '\\1', x.tmp$V10), sub(cre, '\\2', x.tmp$V10))
    
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
   
    # Vectorised interval-containment join (replaces an O(n^2) per-block `which` loop):
    # a gap corresponds to a block when QUERY ranges overlap AND BASE ranges overlap AND
    # direction matches. Base ranges use min/max so the overlap test is direction-agnostic.
    dt.b = data.table::data.table(new = seq_len(nrow(x.tmp)),
                                  q1 = x.tmp$V2, q2 = x.tmp$V3,
                                  b1 = pmin(x.tmp$V4, x.tmp$V5), b2 = pmax(x.tmp$V4, x.tmp$V5),
                                  bd = x.tmp$dir)
    dt.g = data.table::data.table(init = seq_len(nrow(x.gap)),
                                  q1 = x.gap$V2, q2 = x.gap$V3,
                                  b1 = pmin(x.gap$V4, x.gap$V5), b2 = pmax(x.gap$V4, x.gap$V5),
                                  bd = x.gap$dir)
    data.table::setkey(dt.b, q1, q2)
    ov = data.table::foverlaps(dt.g, dt.b, by.x = c('q1','q2'), by.y = c('q1','q2'),
                               type = 'any', nomatch = NULL)
    keep = (ov$b1 <= ov$i.b2) & (ov$i.b1 <= ov$b2) & (ov$bd == ov$i.bd)  # base overlap + same dir
    id.corresp = data.frame(init = ov$init[keep], new = ov$new[keep])
    # Restore the original loop's ordering (by block, then by gap index) so downstream
    # order-sensitive overlap trimming is bit-identical to the pre-optimisation code.
    id.corresp = id.corresp[order(id.corresp$new, id.corresp$init), ]
    if(length(unique(id.corresp$new)) != nrow(x.tmp)) stop('Wrong, no correspondence')
    
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
    x.bw$q.beg = as.numeric(data.table::tstrsplit(x.bw$V1, '|', fixed = TRUE)[[pos.beg.info]]) - 1   # vectorised split
    x.bw$b.beg = as.numeric(data.table::tstrsplit(x.bw$V10, '|', fixed = TRUE)[[pos.beg.info]]) - 1
    
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
    
    rm(x.gap, x.tmp, pos.shift, id.corresp, num.q, num.r, df.gap, idx.bw, idx.remain)
    gc()
    
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
  
  rm(x.sk, x.res, x.bw)
  gc()
  
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
  
  # x.sk1 = x.comb
  
  # x.sk1 = x.res
  # x.sk1 = x.bw
  
  # save(list = ls(), file = "tmp_workspace.RData")
  
  # rownames(x.sk1) = NULL
  # pos.q.occup = rep(0, max.chr.len)
  # for(irow in 1:nrow(x.sk1)){
  #   pp = x.sk1$V4[irow]:x.sk1$V5[irow]
  #   if(sum(pos.q.occup[pp]) > 0) {
  #     pokaz('Non-unique', file=file.log.loop, echo=echo.loop)
  #     stop('non-unique') 
  #   }
  #   # pos.q.occup[pp] = pos.q.occup[pp] + 1
  #   pos.q.occup[pp] = irow
  # }
  # sum(pos.q.occup > 1)
  # sum(pos.q.occup)
  # rm(pos.q.occup)
  # gc()
  
  # ---- Check genomes ----
  
  # save(list = ls(), file = "tmp_workspace.RData")
  # stop('Enough..')
  
  if(check.genomes){
    x.dir = setDir(x.comb, base.len = base.len)
    checkCorrespToGenome(x.dir, query.fas = query.fas.chr,
                         base.fas.fw = base.fas.fw,
                         base.fas.bw = base.fas.bw)
  }
  
  saveRDS(object = x.comb, file = file.aln.full, compress = FALSE)  # gzip default dominated I/O; readRDS reads either
  

  # ---- Checkpoint marker: item fully processed ----
  markDone(item.id, file=file.log.loop, echo=echo.loop)

  # Cleanup variables
  final.vars <- ls()
  new.vars <- setdiff(final.vars, initial.vars)
  rm(list = new.vars)
  gc()
  return(NULL)
}
  

# ***********************************************************************
# ---- Loop  ----


# Rebuild the set of completed items from all worker logs (any core count)
done.set <- getDoneSet(path.log)
if(length(done.set) > 0){
  pokaz('Skip already done:', length(done.set), file=file.log.main, echo=echo.main)
}

if(num.cores == 1){
  assign('.worker.id', 1, envir = .GlobalEnv)   # stable single file core_1.log
  for(f.maj in files.maj){
    loop.function(f.maj, done.set = done.set, echo.loop=echo.loop)
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

    # Assign a stable worker id (1..N) to each worker -> bounded, reused file names
    parallel::clusterApply(myCluster, seq_len(num.cores),
                           function(i) assign('.worker.id', i, envir = .GlobalEnv))

    # Run parallel loop for files in the current batch
    invisible(
      foreach(f.maj = batch.files, .packages = c('crayon')) %dopar% {
        loop.function(f.maj, done.set = done.set, echo.loop = echo.loop)
        NULL
      }
    )
    
    stopCluster(myCluster)  # Stop the cluster after completing the batch
    
    # pokaz('Memory', mem_used())
    # gc()
    # pokaz('Memory after gc', mem_used())
  }
}


# stop('all files tested')

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
