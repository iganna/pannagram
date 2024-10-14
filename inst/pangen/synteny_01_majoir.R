# Alignment-1. Remaining syntenic (major) matches

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/synteny_func.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("--path.blast"),    type = "character", default = NULL, help = "Path to blast results"),
  make_option(c("--path.aln"),      type = "character", default = NULL, help = "Path to the output directory with alignments"),
  
  make_option(c("--ref"),           type = "character", default = NULL, help = "Name of the reference genome"),
  make_option(c("--accessions"),    type = "character", default = NULL, help = "File containing accessions to analyze"),
  make_option(c("--combinations"),  type = "character", default = NULL, help = "File containing combinations to analyze"),
  
  make_option(c("--cores"),         type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option(c("--path.log"),      type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"),     type = "character", default = NULL, help = "Level of log to be shown on the screen")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Variables ----

max.len = 10^6
len.blast = 50

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- Values of parameters ----

# Number of cores
num.cores <- ifelse(!is.null(opt$cores), opt$cores, 30)

path.blast    <- ifelse(!is.null(opt$path.blast), opt$path.blast, stop('Folder with BLAST results is not specified'))
path.aln      <- ifelse(!is.null(opt$path.aln), opt$path.aln, stop('Folder with Alignments is not specified'))
base.acc      <- ifelse(!is.null(opt$ref), opt$ref, stop('Reference genome is not specified'))

# Create folders for the alignment results
if(!dir.exists(path.aln)) dir.create(path.aln)

# Accessions
file.acc <- ifelse(!is.null(opt$accessions), opt$accessions, stop("File with accessions are not specified"))
tmp <- read.table(file.acc, stringsAsFactors = F)
accessions <- tmp[,1]
pokaz('Names of genomes for the analysis:', accessions, 
      file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- Combinations ----

file.combinations <- ifelse(!is.null(opt$combinations), opt$combinations, stop("File with combinations are not specified"))
if (length(readLines(file.combinations)) == 0) {
  files.blast <- list.files(path.blast, pattern = "\\.txt$", full.names = F)

  # Filter files that start with one of the values in accessions
  files.blast <- files.blast[sapply(files.blast, function(x) any(sapply(accessions, function(a) startsWith(x, a))))]
  
  pokaz('All blast files', length(files.blast), file=file.log.main, echo=echo.main)
} else {
  combinations = read.table(file.combinations)
  files.blast = c()
  for(acc in accessions){
    for(i.comb in 1:nrow(combinations)){
      file.tmp = paste0(acc, '_', 
                        combinations[i.comb,1], '_',
                        combinations[i.comb,2], '.txt')
      if(!file.exists(paste0(path.blast, file.tmp))) stop(paste0('NO FILE', file.tmp))
      
      files.blast = c(files.blast, file.tmp)
    }
  }
}

if(length(files.blast) == 0) stop('No BLAST files provided')
pokaz('Number of BLAST-result files:', length(files.blast), file=file.log.main, echo=echo.main)


# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(f.blast, 
                          file.log.loop=NULL,
                          echo.loop=T){
  
  # Output file
  pref.comb <- sub("\\.[^.]*$", "", basename(f.blast))
  file.aln.pre <- paste(path.aln, paste0(pref.comb, '_maj.rds', collapse = ''), sep = '')
  
  # Log files
  file.log.loop = paste0(path.log, 'loop_file_', 
                         pref.comb, # remove the extensions
                         '.log')
  if(!file.exists(file.log.loop)){
    invisible(file.create(file.log.loop))
  }
  
  # ---- Check log Done ----
  if(checkDone(file.log.loop)){
    return(NULL)
  }

  # --- --- --- --- --- --- --- --- --- --- ---
  
  # Parser for BLAST-result file
  parts <- strsplit(pref.comb, "_")[[1]]
  
  base.chr <- parts[length(parts)]
  query.chr <- parts[length(parts) - 1]
  parts = parts[-c(length(parts) - 1, length(parts))]
  acc <- paste0(parts, collapse = '_')
  
  pokaz('Alignment:', acc, query.chr, base.chr, file=file.log.loop, echo=echo.loop)
  
  # ---- Read blast results ----
  # pokaz(paste0(path.blast, f.blast))
  x = readBlast(paste0(path.blast, f.blast))
  if(is.null(x)){
    pokaz('Done.', file=file.log.loop, echo=echo.loop)
    return(NULL)
  }
  
  pokaz('Read blast results finished, number of rows is', nrow(x), file=file.log.loop, echo=echo.loop)
  
  ## ---- Pre-processing ----
  # Save true base coordinate
  x$p.beg <- ifelse(x$V4 < x$V5, x$V4, x$V5)
  x$p.end <- ifelse(x$V4 < x$V5, x$V5, x$V4)
  
  # Set correct position
  start.pos = as.numeric(sapply(strsplit(x[,1], "\\|"), "[", 4)) - 1
  x[,2:3] = x[,2:3] + start.pos
  
  # Set direction
  x$dir = (x$V4 > x$V5) * 1
  # x = glueZero(x)  # DO NOT DO IT
  
  # Set past id
  x$part.id = cumsum(c(T, diff(as.numeric(as.factor(x$V1))) != 0))
  
  # ---- Major skeleton ----
  
  idx.maj = which(!c(F, x$V1[-1] == x$V1[-nrow(x)]))
  x.major = cbind(x[idx.maj,setdiff(colnames(x), c('V8', 'V9'))], idx.maj)
  
  pokaz('Number of rows in the synteny', nrow(x.major), file=file.log.loop, echo=echo.loop)
  if(nrow(x.major) == 0) {
    pokaz('Done.', file=file.log.loop, echo=echo.loop)
    return(NULL)
  }
  
  # Sort to set the id
  x.major = x.major[order(x.major$V2),]  # not needed
  x.major$id = 1:nrow(x.major)
  
  # Sort by the position in the reference
  x.major = x.major[order(-x.major$V7),]
  x.major = x.major[order(x.major$p.beg),]
  
  # If complete overlap - remove shortest
  idx.overlap = which(x.major$p.end[-1] <= x.major$p.end[-nrow(x.major)]) + 1
  if(length(idx.overlap) > 0){
    x.major = x.major[-idx.overlap,]  
  }
  
  # Remove those, which are not in the correct place
  x.major$id = rank(x.major$id)
  idx = which(abs(diff(x.major$id)) > 1)
  idx = intersect(idx, idx + 1)
  
  if(length(idx) > 0){
    x.major = x.major[-idx,,drop=F]
    if(nrow(x.major) == 0){
      pokaz('Done.', file=file.log.loop, echo=echo.loop)
      return(NULL)
    }
  }
  
  # ----  Define blocks in the skeleton ----
  
  # re-arrange IDs
  x.major$id = rank(x.major$id)
  # [1 or 0] - [begin or end] of block
  x.major$block = c(1, abs(x.major$id[-1] - x.major$id[-nrow(x.major)]) != 1) 
  # block ID
  x.major$block.id = cumsum(x.major$block)
  
  # Analyse only beginnings of blocks
  x.block = x.major[x.major$block == 1,]
  x.block$id = rank(x.block$id)
  bl.len = abs(tapply(x.major$p.end, x.major$block.id, max) - tapply(x.major$p.beg, x.major$block.id, min))
  bl.len = bl.len[order(as.numeric(names(bl.len)))]
  x.block$len = bl.len
  
  # Remain blocks, only if they are in a correct place and long enough
  while(T){
    x.block$id = rank(x.block$id)
    idx = which(abs(diff(x.block$id)) > 1)
    idx = unique(c(idx, idx + 1))
    idx = idx[idx < nrow(x.block)]
    # idx = intersect(idx, idx +1)
    if(length(idx) == 0) break
    if(min(x.block$len[idx]) > 20000) break
    idx.remove = idx[x.block$len[idx] == min(x.block$len[idx])][1]
    pokaz('- Remove block', idx.remove, ';length:', x.block$len[idx.remove], file=file.log.loop, echo=echo.loop)
    x.block = x.block[-idx.remove,]
    # rownames(x.block) = NULL
  }
  if(is.unsorted(x.major$p.beg)) pokazAttention('1!!', file=file.log.loop, echo=echo.loop)
  
  # Remove wrong blocks completely
  remain.block = x.block$block.id
  if(length(remain.block) > 0){
    x.major = x.major[x.major$block.id %in% remain.block, , drop=F]
    if(nrow(x.major) == 0) {
      pokaz('Done.', file=file.log.loop, echo=echo.loop)
      return(NULL)
    }
  }
  
  # Check remained blocks
  x.major$id = rank(x.major$id)
  x.major$block = c(1, abs(x.major$id[-1] - x.major$id[-nrow(x.major)]) != 1) 
  x.major$block = x.major$block + c(1, abs(x.major$dir[-1] != x.major$dir[-nrow(x.major)])) 
  x.major$block = (x.major$block > 0) * 1
  
  x.major$block.id = cumsum(x.major$block)  # block ID
  
  if(length(unique((x.major$dir))) > 1){  # different dirertions exist
    cnt = table(x.major$block.id, x.major$dir)
    if(sum(cnt[,1] * cnt[,2]) != 0) stop('Blocks in x.major are wrongly defined')
  }
  
  # ---- Filtration ----
  x = x[x.major$idx.maj,]
  
  x$block.id = x.major$block.id
  # saveRDS(x, file.aln.pre, compress = F)
  
  # ---- Remove short overlaps: twice, because from "both sides" ----
  for(i.tmp in 1:2){
    x = cleanBigOverlaps(x)
    x = cutSmallOverlaps(x)
  }
  
  # ---- Check uniquness of occupancy ----
  # pos.q.occup = rep(0, base.len)
  # for(irow in 1:nrow(x)){
  #   # pos.q.occup[x$V2[irow]:x$V3[irow]] = pos.q.occup[x$V2[irow]:x$V3[irow]] + 1
  #   # if(sum(pos.q.occup[x$V4[irow]:x$V5[irow]]) > 0) stop()
  #   pos.q.occup[x$V4[irow]:x$V5[irow]] = pos.q.occup[x$V4[irow]:x$V5[irow]] + 1
  # }
  # if(sum(pos.q.occup > 1) > 0){
  #   stop('Overlaps in base are remained')
  # } 
  # pokaz('Occupancy of base', sum(pos.q.occup), file=file.log.loop, echo=echo.loop)
  
  # Save
  saveRDS(x, file.aln.pre, compress = F)
  
  rmSafe(x)
  rmSafe(x.major)
  rmSafe(base.fas.bw)
  rmSafe(base.fas.fw)
  rmSafe(query.fas.chr)
  gc()
  
  pokaz('Done.', file=file.log.loop, echo=echo.loop)
  return(NULL)
}

# ***********************************************************************
# ---- Loop  ----


if(num.cores == 1){

  for(f.blast in files.blast){
    loop.function(f.blast, echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(f.blast = files.blast, 
                .packages=c('crayon'), 
                .verbose = F)  %dopar% { 
                  loop.function(f.blast, echo.loop=echo.loop)
                }
  stopCluster(myCluster)
}

warnings()

pokaz('Done.', file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- Manual testing ----



