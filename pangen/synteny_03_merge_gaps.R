suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
})


source("utils/utils.R")
source("pangen/synteny_func.R")
# source("visualisation/visualisation.R")

# pokazStage('Step 6. Alignment-2. Fill the gaps between synteny blocks')

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("--path.chr"), type="character", default=NULL, 
              help="path to query chomosome fasta files", metavar="character"),  
  make_option(c("--path.aln"), type="character", default=NULL, 
              help="path to the output directory with alignments", metavar="character"),
  make_option(c("--ref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("--path.gaps"), type="character", default=NULL, 
              help="prefix of the directory with gaps", metavar="character"),
  make_option(c("--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("--path.log"), type = "character", default = NULL,
              help = "Path for log files", metavar = "character"),
  make_option(c("--log.level"), type = "character", default = NULL,
              help = "Level of log to be shown on the screen", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Logging ----

source('utils/chunk_logging.R') # a common code for all R logging

# ---- Values of parameters ----

# Number of cores
num.cores <- ifelse(!is.null(opt$cores), opt$cores, 30)


if (!is.null(opt$path.chr)) path.chr <- opt$path.chr
if (!is.null(opt$path.aln)) path.aln <- opt$path.aln
if (!is.null(opt$ref)) base.acc <- opt$ref
if (!is.null(opt$path.gaps)) path.gaps <- opt$path.gaps


if(!dir.exists(path.aln)) dir.create(path.aln)
if(!dir.exists(path.gaps)) dir.create(path.gaps)


# ***********************************************************************
# ---- Preparation ----

files.maj <- list.files(path.aln, pattern = "\\maj.rds$")
pokaz('Number of alignment files:', length(files.maj), file=file.log.main, echo=echo.main)
# if(length(files.maj) == 0) stop('No alignment files provided')



# ***********************************************************************
# ---- MAIN program body ----


loop.function <- function(f.maj,
                          echo.loop=T,
                          file.log.loop=NULL){
  
  # Remove extensions
  pref.comb <- sub("\\maj.rds$", "", f.maj)
  
  # Log files
  if (is.null(file.log.loop)){
    file.log.loop = paste0(path.log, 'loop_file_', 
                           pref.comb,
                           '.log')
    invisible(file.create(file.log.loop))
  }
  
  # Parser for BLAST-result file
  parts <- strsplit(pref.comb, "_")[[1]]
  
  base.chr <- parts[length(parts) - 1]
  query.chr <- parts[length(parts)]
  parts = parts[-c(length(parts) - 1, length(parts))]
  acc <- paste0(parts, collapse = '_')
  
  pokaz(acc, query.chr, base.chr, file=file.log.loop, echo=echo.loop)
  
  pref.comb = paste0(acc, '_', query.chr, '_', base.chr, collapse = '')
  
  # If the previous file was not created - next
  file.aln.pre <- paste(path.aln, paste0(pref.comb, '_maj.rds', collapse = ''), sep = '')
  if(!file.exists(file.aln.pre)) {
    pokaz('No file', file.aln.pre, file=file.log.loop, echo=echo.loop)
    pokaz('Done.', file=file.log.loop, echo=echo.loop)
    return(NULL)
  }
  
  # If the full file has been already created - next
  file.aln.full <- paste(path.aln, paste0(pref.comb,  '_full.rds', collapse = ''), sep = '')
  if(file.exists(file.aln.full)) {
    pokaz('File exist', file.aln.full, file=file.log.loop, echo=echo.loop)
    pokaz('Done.', file=file.log.loop, echo=echo.loop)
    return(NULL)
  }
  
  pokaz('Alignment:', acc, query.chr, base.chr, file=file.log.loop, echo=echo.loop)
  
  # # ---- Read genomes ----
  # 
  # # Read reference sequences
  # base.file = paste0(base.acc, '_chr', base.chr , '.fasta', collapse = '')
  # pokaz('Base:', base.file, file=file.log.loop, echo=echo.loop)
  # base.fas.fw = readFastaMy(paste(path.chr, base.file, sep = ''))
  # base.fas.fw = seq2nt(base.fas.fw)
  # base.fas.bw = revCompl(base.fas.fw)
  # base.len = length(base.fas.bw)
  # 
  # # Read query sequences
  # query.file = paste(acc, '_chr',query.chr, '.fasta', sep = '')
  # pokaz('Query:', query.file, file=file.log.loop, echo=echo.loop)
  # 
  # query.fas.chr = readFastaMy(paste(path.chr, query.file, sep = ''))
  # query.fas.chr = seq2nt(query.fas.chr)
  # query.len = length(query.fas.chr)
  
  # ---- Read Major ----
  
  pokaz('Read the skeleton alignment..', file=file.log.loop, echo=echo.loop)
  x.sk = readRDS(file.aln.pre)
  
  # TODO: put real lengths, not approximations
  max.chr.len = max(max(x.sk$V4), max(x.sk$V5)) + 10^6
  
  # ---- Read Gaps between normal blocks ----
  
  complexity.threshold = 200  # Max number of blast hits between two synteny blocks
  
  file.gaps.out = paste0(path.gaps,
                         'acc_', acc, 
                         '_qchr_', query.chr, '_bchr_', base.chr, '_out.txt', collapse = '')
  
  pokaz('gap file', file.gaps.out, file=file.log.loop, echo=echo.loop)
  

  if(file.exists(file.gaps.out)){
    pokaz('Read blast of good gaps..', file=file.log.loop, echo=echo.loop)
    x.gap = readBlast(file.gaps.out)
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
    pokaz('nrow', nrow(x.gap), file=file.log.loop, echo=echo.loop)
    if(nrow(x.gap) == 0) x.gap = NULL
  }  # KOSTYL
  
  if(!is.null(x.gap)){  # KOSTYL
    
    pos.beg.info = 2  # Position in sequence's name, where the genome's begin is started
    
    x.gap$q.beg = as.numeric(sapply(x.gap$V1, function(s) strsplit(s, '\\|')[[1]][pos.beg.info])) - 1
    x.gap$b.beg = as.numeric(sapply(x.gap$V10, function(s) strsplit(s, '\\|')[[1]][pos.beg.info])) - 1
    
    x.gap$q.end = as.numeric(sapply(x.gap$V1, function(s) strsplit(s, '\\|')[[1]][pos.beg.info+1]))
    x.gap$b.end = as.numeric(sapply(x.gap$V10, function(s) strsplit(s, '\\|')[[1]][pos.beg.info+1]))
    
    x.gap$V2 = x.gap$V2 + x.gap$q.beg
    x.gap$V3 = x.gap$V3 + x.gap$q.beg
    x.gap$V4 = x.gap$V4 + x.gap$b.beg
    x.gap$V5 = x.gap$V5 + x.gap$b.beg
    x.gap$dir = (x.gap$V4 > x.gap$V5) * 1
    
    x.gap = glueZero(x.gap)
    x.gap$idx = 1:nrow(x.gap)
    
    if(sum(x.gap$V3 > x.gap$q.end) > 0) stop('query')
    if(sum(x.gap$V5 > x.gap$b.end) > 0) stop('base')
    
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
      df.gap = x.tmp[,c('V2', 'V3', 'V4', 'V5', 'idx', 'dir')]
      tmp = df.gap$V4[df.gap$dir == 1]
      df.gap$V4[df.gap$dir == 1] = df.gap$V5[df.gap$dir == 1]
      df.gap$V5[df.gap$dir == 1] = tmp
      df.gap$dir = 0
      
      # Transform positions to start from 1
      x.min = min(df.gap$V2)
      y.min = min(df.gap$V4)
      df.gap[,c('V2', 'V3')] = df.gap[,c('V2', 'V3')] - x.min + 1
      df.gap[,c('V4', 'V5')] = df.gap[,c('V4', 'V5')] - y.min + 1
      df.gap$len.y = df.gap$V5 - df.gap$V4  + 1
      
      # ---- Greedy loop ----
      
      overlap.cutoff = 0.2
      
      # New system of coordinates on Y axis
      y.ticks = sort(unique(c(df.gap$V4, df.gap$V5+1)))
      cell.y.pos = data.frame(beg = y.ticks[-length(y.ticks)], end = y.ticks[-1]-1)
      cell.y.pos$len = cell.y.pos$end - cell.y.pos$beg + 1
      pos.y.attr = rep(0, max(cell.y.pos$end))
      
      df.gap$X4 = 0
      df.gap$X5 = 0
      y.list = list()
      for(irow in 1:nrow(df.gap)){
        pos.tmp = which((cell.y.pos$beg >= df.gap$V4[irow]) & (cell.y.pos$end <= df.gap$V5[irow]))
        df.gap$X4[irow] = min(pos.tmp)
        df.gap$X5[irow] = max(pos.tmp)
      }
      
      df.gap$len.x = df.gap$V3 - df.gap$V2 + 1
      df.gap$len.y = abs(df.gap$V4 - df.gap$V5) + 1
      
      pos.y.occ = rep(0, nrow(cell.y.pos))
      pos.x = 0
      idx.added <- c()
      
      df.gap = df.gap[order(df.gap$V2),]
      
      while(T){
        d.x <- df.gap$V3 - pos.x
        idx.next <- which((d.x > 0) & (d.x / df.gap$len.x > (1-overlap.cutoff)))
        
        
        if(length(idx.next) == 0) break
        
        i.next.found = 0
        for(i.next in idx.next) {
          
          # if(i.next == 5) stop()
          
          pos.y.i.next <- df.gap$X4[i.next]:df.gap$X5[i.next]
          pos.y.i.next.overlap = pos.y.i.next[pos.y.occ[pos.y.i.next] > 0]
          
          y.overlap <- sum(cell.y.pos$len[pos.y.i.next.overlap]) / df.gap$len.y[i.next] 
          if(y.overlap >= overlap.cutoff) next
          
          i.next.found = i.next
          break
        }
        
        if(i.next.found == 0) break
        
        idx.added = c(idx.added, i.next.found)
        pos.y.occ[pos.y.i.next] <- pos.y.occ[pos.y.i.next] + 1
        pos.x <- df.gap$V3[i.next]
        
      }      
      
      # ---- Remaining ----
      idx.good = c(idx.good, df.gap$idx[idx.added])
      
    }
    
    if(length(idx.good) > 0){
      x.res = x.gap[idx.good,]  
      # Clean overlaps from both (base and query) sides
      x.res = cleanOverlaps(x.res)
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
  
  
  # ---- Read additional alignments ----
  
  file.gaps.out = paste0(path.gaps,
                         'acc_', acc, 
                         '_qchr_', query.chr, '_bchr_', base.chr, '_residual_out.txt', collapse = '')
  x.gap = NULL
  if(file.exists(file.gaps.out)){
    # Read blast results on "between blocks"
    pokaz('Read blast of "bad" gaps..', file.gaps.out, file=file.log.loop, echo=echo.loop)
    x.gap = readBlast(file.gaps.out)    
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

  
  if(!is.null(x.gap) & F) {
   
    
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
  pos.q.occup = rep(0, max.chr.len)
  x.sk1 = x.comb
  # x.sk1 = x.res
  # x.sk1 = x.bw
  for(irow in 1:nrow(x.sk1)){
    pp = x.sk1$V4[irow]:x.sk1$V5[irow]
    if(sum(pos.q.occup[pp]) > 0) {
      pokaz('Non-unique', file=file.log.loop, echo=echo.loop)
      stop('non-unique') 
    }
    pos.q.occup[pp] = pos.q.occup[pp] + 1
  }
  sum(pos.q.occup > 1)
  sum(pos.q.occup)
  
  # ---- Check genomes ---- 
  # x.dir = setDir(x.comb, base.len = base.len)
  # checkCorrespToGenome(x.dir, query.fas = query.fas.chr,
  #                      base.fas.fw = base.fas.fw,
  #                      base.fas.bw = base.fas.bw)
  
  saveRDS(object = x.comb, file = file.aln.full) 
  
  suppressWarnings(rm(list = c("x.sk", "x.res", "x.bw", "base.fas.bw", "base.fas.fw", "query.fas.chr")))
  gc()
  
  pokaz('Done.', file=file.log.loop, echo=echo.loop)
  return(NULL)
}
  

# ***********************************************************************
# ---- Loop  ----


if(num.cores == 1){
  file.log.loop = paste0(path.log, 'loop_all.log')
  invisible(file.create(file.log.loop))
  for(f.maj in files.maj){
    loop.function(f.maj,
                  file.log.loop = file.log.loop, 
                  echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(f.maj = files.maj, 
                .packages=c('crayon'))  %dopar% { 
                  loop.function(f.maj, 
                                echo.loop=echo.loop)
                }
  stopCluster(myCluster)
}

warnings()

pokaz('Done.',
      file=file.log.main, echo=echo.main)

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

