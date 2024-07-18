# Step 2. Chromosomes into parts

suppressMessages({
  library("optparse")
  library("foreach")
  library(doParallel)
})

source("utils/utils.R")

# ***********************************************************************
# ---- Command line arguments ----
args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("--all.chr"), type = "logical", default = FALSE, 
              help = "Flag to use all chromosomes, not only the provided number", metavar = "logical"),
  make_option(c("--n.chr"), type = "character", default = NULL, 
              help = "number of chromosomes", metavar = "character"),
  make_option(c("--part.len"), type = "character", default = NULL, 
              help = "number of base pairs in the part file", metavar = "character"),
  make_option(c("--part.step"), type = "character", default = NULL, 
              help = "number of base pairs as a step", metavar = "character"),
  make_option(c("--path.chr"), type = "character", default = NULL, 
              help = "pathway to the chromosome directory", metavar = "character"),
  make_option(c("--path.parts"), type = "character", default = NULL, 
              help = "pathway to the parts directory", metavar = "character"),
  make_option(c("--purge_reps"), type = "character", default = NULL, 
              help = "flag to keep or not repeats", metavar = "character"),
  make_option(c("--rev"), type = "character", default = NULL, 
              help = "flag make the reverce sequences", metavar = "character"),
  make_option(c("--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("--path.log"), type = "character", default = NULL,
              help = "Path for log files", metavar = "character"),
  make_option(c("--log.level"), type = "character", default = NULL,
              help = "Level of log to be shown on the screen", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# pokaz(opt)

# ***********************************************************************
# ---- Logging ----

source('utils/chunk_logging.R') # a common code for all R logging

# ---- Values of parameters ----

# Number of cores
num.cores <- ifelse(!is.null(opt$cores), opt$cores, 30)

# Ensure the number of chromosomes is specified and set it
all.chr <- ifelse(!is.null(opt$all.chr), as.logical(opt$all.chr), F)
if(all.chr){
  n.chr <- NULL
} else {
  n.chr <- ifelse(!is.null(opt$n.chr), as.numeric(opt$n.chr), 
                  stop("The input number of chromosomes 'n.chr' must be specified!"))  
}
pokaz('Number of chromosomes:', n.chr, file=file.log.main, echo=echo.main)

# Set chromosome and parts paths
path.chr <- ifelse(!is.null(opt$path.chr), opt$path.chr, stop("The chromosome path 'path.chr' must be specified!"))
path.parts <- ifelse(!is.null(opt$path.parts), opt$path.parts, stop("The parts path 'path.parts' must be specified!"))
if(!dir.exists(path.parts)) dir.create(path.parts)

# Common attributes
len.parts <- ifelse(!is.null(opt$part.len), as.numeric(opt$part.len), 5000)

if(!is.null(opt$part.step)){
  len.step = as.numeric(opt$part.step)
} else {
  len.step = NULL
}
# len.step <- ifelse(!is.null(opt$part.step), as.numeric(opt$part.step), NULL)

# Purge repeats by the complexity
if(is.null(opt$purge_reps)){
  purge_reps = F
} else {
  purge_reps = as.logical(opt$purge_reps)
  if(is.na(purge_reps)) stop('Wrong flag for purging repeats')
}


flag.rev <- as.numeric(ifelse(!is.null(opt$rev), opt$rev, 0))

if(flag.rev == 1){
  pokaz("Mirror universe", 
        file=file.log.main, echo=echo.main)
}

# ***********************************************************************
# ---- Preparation ----

pokaz('Path with chromosomes:', path.chr, 
      file=file.log.main, echo=echo.main)
files.query = list.files(path = path.chr, pattern = paste0('\\.', 'fasta', '$', collapse = '') )
files.query <- sub("\\.fasta$", "", files.query)
query.name = unique(sapply(files.query, function(s) strsplit(s, '_chr')[[1]][1]))
pokaz('Names of genomes for the analysis:', query.name, 
      file=file.log.main, echo=echo.main)

if(length(query.name) == 0){
  stop('Wrong names of chromosomal files or files are not provided')
}

if(all.chr){
  combinations = data.frame(acc = sapply(files.query, function(s) strsplit(s, '_chr')[[1]][1]),
                            i.chr = sapply(files.query, function(s) strsplit(s, '_chr')[[1]][2]))
} else {
  combinations <- expand.grid(acc = query.name, i.chr = 1:n.chr)  
}



# ***********************************************************************
# ---- MAIN program body ----
loop.function <- function(i.comb, 
                          echo.loop=T, 
                          file.log.loop=NULL){
  
  acc <- combinations$acc[i.comb]
  i.chr <- combinations$i.chr[i.comb]
  
  # Log files
  if (is.null(file.log.loop)){
    file.log.loop = paste0(path.log, 'loop_', i.comb, '_acc_', acc,'.log')
    invisible(file.create(file.log.loop))
  }
  
  #' --- --- --- --- --- --- --- --- --- --- ---
  
  file.in = paste0(path.chr, acc, '_chr', i.chr, '.fasta', collapse = '')
  
  file.out = paste0(path.parts, acc, '_chr', i.chr, '.fasta', collapse = '')
  if( file.exists(file.out)) {
    return(NULL)
  }
  pokaz('File:', file.in, 
        file=file.log.loop, echo=echo.loop)
  q.fasta = readFastaMy(file.in)[1]
  q.fasta = toupper(q.fasta)
  
  if(flag.rev == 1){  # Mirror universe
    q.fasta = nt2seq(rev(seq2nt(q.fasta)))
  }
  
  s = splitSeq(q.fasta, n=len.parts)
  len.chr = nchar(q.fasta)
  pos.beg = seq(1, len.chr, len.parts)
  
  if(!is.null(len.step)){
    s = c(s, 
          splitSeq(q.fasta, n=len.parts, step = len.step))
    pos.beg = c(pos.beg, 
                seq(1 + len.step, len.chr, len.parts))
    
    idx.order = order(pos.beg)
    s = s[idx.order]
    pos.beg = pos.beg[idx.order]
  }
  
  if(length(s) != length(pos.beg)) {
    pokaz('Problem with chunks', 
          file=file.log.loop, echo=echo.loop)
    stop('Problem with chunks') 
  }
  names(s) = paste('acc_', acc, '|chr_', i.chr, '|part_', 1:length(s), '|', pos.beg, sep='')
  
  if(purge_reps){  # Filter out repeats
    file.out = paste0(path.parts, acc, '_chr', i.chr, '.fasta', collapse = '')
    
    file.out.rest = paste0(path.parts, acc, '_chr', i.chr, '.rest', collapse = '')
    
    seqs.score = sapply(s, repeatScore)
    
    writeFastaMy(s[seqs.score <= 0.2], file.out)
    writeFastaMy(s[seqs.score > 0.2], file.out.rest)
  } else {
    file.out = paste0(path.parts, acc, '_chr', i.chr, '.fasta', collapse = '')
    writeFastaMy(s, file.out)
    
  }
  
  rmSafe(q.fasta)
  rmSafe(s)
  rmSafe(pos.beg)
  rmSafe(seqs.score)
  pokaz('Done.',
        file=file.log.loop, echo=echo.loop)
}

# ***********************************************************************
# ---- Loop  ----

if(num.cores == 1){
  file.log.loop = paste0(path.log, 'loop_all.log')
  invisible(file.create(file.log.loop))
  for(i.comb in 1:nrow(combinations)){
    loop.function(i.comb, 
                  file.log.loop = file.log.loop, 
                  echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK")
  registerDoParallel(myCluster)
  
  tmp.output = foreach(i.comb = 1:nrow(combinations), 
                       .packages=c('crayon'),
                       .export = c('n.chr')) %dopar% {
    loop.function(i.comb,
                  echo.loop=echo.loop)
  }
  stopCluster(myCluster)
}

warnings()

pokaz('Done.',
      file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- Manual testing ----

if(F){
  len.parts = 5000  # lengths of parts
  n.chr = 5  # number of chromosomes
  path.chr = '../pb_databases/'
  path.parts = '../pb_parts/'
}

#Rscript query_to_parts.R -n 5 -t fasta --path.chr ../pb_chromosomes/ -b 5000 --path.parts ../pb_parts/
# Rscript query_to_parts.R -n 8 -t fasta --path.chr ../ly_chromosomes/ -b 5000 --path.parts ../ly_parts/
# Rscript query_to_parts.R -n 1 -t fasta --path.chr ../rhiz_chromosomes/ -b 5000 --path.parts ../rhiz_parts/ 










