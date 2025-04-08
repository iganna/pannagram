# Step 2. Chromosomes into parts

suppressMessages({
  library("optparse")
  library("foreach")
  library(doParallel)
})

source(system.file("utils/utils.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----
args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("--path.chr"),   type = "character", default = NULL, help = "Path to the chromosome directory"),
  make_option(c("--path.parts"), type = "character", default = NULL, help = "Path to the parts directory"),
  make_option(c("--accessions"), type = "character", default = NULL, help = "File containing accessions to analyze"),
  
  make_option(c("--n.chr"),      type = "integer", default = 0,    help = "Number of chromosomes"),
  make_option(c("--part.len"),   type = "integer", default = 5000, help = "Length of each part file in bp"),
  make_option(c("--part.step"),  type = "integer", default = 0,    help = "Step size in bp between parts"),
  
  make_option(c("--purge.reps"), type = "logical", default = FALSE, help = "Flag to specify whether to remove repeats"),
  make_option(c("--rev"),        type = "logical", default = FALSE, help = "Flag to reverse sequences"),
  
  make_option(c("--cores"),     type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option(c("--path.log"),  type = "character", default = NULL, help = "Path for the log files"),
  make_option(c("--log.level"), type = "integer", default = NULL, help = "Logging level to display on screen")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# pokaz(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- Values of parameters ----

# Set chromosome and parts paths
path.chr <- ifelse(!is.null(opt$path.chr), opt$path.chr, stop("The chromosome path 'path.chr' must be specified!"))
path.parts <- ifelse(!is.null(opt$path.parts), opt$path.parts, stop("The parts path 'path.parts' must be specified!"))
if(!dir.exists(path.parts)) dir.create(path.parts)

# Number of cores
num.cores <- opt$cores

# Ensure the number of chromosomes is specified and set it
n.chr <- opt$n.chr  
if(n.chr == 0){
  pokaz("All chromosomes are considered", file=file.log.main, echo=echo.main)
} else {
  pokaz('Number of chromosomes:', n.chr, file=file.log.main, echo=echo.main)  
}

# Accessions to analyse
file.acc <- ifelse(!is.null(opt$accessions), opt$accessions, stop("File with accessions are not specified"))
tmp <- read.table(file.acc, stringsAsFactors = F)
accessions <- as.character(tmp[,1])
pokaz('Names of genomes for the analysis:', accessions, 
      file=file.log.main, echo=echo.main)

# ---- Common attributes ----

# Length of parts
len.parts <- opt$part.len
len.step <- opt$part.step

# Purge repeats by the complexity
purge.reps = opt$purge.reps

# Mirror universe
flag.rev <- opt$rev
if(flag.rev){
  pokaz("Mirror universe!!!", 
        file=file.log.main, echo=echo.main)
}

# ***********************************************************************
# ---- Prepare combinations ----

if(n.chr == 0){
  # for every accession read the "length" file and cound the nu,ber of chromosomes
  combinations = c()
  for(acc in accessions){
    file.acc.len = paste0(path.chr, acc, '_chr_len.txt', collapse = '')
    acc.len = read.table(file.acc.len, header = 1)
    
    combinations = rbind(combinations,
                         data.frame(acc = acc, i.chr = 1:nrow(acc.len)))
  }
} else {
  combinations <- expand.grid(acc = accessions, i.chr = 1:n.chr)  
}

# ***********************************************************************
# ---- MAIN program body ----
loop.function <- function(i.comb, 
                          echo.loop=T, 
                          file.log.loop=NULL){
  
  acc <- combinations$acc[i.comb]
  i.chr <- combinations$i.chr[i.comb]
  
  # Log files
  file.log.loop = paste0(path.log, 'loop_acc_', acc, '_chr', i.chr, '.log')
  if(!file.exists(file.log.loop)){
    invisible(file.create(file.log.loop))
  }
  
  # Check log Done
  if(checkDone(file.log.loop)){
    return(NULL)
  }
  
  # ***********************************
  pokaz('New attempt:', file=file.log.loop, echo=echo.loop)
  
  # Files
  file.in = paste0(path.chr, acc, '_chr', i.chr, '.fasta', collapse = '')
  file.out = paste0(path.parts, acc, '_chr', i.chr, '.fasta', collapse = '')
  
  pokaz('File chromosomal:', file.in, 
        file=file.log.loop, echo=echo.loop)
  if(!file.exists(file.in)){
    pokaz('Chr file not exist', file=file.log.loop, echo=echo.loop)
    return(NULL)
  }
  q.fasta = readFastaMy(file.in)[1]
  q.fasta = toupper(q.fasta)
  
  if(flag.rev){  # Mirror universe
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
  
  if(purge.reps){  # Filter out repeats
    pokaz('Filterout repeats', file=file.log.loop, echo=echo.loop)
    
    file.out = paste0(path.parts, acc, '_chr', i.chr, '.fasta', collapse = '')
    
    file.out.rest = paste0(path.parts, acc, '_chr', i.chr, '.rest', collapse = '')
    file.out.masking = paste0(path.chr, 'mask_', acc, '_chr', i.chr, '.rds', collapse = '')
    
    seqs.score = sapply(s, repeatScore)
    
    if (sum(seqs.score <= 0.2) > 0){
      writeFastaMy(s[seqs.score <= 0.2], file.out)  
    } else {
      pokaz('No good sequences sequences.', file=file.log.loop, echo=echo.loop)
    }
    
    if (sum(seqs.score > 0.2) > 0){
      
      s.repeat = s[seqs.score > 0.2]
      writeFastaMy(s.repeat, file.out.rest)
      
      # Masking positions
      pos.repeat = as.numeric(sapply(names(s.repeat), function(s) strsplit(s, '\\|')[[1]][4]))
      pos.masking = data.frame(beg = pos.repeat, end = pos.repeat + len.parts - 1)
      
      saveRDS(pos.masking, file.out.masking)
      
    } else {
      pokaz('No rest sequences.', file=file.log.loop, echo=echo.loop)
    }
    
  } else {
    file.out = paste0(path.parts, acc, '_chr', i.chr, '.fasta', collapse = '')
    writeFastaMy(s, file.out)
    
  }
  
  rmSafe(q.fasta)
  rmSafe(s)
  rmSafe(pos.beg)
  rmSafe(seqs.score)
  pokaz('Done.', file=file.log.loop, echo=echo.loop)
}

# ***********************************************************************
# ---- Loop  ----

if(num.cores == 1){
  for(i.comb in 1:nrow(combinations)){
    loop.function(i.comb, echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK")
  registerDoParallel(myCluster)
  
  tmp.output = foreach(i.comb = 1:nrow(combinations), 
                       .packages=c('crayon',
                                   'stringi'),  # for purging repeats
                       .export = c('n.chr')) %dopar% {
    loop.function(i.comb, echo.loop=echo.loop)
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










