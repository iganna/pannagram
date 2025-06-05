# Step 1. Genomes into chromosomes

suppressMessages({
  library("optparse")
  library("foreach")
  library(doParallel)
})

source(system.file("utils/utils.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
  make_option(c("--path.chr"),    type = "character", default = NULL, help = "Path to the input directory"),
  make_option(c("--path.orf"),   type = "character", default = NULL, help = "Path to the output directory"),
  
  make_option(c("--accessions"), type = "character", default = NULL, help = "Files with accessions to analyze"),
  make_option(c("--n.chr"),      type = "integer",   default = 0, help = "Number of chromosomes"),
  
  make_option(c("--cores"),      type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option(c("--path.log"),   type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"),  type = "numeric",   default = 0,    help = "Level of log to be shown on the screen")
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# print(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- Values of other parameters ----

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

# Set input and output paths
path.chr <- ifelse(!is.null(opt$path.chr), opt$path.chr, stop("The input path 'path.in' must be specified!"))
path.orf   <- ifelse(!is.null(opt$path.orf), opt$path.orf, stop("The chromosome-out path 'path.orf' must be specified!"))
if(!dir.exists(path.orf)) dir.create(path.orf)

# ***********************************************************************
# ---- Prepare combinations ----

pokaz('Path with chromosomes:', path.chr, file=file.log.main, echo=echo.main)

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

loop.function <- function(i.comb, echo.loop=T){
  
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
  
  file.in = paste0(path.chr, acc, '_chr', i.chr, '.fasta', collapse = '')
  file.out = paste0(path.orf, acc, '_chr', i.chr, '_orf.fasta', collapse = '')
  
  q.fasta = readFasta(file.in)
  
  pokaz('File out:', file.out,
                 file=file.log.loop, echo=echo.loop)
  
  s = seq2nt(q.fasta)
  s = s[s != 'n']
  s = s[s != 'N']
  s = nt2seq(s) 
  
  orf.res = orfFinder(s)
  
  s.orf = orf.res$orf
  names(s.orf) = paste(acc, '_Chr', i.chr, names(s.orf), sep='')
  
  writeFastaMy(s.orf, file=file.out, append=F)
  
  rm(q.fasta)
  
  pokaz('Done.',
        file=file.log.loop, echo=echo.loop)
  return(NULL)
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

# Rscript query_to_chr.R -n 5 -t fasta --path.in ../pb_genomes/ --path.orf ../pb_chromosomes/
# Rscript query_to_chr.R -n 8 -t fasta --path.in ../lyrata/ --path.orf ../ly_chromosomes/    
# Rscript query_to_chr.R -n 1 -t fasta --path.in ../rhizobia/ --path.orf ../rhiz_chromosomes/ -s T



