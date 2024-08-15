# Step 1. Genomes into chromosomes

suppressMessages({
  library("optparse")
  library("foreach")
  library(doParallel)
})

source("utils/utils.R")

# ***********************************************************************
# ---- Command line arguments ----

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
  make_option(c("--all.chr"), type = "logical", default = FALSE, 
              help = "Flag to use all chromosomes, not only the provided number", metavar = "logical"),
  make_option(c("--n.chr"), type = "numeric", default = NULL, 
              help = "Number of chromosomes", metavar = "numeric"),
  make_option(c("--path.in"), type = "character", default = NULL, 
              help = "Path to the input directory", metavar = "character"),
  make_option(c("--path.out"), type = "character", default = NULL, 
              help = "Path to the output directory", metavar = "character"),
  make_option(c("--sort"), type = "logical", default = FALSE, 
              help = "Sort chromosomes by lengths or not", metavar = "logical"),
  make_option(c("--cores"), type = "integer", default = 1, 
              help = "Number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("--accessions"), type = "character", default = NULL,
              help = "Files with accessions to analyze", metavar = "character"),
  make_option(c("--path.log"), type = "character", default = NULL,
              help = "Path for log files", metavar = "character"),
  make_option(c("--log.level"), type = "numeric", default = 0,
              help = "Level of log to be shown on the screen", metavar = "numeric")
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

print(opt)

# ***********************************************************************
# ---- Logging ----

source('utils/chunk_logging.R') # a common code for all R logging

# ---- Values of other parameters ----

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

# Accessions to analyse
acc.anal <- opt$accessions
if(!is.null(acc.anal)){
  if(acc.anal == 'NULL'){
    acc.anal <- NULL
  } else {
    if (!file.exists(acc.anal)) {
      pokazAttention('File', acc.anal, 'does NOT exists, so no accession filtration is applied.', 
                     file=file.log.main, echo=echo.main)
      acc.anal <- NULL
    } else {
      tmp <- read.table(acc.anal, stringsAsFactors = F)
      acc.anal <- tmp[,1]
    }
  }
}

# Set input and output paths
path.query <- ifelse(!is.null(opt$path.in), opt$path.in, stop("The input path 'path.in' must be specified!"))
path.chr   <- ifelse(!is.null(opt$path.out), opt$path.out, stop("The chromosome-out path 'path.out' must be specified!"))
if(!dir.exists(path.chr)) dir.create(path.chr)

# Decide whether to sort by length based on provided input or default to FALSE
sort.by.lengths <- ifelse(!is.null(opt$sort), as.logical(opt$sort), FALSE)

# ***********************************************************************
# ---- Sort chromosomal lengths ----

# if(!sort.by.lengths){
#   msg = 'If you use -one2one option: please be sure that all chromosomes in files are sorted in the same order' # or use \"-s T\" flag'
#   pokazAttention(msg, file=file.log.main, echo=echo.main)
# } else {
#   msg = 'Chromosomes will be sorted by their length'
#   pokazAttention(msg, file=file.log.main, echo=echo.main)
# }

# ***********************************************************************
# ---- Preparation ----

# Processor of input genome files
pokaz('Path with genomes:', path.query, file=file.log.main, echo=echo.main)

# Set accepted genome file types
query.types <- c('fasta', 'fna', 'fa', 'fas') #  'ffn', 'faa', 'frn'
pokazAttention('Only the following extensions will be considered:', query.types, file=file.log.main, echo=echo.main)

# List and filter genome files in the specified path based on the accepted file types
search.pattern <- paste0(".*\\.(?:", paste(query.types, collapse="|"), ")$")
query.name <- basename(list.files(path.query, pattern = search.pattern, full.names = TRUE))
if(length(query.name) == 0) stop('No accessions is provided for the analysys.')
query.name <- data.frame(file=query.name, acc=gsub("(\\.[^.]+)$", "", query.name), stringsAsFactors = F)

# Optional: Filter based on a list of accession numbers, if provided
pokaz('query.name', query.name)
pokaz('acc.anal', acc.anal)
if(!is.null(acc.anal)){
  acc.anal <- gsub("(\\.[^.]+)$", "", basename(acc.anal))
  # pokaz('Accessions in the folder', query.name, file=file.log.main, echo=echo.main)
  # pokaz('Accessions in the filtration file', acc.anal, file=file.log.main, echo=echo.main)
  query.name <- query.name[query.name$acc %in% acc.anal,,drop=F]
}

# Final check and display of genome names for analysis
if((nrow(query.name) == 0) || (length(query.name) == 0)) stop('No accessions is provided for the analysys.')
pokaz('Names of genomes for the analysis:', query.name$acc, 
      file=file.log.main, echo=echo.main)


# ***********************************************************************
# ---- MAIN program body ----

file.log.pref <- ''

loop.function <- function(i.acc, 
                          echo.loop=T, 
                          file.log.loop=NULL){
  
  acc = query.name$acc[i.acc]
  
  # Log files
  if (is.null(file.log.loop)){
    file.log.loop = paste0(path.log, 'loop_', i.acc, '_acc_', acc,'.log')
    invisible(file.create(file.log.loop))
  }
  
  #' --- --- --- --- --- --- --- --- --- --- ---
  
  # Don't run if the chromosomal files exist
  if(!all.chr){  # if the chromosome number is an important parameter, not ALL_CHROMOSOMES
    n.exist = 0
    for(i.chr in 1:n.chr){
      file.out = paste0(path.chr, acc, '_chr', i.chr, '.fasta', collapse = '')
      if(file.exists(file.out)){
        n.exist = n.exist + 1
      }
    }
    # pokaz(n.exist, n.chr,
    #       file=file.log.loop, echo=echo.loop)
    if(n.exist == n.chr){  # If chromosomal files are already formed
      return(NULL)
    }
  }
 
  # ***********************************
  # When no chromosome-files were produced before
  
  pokaz('Accession', acc,
        file=file.log.loop, echo=echo.loop)
  q.fasta = readFastaMy(paste0(path.query, query.name$file[i.acc], collapse = ''))
  
  
  if(all.chr){ # if to analyse all chromosomes
    n.chr = length(q.fasta)
  }
  if(length(q.fasta) < n.chr){
    pokazAttention('Accession', acc, 'was not analysed, not enough chromosomes in the genome.\n
                   Exist:', length(q.fasta), 'Requeired:', n.chr, 
                   file=file.log.loop, echo=echo.loop)
    return(NULL)
  }
  
  # Chromosomal lengths
  file.acc.len = paste0(path.chr, acc, '_chr_len.txt', collapse = '')
  df.chr.lengths = data.frame(acc = acc,
                              chr = 1:length(q.fasta),
                              len = nchar(q.fasta), 
                              name = gsub(" ", "_", names(q.fasta)))
  df.chr.lengths = df.chr.lengths[1:n.chr,]
  write.table(df.chr.lengths, file.acc.len, sep = '\t', col.names = T, row.names = F, quote = F)
  
  
  if(sort.by.lengths){
    q.len = sapply(q.fasta, nchar)
    pokaz('Lengths', q.len)
    q.fasta = q.fasta[order(-q.len)]
  }
  
  pokaz('Chromosomes', names(q.fasta)[1:(n.chr)], 'will be processed',
                 file=file.log.loop, echo=echo.loop)
  
  if(length(q.fasta) > n.chr){
    pokaz('Chromosomes', names(q.fasta)[(n.chr+1):length(q.fasta)], 'will NOT be processed',
                   file=file.log.loop, echo=echo.loop)
  }
  
  for(i.chr in 1:n.chr){
    # acc.s = gsub('_', '-', acc)
    acc.s = acc
    file.out = paste0(path.chr, acc.s, '_chr', i.chr, '.fasta', collapse = '')
    if(file.exists(file.out)){
      next
    }
    pokaz('File out:', file.out,
                   file=file.log.loop, echo=echo.loop)
    
    s = toupper(q.fasta[i.chr])
    writeFastaMy(s, file=file.out, append=F, seq.names = paste(acc.s, '_Chr', i.chr , sep=''))
  }
  
  rm(q.fasta)
  
  pokaz('Done.',
        file=file.log.loop, echo=echo.loop)
  return(NULL)
}
  
# ***********************************************************************
# ---- Loop  ----
  

if(num.cores == 1){
  file.log.loop = paste0(path.log, 'loop_all.log')
  invisible(file.create(file.log.loop))
  for(i.acc in 1:nrow(query.name)){
    loop.function(i.acc, 
                  file.log.loop = file.log.loop, 
                  echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(i.acc = 1:nrow(query.name), 
                .packages=c('crayon'),
                .export = c('n.chr')) %dopar% {
                                     loop.function(i.acc, 
                                                   echo.loop=echo.loop)
                                   }
  stopCluster(myCluster)
}

warnings()

pokaz('Done.',
      file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- Manual testing ----

# Rscript query_to_chr.R -n 5 -t fasta --path.in ../pb_genomes/ --path.out ../pb_chromosomes/
# Rscript query_to_chr.R -n 8 -t fasta --path.in ../lyrata/ --path.out ../ly_chromosomes/    
# Rscript query_to_chr.R -n 1 -t fasta --path.in ../rhizobia/ --path.out ../rhiz_chromosomes/ -s T



