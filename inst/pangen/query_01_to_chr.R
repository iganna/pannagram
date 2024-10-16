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
  make_option(c("--path.in"),   type = "character", default = NULL, help = "Path to the input directory"),
  make_option(c("--path.out"),  type = "character", default = NULL, help = "Path to the output directory"),
  make_option(c("--accessions"),type = "character", default = NULL, help = "File with accessions to analyze"),
  
  make_option(c("--n.chr"),     type = "integer",   default = 0,    help = "Number of chromosomes"),
  
  make_option(c("--cores"),     type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option(c("--path.log"),  type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"), type = "integer",   default = NULL, help = "Log level to display on the screen")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# print(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ***********************************************************************
# ---- Values of other parameters ----

# Number of cores
num.cores <- ifelse(!is.null(opt$cores), opt$cores, 30)

# Ensure the number of chromosomes is specified and set it
n.chr <- ifelse(!is.null(opt$n.chr), as.numeric(opt$n.chr), 
                stop("The input number of chromosomes 'n.chr' must be specified!"))  
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
path.query <- ifelse(!is.null(opt$path.in), opt$path.in, stop("The input path 'path.in' must be specified!"))
path.chr   <- ifelse(!is.null(opt$path.out), opt$path.out, stop("The chromosome-out path 'path.out' must be specified!"))
if(!dir.exists(path.chr)) dir.create(path.chr)

# ***********************************************************************
# ---- Preparation ----

# Processor of input genome files
pokaz('Path with genomes:', path.query, file=file.log.main, echo=echo.main)

# Set accepted genome file types
query.types <- c('fasta', 'fna', 'fa', 'fas') #  'ffn', 'faa', 'frn'
pokazAttention('Only the following extensions will be considered:', query.types, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(acc, echo.loop=T){
  
  
  # Check log-files
  file.acc.len = paste0(path.chr, acc, '_chr_len.txt', collapse = '')
  if(file.exists(file.acc.len)){
    chr.len = read.table(file.acc.len, stringsAsFactors = F, header = 1)
    if(n.chr == 0){
      n.log.files = nrow(chr.len)
    } else {
      n.log.files = n.chr
    }
    
    flag.exist = 0
    for(i.chr in 1:n.log.files){
      file.log.loop = paste0(path.log, 'loop_acc_', acc,'_', i.chr, '.log')
      if(file.exists(file.log.loop)){
        if(checkDone(file.log.loop)){
          flag.exist = flag.exist +1
        }
      }
    }
    
    if(flag.exist == n.log.files){
      return(NULL)
    }
    
  }
 
  # ***********************************
  
  # Get the genome file
  file.genome = c()
  for(s.ext in query.types){
    file.genome = paste0(path.query, acc, '.', s.ext, collapse = '')
    if(file.exists(file.genome)) break
  }
  if(!file.exists(file.genome)){
    pokazAttention('Genome file for the accession', acc, 'does not exist')
    return(NULL)
  }

  pokaz('Accession', acc, file=file.log.main, echo=echo.loop)
  
  # Read the genome file
  q.fasta = readFastaMy(file.genome)
  
  # Check the chromosome number
  if(n.chr == 0){ # if to analyse all chromosomes
    n.chr = length(q.fasta)
  }
  if(length(q.fasta) < n.chr){
    pokazAttention('Accession', acc, 'was not analysed, not enough chromosomes in the genome.\n
                   Exist:', length(q.fasta), 'Requeired:', n.chr, 
                   file=file.log.main, echo=echo.loop)
    return(NULL)
  }
  
  # Save chromosomal lengths
  df.chr.lengths = data.frame(acc = acc,
                              chr = 1:length(q.fasta),
                              len = nchar(q.fasta), 
                              name = gsub(" ", "_", names(q.fasta)))
  df.chr.lengths = df.chr.lengths[1:n.chr,]
  write.table(df.chr.lengths, file.acc.len, sep = '\t', col.names = T, row.names = F, quote = F)
  
  pokaz('Chromosomes', names(q.fasta)[1:(n.chr)], 'will be processed',
                 file=file.log.main, echo=echo.loop)
  
  if(length(q.fasta) > n.chr){
    pokaz('Chromosomes', names(q.fasta)[(n.chr+1):length(q.fasta)], 'will NOT be processed',
                   file=file.log.main, echo=echo.loop)
  }
  
  # Write every chromosome to a separate file
  for(i.chr in 1:n.chr){
    
    # ---- Log files ----
    file.log.loop = paste0(path.log, 'loop_acc_', acc,'_', i.chr, '.log')
    if(!file.exists(file.log.loop)){
      invisible(file.create(file.log.loop))
    }
    
    # Check log Done
    if(checkDone(file.log.loop)){
      next
    }
    
    # ---- Write chromosome ----
    file.out = paste0(path.chr, acc, '_chr', i.chr, '.fasta', collapse = '')

    pokaz('File out:', file.out,
                   file=file.log.loop, echo=echo.loop)
    
    writeFastaMy(toupper(q.fasta[i.chr]), 
                 file=file.out, append=F, 
                 seq.names = paste0(acc, '_Chr', i.chr ))
    
    pokaz('Done.',
          file=file.log.loop, echo=echo.loop)
    
  }
  
  rm(q.fasta)
  return(NULL)
}
  
# ***********************************************************************
# ---- Loop  ----


if(num.cores == 1){
  for(acc in accessions){
    loop.function(acc, echo.loop=echo.loop)
  }
} else {
  # Initialise clusters
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(acc = accessions, 
                .packages=c('crayon'),
                .export = c('n.chr')) %dopar% {
                                     loop.function(acc, echo.loop=echo.loop)
                                   }
  stopCluster(myCluster)
}

warnings()

pokaz('Done.',
      file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- Manual testing ----

