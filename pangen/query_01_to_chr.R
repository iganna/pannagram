suppressMessages({
  library("optparse")
  # library(Biostrings)
  # library("seqinr")       # read.fasta
  library("foreach")
  library(doParallel)
})

source("utils/utils.R")

# pokazStage('Step 1. Genomes into chromosomes')

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("--all.chr"), type = "character", default = NULL, 
              help = "Flag to use all chromosomes, not only the provided number", metavar = "character"),
  make_option(c("--n.chr"), type = "character", default = NULL, 
              help = "number of chromosomes", metavar = "character"),
  make_option(c("--path.in"), type = "character", default = NULL, 
              help = "pathway to the input directory", metavar = "character"),
  make_option(c("--path.out"), type = "character", default = NULL, 
              help = "pathway to the output directory", metavar = "character"),
  make_option(c("--sort"), type = "character", default = NULL, 
              help = "sort chromosomes by lengths or not", metavar = "character"),
  make_option(c("--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("--acc.anal"), type = "character", default = NULL,
              help = "files with accessions to analyze", metavar = "character")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# print(opt)

# ***********************************************************************
# ---- Values of parameters ----

# Number of cores
num.cores <- ifelse(!is.null(opt$cores), opt$cores, 30)

# Ensure the number of chromosomes is specified and set it
all.chr <- ifelse(!is.null(opt$all.chr), as.logical(opt$all.chr), F)
if(all.chr){
  n.chr = NULL
} else {
  n.chr <- ifelse(!is.null(opt$n.chr), as.numeric(opt$n.chr), stop("The input number of chromosomes 'n.chr' must be specified!"))  
}

# Accessions to analyse
acc.anal <- opt$acc.anal

if(!is.null(acc.anal)){
  if (!file.exists(acc.anal)) {
    acc.anal = NULL
    pokazAttention('File', acc.anal, 'does NOT exists, so no accession filtration is applied.')
  } else {
    tmp = read.table(acc.anal, stringsAsFactors = F)
    acc.anal = tmp[,1]
  }
}

pokaz('Number of chromosomes:', n.chr)

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
#   pokazAttention(msg)
# } else {
#   msg = 'Chromosomes will be sorted by their length'
#   pokazAttention(msg)
# }

# ***********************************************************************
# ---- Preparation ----

# Processor of input genome files
pokaz('Path with genomes:', path.query)

# Set accepted genome file types
query.types = c('fasta', 'fna', 'fa', 'fas') #  'ffn', 'faa', 'frn'
pokazAttention('Only the following extensions will be considered:', query.types)

# List and filter genome files in the specified path based on the accepted file types
search.pattern <- paste0(".*\\.(?:", paste(query.types, collapse="|"), ")$")
query.name <- basename(list.files(path.query, pattern = search.pattern, full.names = TRUE))
if(length(query.name) == 0) stop('No accessions is provided for the analysys.')
query.name <- data.frame(file=query.name, acc=gsub("(\\.[^.]+)$", "", query.name))

# Optional: Filter based on a list of accession numbers, if provided
if(!is.null(acc.anal)){
  acc.anal = gsub("(\\.[^.]+)$", "", basename(acc.anal))
  # pokaz('Accessions in the folder', query.name)
  # pokaz('Accessions in the filtration file', acc.anal)
  query.name = query.name[query.name$acc %in% acc.anal,,drop=F]
}

# Final check and display of genome names for analysis
if((nrow(query.name) == 0) || (length(query.name) == 0)) stop('No accessions is provided for the analysys.')
pokaz('Names of genomes for the analysis:', query.name$acc)


# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(i.acc, echo = T){
  
  acc = query.name$acc[i.acc]
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
    if(echo) pokaz(n.exist, n.chr)
    if(n.exist == n.chr){  # If chromosomal files are already formed
      # if(for.flag) next
      return(NULL)
    }
  }
 
  # ***********************************
  
  if(echo) pokaz('Accession', acc)
  q.fasta = readFastaMy(paste0(path.query, query.name$file[i.acc], collapse = ''))
  
  if(all.chr){ # if to analyse all chromosomes
    n.chr = length(q.fasta)
  }
  if(length(q.fasta) < n.chr){
    pokazAttention('Accession', acc, 'was not analysed, not enough chromosomes in the genome.\n
                   Exist:', length(q.fasta), 'Requeired:', n.chr)
    return(NULL)
  }
  
  if(sort.by.lengths){
    q.len = sapply(q.fasta, nchar)
    if(echo) pokaz('Lengths', q.len)
    q.fasta = q.fasta[order(-q.len)]
  }
  
  if(echo) pokaz('Chromosomes', names(q.fasta)[1:(n.chr)], 'will be processed')
  
  if(length(q.fasta) > n.chr){
    if(echo) pokaz('Chromosomes', names(q.fasta)[(n.chr+1):length(q.fasta)], 'will NOT be processed')
  }
  
  for(i.chr in 1:n.chr){
    # acc.s = gsub('_', '-', acc)
    acc.s = acc
    file.out = paste0(path.chr, acc.s, '_chr', i.chr, '.fasta', collapse = '')
    if(file.exists(file.out)){
      # if(for.flag) next
      return(NULL)
    }
    if(echo) pokaz('File out', file.out)
    
    s = toupper(q.fasta[i.chr])
    writeFastaMy(s, file=file.out, append=F, seq.names = paste(acc.s, '_Chr', i.chr , sep=''))
    
    # write(paste('>', acc.s, '_Chr', i.chr , sep=''), file=file.out, append=F)
    # write(s, file=file.out, append=T)
    # write('\n', file=file.out, append=T)
  }
  
  rm(q.fasta)
}
  
# ***********************************************************************
# ---- Loop  ----
  

if(num.cores == 1){
  for(i.acc in 1:nrow(query.name)){
    loop.function(i.acc)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(i.acc = 1:nrow(query.name), 
                .packages=c('crayon'),
                .export = c('n.chr')) %dopar% {
                                     loop.function(i.acc)
                                   }
  stopCluster(myCluster)
}

warnings()

# ***********************************************************************
# ---- Manual testing ----

# Rscript query_to_chr.R -n 5 -t fasta --path.in ../pb_genomes/ --path.out ../pb_chromosomes/
# Rscript query_to_chr.R -n 8 -t fasta --path.in ../lyrata/ --path.out ../ly_chromosomes/    
# Rscript query_to_chr.R -n 1 -t fasta --path.in ../rhizobia/ --path.out ../rhiz_chromosomes/ -s T



