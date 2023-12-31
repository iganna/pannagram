suppressMessages({
  library("optparse")
  library(Biostrings)
  # library("seqinr")       # read.fasta
  library("foreach")
  library(doParallel)
})

source("utils/utils.R")

pokazStage('Step 1. Genomes into chromosomes')

# Rscript query_to_chr.R -n 5 -t fasta --path.in ../pb_genomes/ --path.out ../pb_chromosomes/
# Rscript query_to_chr.R -n 8 -t fasta --path.in ../lyrata/ --path.out ../ly_chromosomes/    
# Rscript query_to_chr.R -n 1 -t fasta --path.in ../rhizobia/ --path.out ../rhiz_chromosomes/ -s T

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("-n", "--n.chr"), type = "character", default = NULL, 
              help = "number of chromosomes", metavar = "character"),
  make_option(c("-i", "--path.in"), type = "character", default = NULL, 
              help = "pathway to the input directory", metavar = "character"),
  make_option(c("-o", "--path.out"), type = "character", default = NULL, 
              help = "pathway to the output directory", metavar = "character"),
  make_option(c("-s", "--sort"), type = "character", default = NULL, 
              help = "sort chromosomes by lengths or not", metavar = "character"),
  make_option(c("-b", "--bp"), type = "character", default = NULL, 
              help = "number of base pairs in the part file", metavar = "character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("-a", "--acc.anal"), type = "character", default = NULL,
              help = "what axes to analyze", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#' ----------------------------------------------------------------------

# Set the number of cores for parallel processing
num_cores <- ifelse(!is.null(opt$cores), opt$cores, 30)
myCluster <- makeCluster(num_cores, type = "PSOCK") 
registerDoParallel(myCluster) 

# Ensure the number of chromosomes is specified and set it
n.chr <- ifelse(!is.null(opt$n.chr), as.numeric(opt$n.chr), stop("The input number of chromosomes 'n.chr' must be specified!"))
acc.anal <- ifelse(!is.null(opt$acc.anal), (opt$acc.anal), NULL)
if(acc.anal == 'NULL') acc.anal = NULL
if(!is.null(acc.anal)){
  if (!file.exists(acc.anal)) {
    acc.anal = NULL
    pokazAttention('File', acc.anal, 'does NOT exists, so no accession filtration is applied.')
  } else {
    tmp = read.table(acc.anal, stringsAsFactors = F)
    acc.anal = tmp[,1]
  }
}

# Set input and output paths
path.query <- ifelse(!is.null(opt$path.in), opt$path.in, stop("The input path 'path.in' must be specified!"))
path.chr   <- ifelse(!is.null(opt$path.out), opt$path.out, stop("The chromosome-out path 'path.out' must be specified!"))
if(!dir.exists(path.chr)) dir.create(path.chr)

# Common attributes
len.parts  <- ifelse(!is.null(opt$bp), as.numeric(opt$bp), 5000)

# Decide whether to sort by length based on provided input or default to FALSE
sort.by.lengths <- ifelse(!is.null(opt$sort), as.logical(opt$sort), FALSE)
pokaz('sort_chr_by_length', sort.by.lengths)


#' ----------------------------------------------------------------------

if(!sort.by.lengths){
  msg = 'Please be sure that all chromosomes in files are sorted in the same order' # or use \"-s T\" flag'
  pokazAttention(msg)
} else {
  msg = 'Chromosomes will be sorted by their length'
  pokazAttention(msg)
}

#' ----------------------------------------------------------------------
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
  pokaz('Accessions in the folder', query.name)
  pokaz('Accessions in the filtration file', acc.anal)
  query.name = query.name[query.name$acc %in% acc.anal]
}

# Final check and display of genome names for analysis
if(nrow(query.name) == 0) stop('No accessions is provided for the analysys.')
pokaz('Names of genomes for the analysis:', query.name$acc)

#' ----------------------------------------------------------------------

for.flag = F
tmp = foreach(i.acc = 1:nrow(query.name), .packages=c('stringr','Biostrings', 'seqinr', 'crayon')) %dopar% {
# for.flag = T
# for(i.acc in 1:nrow(query.name)){
  
  acc = query.name$acc[i.acc]
  #' --- --- --- --- --- --- --- --- --- --- ---
  # Don't run if the chromosomes exist
  n.exist = 0
  for(i.chr in 1:n.chr){
    # acc.s = gsub('_', '-', acc)
    acc.s = acc
    file.out = paste0(path.chr, acc.s, '_chr', i.chr, '.fasta', collapse = '')
    if(file.exists(file.out)){
      n.exist = n.exist + 1
    }
  }
  if(n.exist == n.chr){
    if(for.flag) next
    return(NULL)
  }
  #' --- --- --- --- --- --- --- --- --- --- ---
  
  pokaz('Number of existing ')
  pokaz('Accession', acc)
  q.fasta = readFastaMy(paste0(path.query, query.name$file[i.acc], collapse = ''))
  
  if(length(q.fasta) < n.chr){
    if(for.flag) next
    return(NULL)
  }
  
  if(sort.by.lengths){
    q.len = sapply(q.fasta, nchar)
    pokaz('Lengths', q.len)
    q.fasta = q.fasta[order(-q.len)]
  }
  
  pokaz('Chromosomes', names(q.fasta)[1:(n.chr)], 'will be processed')
  
  if(length(q.fasta) > n.chr){
    pokaz('Chromosomes', names(q.fasta)[(n.chr+1):length(q.fasta)], 'will NOT be processed')
  }
  
  for(i.chr in 1:n.chr){
    # acc.s = gsub('_', '-', acc)
    acc.s = acc
    file.out = paste0(path.chr, acc.s, '_chr', i.chr, '.fasta', collapse = '')
    if(file.exists(file.out)){
      if(for.flag) next
      return(NULL)
    }
    pokaz('File out', file.out)
    
    s = toupper(q.fasta[i.chr])
    writeFastaMy(s, file=file.out, append=F, seq.names = paste(acc.s, '_Chr', i.chr , sep=''))
    
    # write(paste('>', acc.s, '_Chr', i.chr , sep=''), file=file.out, append=F)
    # write(s, file=file.out, append=T)
    # write('\n', file=file.out, append=T)
  }
  
  rm(q.fasta)
}

stopCluster(myCluster)





