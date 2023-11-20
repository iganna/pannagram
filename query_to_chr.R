suppressMessages({
  library("optparse")
  library(Biostrings)
  library("seqinr")       # read.fasta
  library("foreach")
  library(doParallel)
})

source("utils.R")

pokazStage('Genomes into chromosomes')

# Rscript query_to_chr.R -n 5 -t fasta --path.in ../pb_genomes/ --path.out ../pb_chromosomes/
# Rscript query_to_chr.R -n 8 -t fasta --path.in ../lyrata/ --path.out ../ly_chromosomes/    
# Rscript query_to_chr.R -n 1 -t fasta --path.in ../rhizobia/ --path.out ../rhiz_chromosomes/ -s T

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("-n", "--n.chr"), type = "character", default = NULL, 
              help = "number of chromosomes", metavar = "character"),
  make_option(c("-t", "--type"), type = "character", default = NULL, 
              help = "type of fasta files", metavar = "character"),
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
query.type <- ifelse(!is.null(opt$type), opt$type, ".fasta")

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

pokaz('Path with genomes:', path.query)
files.query = list.files(path = path.query, pattern = paste0('\\.', query.type, '$', collapse = '') )
query.name = gsub(paste0('*.', query.type, collapse = ''), "" ,files.query)
if(!is.null(acc.anal)){
  query.name = intersect(query.name, acc.anal)
}
if(length(query.name) == 0) stop('No accessions for the analysys.')
pokaz('Names of genomes:', query.name)

#for.flag = F
#tmp = foreach(acc = query.name, .packages=c('stringr','Biostrings', 'seqinr', 'crayon')) %dopar% {
 for.flag = T
 for(acc in query.name[1:length(query.name)]){
  
  
  #' --- --- --- --- --- --- --- --- --- --- ---
  # Don't run if the chromosomes exist
  n.exist = 0
  for(i.chr in 1:n.chr){
    acc.s = gsub('_', '-', acc)
    file.out = paste0(path.chr, acc.s, '_chr', i.chr, '.fasta', collapse = '')
    print(file.out)
    if(file.exists(file.out)){
      n.exist = n.exist + 1
    }
  }
  pokaz('Number of existing chromosomes', n.exist)
  if(n.exist == n.chr){
    if(for.flag) next
    return(NULL)
  }
  pokaz('New chromosomes will be saved')
  #' --- --- --- --- --- --- --- --- --- --- ---
  
  pokaz('Number of existing ')
  pokaz('Accession', acc)
  q.fasta = readFastaMy(paste0(path.query, acc, '.', query.type, collapse = ''))
  
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
    acc.s = gsub('_', '-', acc)
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





