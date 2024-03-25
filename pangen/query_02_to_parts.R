suppressMessages({
  library("optparse")
  # library(Biostrings)
  # library("seqinr")       # read.fasta
  library("foreach")
  library(doParallel)
  # library(stringi)
})

source("utils/utils.R")

pokazStage('Step 2. Chromosomes into parts')

# ***********************************************************************
# ---- Command line arguments ----
args = commandArgs(trailingOnly=TRUE)


option_list <- list(
  make_option(c("--all.chr"), type = "character", default = NULL, 
              help = "Flag to use all chromosomes, not only the provided number", metavar = "character"),
  make_option(c("--n.chr"), type = "character", default = NULL, 
              help = "number of chromosomes", metavar = "character"),
  make_option(c("--part.len"), type = "character", default = NULL, 
              help = "number of base pairs in the part file", metavar = "character"),
  make_option(c("--path.chr"), type = "character", default = NULL, 
              help = "pathway to the chromosome directory", metavar = "character"),
  make_option(c("--path.parts"), type = "character", default = NULL, 
              help = "pathway to the parts directory", metavar = "character"),
  make_option(c("--filter_rep"), type = "character", default = NULL, 
              help = "flag to keep or not repeats", metavar = "character"),
  make_option(c("--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

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

# Set chromosome and parts paths
path.chr <- ifelse(!is.null(opt$path.chr), opt$path.chr, stop("The chromosome path 'path.chr' must be specified!"))
path.parts <- ifelse(!is.null(opt$path.parts), opt$path.parts, stop("The parts path 'path.parts' must be specified!"))
if(!dir.exists(path.parts)) dir.create(path.parts)

# Common attributes
len.parts <- ifelse(!is.null(opt$part.len), as.numeric(opt$part.len), 5000)
filter_rep <- as.numeric(ifelse(!is.null(opt$filter_rep), as.numeric(opt$filter_rep), 0))

# ***********************************************************************
# ---- Preparation ----

pokaz('Directory with chromosomes:', path.chr)
files.query = list.files(path = path.chr, pattern = paste0('\\.', 'fasta', '$', collapse = '') )
files.query <- sub("\\.fasta$", "", files.query)
query.name = unique(sapply(files.query, function(s) strsplit(s, '_chr')[[1]][1]))
pokaz('Names of genomes:', query.name)

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
loop.function <- function(i.comb, echo = T){
  
  acc <- combinations$acc[i.comb]
  i.chr <- combinations$i.chr[i.comb]
  
  file.in = paste0(path.chr, acc, '_chr', i.chr, '.fasta', collapse = '')
  
  file.out = paste0(path.parts, acc, '_chr', i.chr, '.fasta', collapse = '')
  if( file.exists(file.out)) {
    return(NULL)
  }
  
  q.fasta = readFastaMy(file.in)[1]
  q.fasta = toupper(q.fasta)
  
  s = splitSeq(q.fasta, n=len.parts)
  len.chr = nchar(q.fasta)
  pos.beg = seq(1, len.chr, len.parts)
  
  if(length(s) != length(pos.beg)) stop('Problem with chunks')
  names(s) = paste('acc_', acc, '|chr_', i.chr, '|part_', 1:length(s), '|', pos.beg, sep='')
  
  
  if(filter_rep == 0){
    file.out = paste0(path.parts, acc, '_chr', i.chr, '.fasta', collapse = '')
    writeFastaMy(s, file.out)
  } else {
    
    file.out = paste0(path.parts, acc, '_chr', i.chr, '.fasta', collapse = '')
    
    file.out.rest = paste0(path.parts, acc, '_chr', i.chr, '.rest', collapse = '')
    
    seqs.score = sapply(s, repeatScore)
    
    writeFastaMy(s[seqs.score <= 0.2], file.out)
    writeFastaMy(s[seqs.score > 0.2], file.out.rest)
    
  }
  
  rmSafe(q.fasta)
  rmSafe(s)
  rmSafe(pos.beg)
  rmSafe(seqs.score)
}

# ***********************************************************************
# ---- Loop  ----

if(num.cores == 1){
  for(i.comb in 1:nrow(combinations)){
    loop.function(i.comb)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK")
  registerDoParallel(myCluster)
  
  tmp.output = foreach(i.comb = 1:nrow(combinations), 
                       .packages=c(
                                   # 'stringr',
                                   # 'Biostrings', 
                                   # 'seqinr', 
                                   # 'stringi',
                                  'crayon'
                                   )) %dopar% {
    loop.function(i.comb)
  }
  stopCluster(myCluster)
}

warnings()

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










