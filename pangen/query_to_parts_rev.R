suppressMessages({
  library("optparse")
  library(Biostrings)
  library("seqinr")       # read.fasta
  library("foreach")
  library(doParallel)
})

source("utils/utils.R")

#Rscript query_to_parts.R -n 5 -t fasta --path.chr ../pb_chromosomes/ -b 5000 --path.parts ../pb_parts/
# Rscript query_to_parts.R -n 8 -t fasta --path.chr ../ly_chromosomes/ -b 5000 --path.parts ../ly_parts/
# Rscript query_to_parts.R -n 1 -t fasta --path.chr ../rhiz_chromosomes/ -b 5000 --path.parts ../rhiz_parts/ 

myCluster <- makeCluster(30, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

# len.parts = 5000  # lengths of parts
# n.chr = 5  # number of chromosomes
# query.type = 'fasta'
# path.chr = '../pb_databases/'
# path.parts = '../pb_parts/'

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("-b", "--part.len"), type = "character", default = NULL, 
              help = "number of base pairs in the part file", metavar = "character"),
  make_option(c("-n", "--n.chr"), type = "character", default = NULL, 
              help = "number of chromosomes", metavar = "character"),
  make_option(c("-t", "--type"), type = "character", default = NULL, 
              help = "type of fasta files", metavar = "character"),
  make_option(c("-i", "--path.chr"), type = "character", default = NULL, 
              help = "pathway to the chromosome directory", metavar = "character"),
  make_option(c("-o", "--path.parts"), type = "character", default = NULL, 
              help = "pathway to the parts directory", metavar = "character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Set the number of cores for parallel processing
num_cores <- ifelse(!is.null(opt$cores), opt$cores, 30)
myCluster <- makeCluster(num_cores, type = "PSOCK")
registerDoParallel(myCluster)

# Ensure the number of chromosomes is specified and set it
n.chr <- ifelse(!is.null(opt$n.chr), as.numeric(opt$n.chr), stop("The input number of chromosomes 'n.chr' must be specified!"))

# Set chromosome and parts paths
path.chr <- ifelse(!is.null(opt$path.chr), opt$path.chr, stop("The chromosome path 'path.chr' must be specified!"))
path.parts <- ifelse(!is.null(opt$path.parts), opt$path.parts, stop("The parts path 'path.parts' must be specified!"))
if(!dir.exists(path.parts)) dir.create(path.parts)

# Common attributes
len.parts <- ifelse(!is.null(opt$part.len), as.numeric(opt$part.len), 5000)
query.type <- ifelse(!is.null(opt$type), opt$type, "fasta")


#' ----------------------------------------------------------------------
pokazStage('Chromosomes into parts')

pokaz('Directory with chromosomes:', path.chr)
files.query = list.files(path = path.chr, pattern = paste0('\\.', query.type, '$', collapse = '') )
query.name = unique(sapply(files.query, function(s) strsplit(s, '_chr')[[1]][1]))
pokaz('Names of genomes:', query.name)

if(length(query.name) == 0){
  stop('Wrong names of chromosomal files or files are not provided')
}

for.flag = F
tmp = foreach(acc = query.name, .packages=c('stringr','Biostrings', 'seqinr', 'crayon')) %dopar% {
# for.flag = T
# for(acc in query.name){
  
  for(i.chr in 1:n.chr){
    file.in = paste0(path.chr, acc, '_chr', i.chr, '.fasta', collapse = '')
    
    file.out = paste0(path.parts, acc, '_', i.chr, '.fasta', collapse = '')
    if( file.exists(file.out)) {
      if(for.flag) next
      return(NULL)
    }
    
    
    q.fasta = readFastaMy(file.in)[1]
    q.fasta = nt2seq(rev(seq2nt(q.fasta)))
    
    q.fasta = toupper(q.fasta)
    
    s = splitSeq(q.fasta, n=len.parts)
    len.chr = nchar(q.fasta)
    pos.beg = seq(1, len.chr, len.parts)
    
    if(length(s) != length(pos.beg)) stop('Problem with chunks')
    names(s) = paste('acc_', acc, '|chr_', i.chr, '|part_', 1:length(s), '|', pos.beg, sep='')
    
    file.out = paste0(path.parts, acc, '_', i.chr, '.fasta', collapse = '')
    writeFastaMy(s, file.out)
    
    rm(q.fasta)
    rm(s)
    rm(pos.beg)
  }
}






