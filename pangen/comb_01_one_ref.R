#' How the output files look like:
#'     group           name       otype   dclass      dim
#'         /           accs   H5I_GROUP                  
#'     /accs              0 H5I_DATASET    FLOAT 28940631
#'     /accs          10002 H5I_DATASET    FLOAT 28940631
#'     /accs          10015 H5I_DATASET    FLOAT 28940631

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
})

source("utils/utils.R")
source("pangen/synteny_funcs.R")

# pokazStage('Step 7. Combine reference-based alignments by chromosomes')

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.chr.len"), type="character", default=NULL, 
              help="file with lengths of chromosomes", metavar="character"),
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to consensus directory", metavar="character"),
  make_option(c("--path.aln"), type="character", default=NULL, 
              help="path to the output directory with alignments", metavar="character"),
  make_option(c("--pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("--path.chr"), type="character", default=NULL, 
              help="path to the reference file", metavar="character"),
  make_option(c("--type"), type="character", default=NULL, 
              help="type of fasta files", metavar="character")
)



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);


# ***********************************************************************
# ---- Values of parameters ----

# print(opt)

# Number of cores
num.cores.max = 10
num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))

# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) system(paste('mkdir ', path.cons, sep = ''))

# File with chromosomal lengths, which should be in the consensus dir
path.chr.len = ifelse(!is.null(opt$path.chr.len), opt$path.chr.len, 'chr_len/')
path.chr.len = paste(path.cons, path.chr.len, sep = '')
if(!dir.exists(path.chr.len)) system(paste('mkdir ', path.chr.len, sep = ''))

if (!is.null(opt$path.chr)) path.chr <- opt$path.chr  # to know the chromosomal lengths
if (!is.null(opt$pref)) base.acc.ref <- opt$pref

# Path with alignments
if (!is.null(opt$path.aln)) path.aln <- opt$path.aln

# ***********************************************************************

# ---- Get accession names ----

aln.suff <- "_full.rds"
aln.files <- list.files(path.aln)

pokaz('Paths with alignments', path.aln)
# pokaz('Files', aln.files)

aln.files <- aln.files[grep(paste0(aln.suff, "$"), aln.files)]
# pokaz(aln.files)

accessions <- sapply(aln.files, function(filename){
  parts <- unlist(strsplit(filename, "_", fixed = TRUE))
  name <- paste(parts[1:(length(parts) - 3)], collapse = "_")
  return(name)})
names(accessions) = NULL

accessions <- sort(unique(accessions))
pokaz('Accessions:', accessions)

# ---- Combinations of chromosomes query-base to create the alignments ----

chromosome.pairs <- unique(do.call(rbind, lapply(aln.files, function(filename){
  parts <- unlist(strsplit(filename, "_", fixed = TRUE))
  s.comb <- c(as.numeric(parts[length(parts) - 2]),
              as.numeric(parts[length(parts) - 1]))
  return(s.comb)})))

pokaz('Combinations:', paste(chromosome.pairs[,1], chromosome.pairs[,2], sep = '_'))

# ---- Length of reference chromosomes ----

file.chr.len = paste(path.chr.len, base.acc.ref, '_len.rds', sep = '')
# pokaz('File with chromosomal lengths', file.chr.len)
if(!file.exists(file.chr.len)){
  chr.len = c()
  
  # Define the number of chromosomes in the reference genome by the number of files in the folder, 
  # which match the reference genome name.
  pattern <- paste0("^", base.acc.ref, "_chr(\\d+)\\.fasta$")
  files.ref.chr <- list.files(path = path.chr, pattern = pattern)
  n.chr.ref = length(files.ref.chr)
  
  for(i.chr in 1:n.chr.ref){
    acc = base.acc.ref
    # print(c(i.chr, acc))
    
    # Read base chromosome
    base.file = paste0(base.acc.ref, '_chr', i.chr , '.fasta', collapse = '')
    # pokaz('Base:', base.file)
    base.fas.fw = readFastaMy(paste(path.chr, base.file, sep = ''))
    base.fas.fw = seq2nt(base.fas.fw)
    chr.len = c(chr.len, length(base.fas.fw))
  }
  saveRDS(chr.len, file.chr.len)
} else {
  chr.len = readRDS(file.chr.len)
}

# pokaz('Chromosomal lengths', chr.len)

# ----  Combine correspondence  ----

pokaz('Reference:', base.acc.ref)
max.len.gap = 20000


# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(i.chr.pair, echo = T){
  
  
  query.chr = chromosome.pairs[i.chr.pair, 1]
  base.chr = chromosome.pairs[i.chr.pair, 2]
  
  
  if(echo) pokaz('Combination', query.chr, base.chr)
  
  base.len = chr.len[base.chr]
  
  file.comb = paste(path.cons, 'comb_', query.chr, '_', base.chr,'_ref_',base.acc.ref,'.h5', sep = '')
  if (file.exists(file.comb)) file.remove(file.comb)
  h5createFile(file.comb)
  
  # Path to accessions chunks
  gr.accs <- "accs/"
  # TODO: Check the availability of the group before creating it
  h5createGroup(file.comb, gr.accs)
  
  
  # gr.break = 'break/'
  # h5createGroup(file.comb, gr.break)
  
  idx.break = 0
  # idx.gaps = rep(0, base.len)
  
  for(acc in accessions){
    
    # pokaz('Accession', acc, 'qchr', query.chr, 'bchr', base.chr)
    
    pref.comb = paste0(acc, '_', query.chr, '_', base.chr, collapse = '')
    file.aln.full <- paste(path.aln, paste0(pref.comb,  '_full.rds', collapse = ''), sep = '')
    if(!file.exists(file.aln.full)) next
    
    # Reading the alignment
    x = readRDS(file.aln.full)
    
    # Get query coordinates in base order
    x.corr = getCorresp2BaseSign(x, base.len)
    
    if(sum(duplicated(x.corr[x.corr != 0])) > 0) stop('DUPLICSTIONS')
    
    # Write into file
    suppressMessages({
      h5write(x.corr, file.comb, paste(gr.accs, '', acc, sep = ''))
    })
    
    
    rmSafe(x.corr)
    rmSafe(x)
    rmSafe(v)
    rmSafe(idx.tmp.acc)
    # rmSafe(idx.break.acc)
    
  }
  
  suppressMessages({
    # h5write(idx.break, file.comb, 'breaks_all')
    # h5write(idx.gaps, file.comb, 'gaps_all')
    h5write(base.acc.ref, file.comb, 'ref')
    
    h5write(1:base.len, file.comb, paste(gr.accs, '', base.acc.ref, sep = ''))
    # h5write(NULL, file.comb, paste(gr.break, base.acc.ref, sep = ''))
  })
  
  # rmSafe(idx.break)
  rmSafe(idx.gaps)
  
  H5close()
  gc()
}


# ***********************************************************************
# ---- Loop  ----


if(num.cores == 1){
  for(i.chr.pair in 1:nrow(chromosome.pairs)){
    loop.function(i.chr.pair)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(i.chr.pair = 1:nrow(chromosome.pairs), .packages=c('rhdf5', 'crayon'))  %dopar% { 
                              loop.function(i.chr.pair)
                            }
  stopCluster(myCluster)
}

warnings()


