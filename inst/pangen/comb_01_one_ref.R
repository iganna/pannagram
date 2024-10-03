#' How the output files look like:
#'     group           name       otype   dclass      dim
#'         /           accs   H5I_GROUP                  
#'     /accs              0 H5I_DATASET    FLOAT 28940631
#'     /accs          10002 H5I_DATASET    FLOAT 28940631
#'     /accs          10015 H5I_DATASET    FLOAT 28940631
#'     
#' In combo-files the reference accession column will not have zeros, because it will participate further lika a function.

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/synteny_func.R", package = "pannagram"))

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
              help="type of fasta files", metavar="character"),
  make_option(c("--path.log"), type = "character", default = NULL,
              help = "Path for log files", metavar = "character"),
  make_option(c("--log.level"), type = "character", default = NULL,
              help = "Level of log to be shown on the screen", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

# ***********************************************************************
# ---- Values of parameters ----

# print(opt)

# Number of cores
num.cores.max = 10
num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))

# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) system(paste0('mkdir ', path.cons))

# File with chromosomal lengths, which should be in the consensus dir
path.chr.len = ifelse(!is.null(opt$path.chr.len), opt$path.chr.len, 'chr_len/')
path.chr.len = paste0(path.cons, path.chr.len)
if(!dir.exists(path.chr.len)) system(paste0('mkdir ', path.chr.len))

if (!is.null(opt$path.chr)) path.chr <- opt$path.chr  # to know the chromosomal lengths
if (!is.null(opt$pref)) base.acc.ref <- opt$pref

# Path with alignments
if (!is.null(opt$path.aln)) path.aln <- opt$path.aln

# ***********************************************************************

# ---- Get accession names ----

aln.suff <- "_full.rds"
aln.files <- list.files(path.aln)

pokaz('Paths with alignments', path.aln, file=file.log.main, echo=echo.main)
# pokaz('Files', aln.files, file=file.log.main, echo=echo.main)

aln.files <- aln.files[grep(paste0(aln.suff, "$"), aln.files)]
# pokaz(aln.files)

accessions <- sapply(aln.files, function(filename){
  parts <- unlist(strsplit(filename, "_", fixed = TRUE))
  name <- paste(parts[1:(length(parts) - 3)], collapse = "_")
  return(name)})
names(accessions) = NULL

accessions <- sort(unique(accessions))
pokaz('Accessions:', accessions, file=file.log.main, echo=echo.main)

# ---- Combinations of chromosomes query-base to create the alignments ----

chromosome.pairs <- unique(do.call(rbind, lapply(aln.files, function(filename){
  parts <- unlist(strsplit(filename, "_", fixed = TRUE))
  s.comb <- c(as.numeric(parts[length(parts) - 2]),
              as.numeric(parts[length(parts) - 1]))
  return(s.comb)})))

pokaz('Combinations:', paste(chromosome.pairs[,1], chromosome.pairs[,2], sep = '_'), file=file.log.main, echo=echo.main)

# ---- Length of reference chromosomes ----

file.chr.len = paste0(path.chr.len, base.acc.ref, '_len.rds')
pokaz('File with chromosomal lengths', file.chr.len, file=file.log.main, echo=echo.main)
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
    base.fas.fw = readFastaMy(paste0(path.chr, base.file))
    base.fas.fw = seq2nt(base.fas.fw)
    chr.len = c(chr.len, length(base.fas.fw))
  }
  saveRDS(chr.len, file.chr.len)
} else {
  chr.len = readRDS(file.chr.len)
}

pokaz('Chromosomal lengths:', chr.len, file=file.log.main, echo=echo.main)

# ----  Combine correspondence  ----

pokaz('Reference:', base.acc.ref, file=file.log.main, echo=echo.main)


# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(i.chr.pair, 
                          echo.loop=T, 
                          file.log.loop=NULL){
  
  
  query.chr = chromosome.pairs[i.chr.pair, 1]
  base.chr = chromosome.pairs[i.chr.pair, 2]
  
  # Log files
  if (is.null(file.log.loop)){
    file.log.loop = paste0(path.log, 'loop_file_', 
                           query.chr, '_', base.chr,
                           '.log')
    invisible(file.create(file.log.loop))
  }
  
  # ---- Check log Done ----
  if(checkDone(file.log.loop)){
    return()
  }
  
  pokaz('Combination', query.chr, base.chr, file=file.log.loop, echo=echo.loop)
  pokaz('Chromosomal length', chr.len, file=file.log.loop, echo=echo.loop)
  base.len = chr.len[base.chr]
  
  file.comb = paste0(path.cons, aln.type.ref, query.chr, '_', base.chr,'_ref_',base.acc.ref,'.h5')
  if (file.exists(file.comb)) file.remove(file.comb)
  h5createFile(file.comb)
  
  # Path to accessions chunks
  # TODO: Check the availability of the group before creating it
  h5createGroup(file.comb, gr.accs.e)
  
  
  # gr.break = 'break/'
  # h5createGroup(file.comb, gr.break)
  
  idx.break = 0
  # idx.gaps = rep(0, base.len)
  
  for(acc in accessions){
    
    pokaz('Accession', acc, 'qchr', query.chr, 'bchr', base.chr, 
                   file=file.log.loop, echo=echo.loop)
    
    pref.comb = paste0(acc, '_', query.chr, '_', base.chr, collapse = '')
    file.aln.full <- paste(path.aln, paste0(pref.comb,  '_full.rds', collapse = ''), sep = '')
    if(!file.exists(file.aln.full)) next
    
    pokaz('Alignment file:', file.aln.full, file=file.log.loop, echo=echo.loop)
    
    # Reading the alignment
    x = readRDS(file.aln.full)
    
    pokaz('Base len', base.len, file=file.log.loop, echo=echo.loop)
    # saveRDS(x, 'tmp.rds')
    
    # Get query coordinates in base order
    x.corr = getCorresp2BaseSign(x, base.len)
    
    if(sum(duplicated(x.corr[x.corr != 0])) > 0) stop('DUPLICSTIONS', sum(duplicated(x.corr[x.corr != 0])))
    
    # Write into file
    suppressMessages({
      h5write(x.corr, file.comb, paste0(gr.accs.e, '', acc))
    })
    
    
    rmSafe(x.corr)
    rmSafe(x)
    rmSafe(v)
    rmSafe(idx.tmp.acc)
    # rmSafe(idx.break.acc)
    
  }
  
  suppressMessages({
    
    h5write(base.acc.ref, file.comb, v.ref.name)
    h5write(base.len, file.comb, v.len)
    
    h5write(1:base.len, file.comb, paste0(gr.accs.e, '', base.acc.ref))
  
  })
  
  # rmSafe(idx.break)
  rmSafe(idx.gaps)
  
  H5close()
  gc()
  
  pokaz('Done.', file=file.log.loop, echo=echo.loop)
  return(NULL)
}


# ***********************************************************************
# ---- Loop  ----


if(num.cores == 1){
  # file.log.loop = paste0(path.log, 'loop_all.log')
  # invisible(file.create(file.log.loop))
  for(i.chr.pair in 1:nrow(chromosome.pairs)){
    loop.function(i.chr.pair,
                  # file.log.loop = file.log.loop, 
                  echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(i.chr.pair = 1:nrow(chromosome.pairs), 
                .packages=c('rhdf5', 'crayon'))  %dopar% { 
                  loop.function(i.chr.pair,
                                echo.loop=echo.loop)
                }
  stopCluster(myCluster)
}

warnings()

pokaz('Done.',
      file=file.log.main, echo=echo.main)

