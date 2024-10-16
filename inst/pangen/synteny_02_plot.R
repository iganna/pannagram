### imports
source(system.file("visualisation/visualisation.R", package = "pannagram"))
source(system.file("pangen/synteny_func_plot.R", package = "pannagram"))
source(system.file("utils/utils.R", package = "pannagram"))

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
  
  make_option(c("--path.plot"),   type = "character", default = NULL, help = "Path to the output directory"),
  make_option(c("--path.chr"),    type = "character", default = NULL, help = "Path to chromosomes"),
  make_option(c("--path.aln"),    type = "character", default = NULL, help = "Path to the alignment directory"),
  
  make_option(c("--ref"),           type = "character", default = NULL, help = "Name of the reference genome"),
  make_option(c("--accessions"),    type = "character", default = NULL, help = "File containing accessions to analyze"),
  
  make_option(c("--path.log"),    type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"),   type = "character", default = NULL, help = "Level of log to be shown on the screen"),
  make_option(c("--cores"),       type = "integer",   default = 1,    help = "Number of cores to use")
)



opt <- parse_args(OptionParser(option_list=option_list))

# Logging
source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

if (
    is.null(opt$ref) ||
    is.null(opt$path.chr) ||
    is.null(opt$path.plot) ||
    is.null(opt$path.aln)
) stop("All mandatory arguments must be provided.")

path.chr <- opt$path.chr
path.plot <- opt$path.plot
path.aln <- opt$path.aln
ref <- opt$ref
num.cores <- opt$cores

pokaz('Reference genome:', ref, 
      file=file.log.main, echo=echo.main)


# ***********************************************************************
# ---- Accessions ----

file.acc <- ifelse(!is.null(opt$accessions), opt$accessions, stop("File with accessions are not specified"))
tmp <- read.table(file.acc, stringsAsFactors = F)
accessions <- as.character(tmp[,1])
accessions = setdiff(accessions, ref)
pokaz('Names of genomes for the analysis:', accessions, 
      file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(acc, 
                          echo.loop=T, 
                          file.log.loop=NULL){
  
  # Log files
  file.log.loop = paste0(path.log, 'loop_file_', acc, '.log')
  if(!file.exists(file.log.loop)){
    invisible(file.create(file.log.loop))
  }
  
  # ---- Check log Done ----
  if(checkDone(file.log.loop)){
    return(NULL)
  }
  
  pokaz('Accession', acc, file=file.log.loop, echo=echo.loop)
  
  # Lengths of chromosomes for the accession and reference
  chr.len = c()
  for(id in c(ref, acc)){
    file.chr.len = paste0(path.chr, id, '_chr_len.txt')  
    # pokaz(file.chr.len)
    if(file.exists(file.chr.len)){
      chr.len.id = read.table(file.chr.len, header=T)
      chr.len = rbind(chr.len, chr.len.id)  
    } else {
      for(i.chr in 1:1e10){
        file.chr =  paste0(path.chr, id, '_chr',i.chr,'.fasta')
        # pokaz('Chromosomel file', file.chr)
        if(!file.exists(file.chr)) break
        seq.chr = readFastaMy(file.chr)
        chr.len.tmp = data.frame(acc = id,
                                 chr = i.chr,
                                 len = nchar(seq.chr))
        chr.len = rbind(chr.len, chr.len.tmp)
      }
      write.table(chr.len, file.chr.len, sep = '\t', col.names = T, row.names = F, quote = F)
    }
  }
  
  # ---- Testing ----
  # file.ws = "tmp_workspace.RData"
  # all.local.objects <- ls()
  # save(list = all.local.objects, file = file.ws)
  # pokaz('Workspace is saved in', file.ws, file=file.log.loop, echo=echo.loop)
  # stop('Enough..')
  
  # ---- Plot ----
  # Get ggplot with the synteny
  p <- plotSynAllChr(path.aln,
                     acc=acc,
                     ref=ref,
                     chr.len=chr.len)
  
  # Save
  pdf.name <- paste0(ref, "-", id)
  savePDF(p, path = path.plot, name = pdf.name)
  
}

if(num.cores == 1){
  for(acc in accessions){
    loop.function(acc,
                  echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(acc = accessions, 
                .packages=c('crayon', 'ggplot2'), 
                .verbose = F)  %dopar% { 
                  loop.function(acc,
                                echo.loop=echo.loop)
                }
  stopCluster(myCluster)
}


pokaz('Done.',
      file=file.log.main, echo=echo.main)

