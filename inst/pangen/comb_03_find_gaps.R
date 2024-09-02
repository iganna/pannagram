# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
})

source(system.file("utils/utils.R", package = "pannagram"))
# source("pangen/synteny_funcs.R")

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to consensus directory", metavar="character"),
  make_option(c("-p", "--ref.pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("--path.log"), type = "character", default = NULL,
              help = "Path for log files", metavar = "character"),
  make_option(c("--log.level"), type = "character", default = NULL,
              help = "Level of log to be shown on the screen", metavar = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ***********************************************************************
# ---- Values of parameters ----

# Set the number of cores for parallel processing
num.cores.max = 10
num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))


# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop('Consensus folder doesn’t exist')

# Reference genome
if (is.null(opt$ref.pref)) {
  stop("Error: ref.pref is NULL")
} else {
  ref.pref <- opt$ref.pref
}

# ***********************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

# Testing
if(F){
  library(rhdf5)
  path.cons = './'
  ref.pref = '0'
  options("width"=200, digits=10)
}


s.pattern <- paste0("^", 'res_', ".*", '_ref_', ref.pref)
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub("res_", "", files)
pref.combinations <- sub("_ref.*$", "", pref.combinations)

pokaz('Reference:', ref.pref, file=file.log.main, echo=echo.main)
pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ----  Combine correspondence  ----

gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'
max.len.gap = 20000

gr.blocks = 'blocks/'


# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(s.comb,
                          echo.loop=T, 
                          file.log.loop=NULL){
  
  # Log files
  if (is.null(file.log.loop)){
    file.log.loop = paste0(path.log, 'loop_file_', 
                           s.comb,
                           '.log')
    invisible(file.create(file.log.loop))
  }
  
  # --- --- --- --- --- --- --- --- --- --- ---
  
  file.comb = paste0(path.cons, 'res_', s.comb,'_ref_',ref.pref,'.h5')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  
  # Create group for blocks
  suppressMessages({ h5createGroup(file.comb, gr.blocks) })
  
  idx.break = c()
  for(acc in accessions){
    
    pokaz('Accession', acc, 'combination', s.comb, file=file.log.loop, echo=echo.loop)
    
    v.init = h5read(file.comb, paste0(gr.accs.e, acc))
    if(length(v.init) == 0) {
      stop(sprintf("Skipping empty accession: %s\n", acc))
    }
    v = v.init
    
    
    # ----  Find breaks  ----
    
    # Find blocks of additional breaks
    v = cbind(v, 1:length(v))                       # column 2 - in ref-based coordinates
    
    v = v[!is.na(v[,1]),]                           # column 1 - existing coordinates of accessions
    v = v[v[,1] != 0,]                              # column 1 - existing coordinates of accessions
    
    if(nrow(v) == 0) {
      stop(sprintf("Skipping empty v for accession: %s\n", acc))
    }
    v = cbind(v, 1:nrow(v))                       # 3 - ranked order in ref-based coordinates
    v = cbind(v, rank(abs(v[,1])) * sign(v[,1]))  # 4 - signed-ranked-order in accessions coordinates 
    
    # Save blocks
    idx.block.tmp = which(abs(diff(v[,4])) != 1)
    idx.block.df = data.frame(beg = v[c(1,idx.block.tmp+1), 2], end = v[c(idx.block.tmp, nrow(v)), 2])

    v.block = rep(0, length(v.init))
    for(i.bl in 1:nrow(idx.block.df)){
      v.block[idx.block.df$beg[i.bl]:idx.block.df$end[i.bl]] = i.bl
    }
    
    suppressMessages({
      h5write(v.block, file.comb, paste0(gr.blocks, acc))
    })

    # v = v[order(v[,1]),]  # not necessary
    
    # with the absence, but neighbouring
    idx.tmp = which( (abs(diff(v[,4])) == 1) &  # Neighbouring in accession-based order
                       (abs(diff(abs(v[,3])) == 1)) &  # Neighbouring in ref-based order
                       (abs(diff(v[,1])) <= max.len.gap) &  # Filtering by length in accession coordinates
                       (abs(diff(v[,2])) <= max.len.gap) &  # Filtering by length in reference coordinates
                       (abs(diff(v[,1])) > 1))  # NOT neighbouring in accession-specific coordinates
    
    # Fix (beg < end) order
    if(length(idx.tmp) == 0) {
      stop(sprintf("Skipping empty idx.tmp for accession: %s\n", acc))
    }
    idx.tmp.acc = data.frame(beg = v[idx.tmp,2], end = v[idx.tmp+1,2], acc = acc)
    idx.ord = which(idx.tmp.acc$beg > idx.tmp.acc$end)
    if(length(idx.ord) > 0){
      tmp = idx.tmp.acc$beg[idx.ord]
      idx.tmp.acc$beg[idx.ord] = idx.tmp.acc$end[idx.ord]
      idx.tmp.acc$end[idx.ord] = tmp
    }
    # idx.tmp.acc = idx.tmp.acc[order(idx.tmp.acc$beg),]  # order ONLY if ordered before
    
    # Remove overlaps
    idx.overlap = which( (idx.tmp.acc$beg[-1] - idx.tmp.acc$end[-nrow(idx.tmp.acc)]) <= 3)
    
    i.cnt = 0
    if(length(idx.overlap) > 0){
      j.ov = 0
      for(i.ov in idx.overlap){
        if(i.ov <= j.ov) next
        j.ov = i.ov + 1
        while(j.ov %in% idx.overlap){
          j.ov = j.ov + 1
        }
        # print(c(i.ov, j.ov))
        i.cnt = i.cnt + 1
        idx.tmp.acc$end[i.ov] = idx.tmp.acc$end[j.ov]
      }
      idx.tmp.acc = idx.tmp.acc[-(idx.overlap+1),]
    }
    
    
    # Save breaks
    idx.break = rbind(idx.break, idx.tmp.acc)
    
    rmSafe(x.corr)
    rmSafe(x)
    rmSafe(v)
    rmSafe(v.init)
    rmSafe(idx.tmp.acc)
    rmSafe(idx.break.acc)
    
  }
  
  file.breaks = paste0(path.cons, 'breaks_', s.comb,'_ref_',ref.pref,'.rds')
  saveRDS(idx.break, file.breaks)
  
  rmSafe(idx.break)
  
  H5close()
  gc()
  
  pokaz('Done.', file=file.log.loop, echo=echo.loop)
  return(NULL)
}


# ***********************************************************************
# ---- Loop  ----

if(num.cores == 1){
  file.log.loop = paste0(path.log, 'loop_all.log')
  invisible(file.create(file.log.loop))
  for(s.comb in pref.combinations){
    loop.function(s.comb,
                  file.log.loop = file.log.loop, 
                  echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(s.comb = pref.combinations, 
                .packages=c('rhdf5', 'crayon'))  %dopar% { 
                  loop.function(s.comb,echo.loop=echo.loop)
                }
  stopCluster(myCluster)
}

warnings()

pokaz('Done.',
      file=file.log.main, echo=echo.main)
