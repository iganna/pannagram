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

source(system.file("utils/utils.R", package = "pannagram"))
# source("pangen/synteny_funcs.R")



# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to consensus directory", metavar="character"),
  make_option(c("-r", "--ref0"), type="character", default=NULL, 
              help="reference file 1", metavar="character"),
  make_option(c("-p", "--ref1"), type="character", default=NULL, 
              help="reference file 2", metavar="character"),
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

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores.max = 10
num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))

# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop('Consensus folder doesn’t exist')

# Reference genomes
if (is.null(opt$ref0)) {
  stop("opt$pref is NULL")
} else {
  ref0 <- opt$ref0
}
if (is.null(opt$ref1)) {
  stop("opt$pref is NULL")
} else {
  ref1 <- opt$ref1
}

# ***********************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

pattern <- paste0("^",aln.type.ref,"[0-9]+_[0-9]+_ref_.*\\.h5$")
combo_files <- list.files(path = path.cons, pattern = pattern, full.names = F)

extract_xy <- function(filename) {
  parts <- strsplit(filename, "_")[[1]]
  x <- parts[2]
  y <- parts[3]
  return(paste(x, y, sep = '_'))
}
pref.combinations <- unique(sapply(combo_files, extract_xy))

if(length(pref.combinations) == 0) {
  stop('No files with the ref-based alignments are found')
}

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ----  Combine correspondence  ----

pokaz('References:', ref0, ref1, file=file.log.main, echo=echo.main)

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
  
  # ---- Check log Done ----
  if(checkDone(file.log.loop)){
    return()
  }
  
  # --- --- --- --- --- --- --- --- --- --- ---
  file.comb0 = paste0(path.cons, aln.type.ref ,s.comb,'_ref_',ref0,'.h5')
  file.comb1 = paste0(path.cons, aln.type.ref ,s.comb,'_ref_',ref1,'.h5')
  
  # Combined file. If it exists, then use it for the growing correspondence
  file.res = paste0(path.cons, aln.type.comb, s.comb, '.h5')
  if(file.exists(file.res)){
    file.comb0 = file.res
    
    # idx.trust
    idx.trust = h5read(file.comb0, v.idx.trust)
  } else {
    h5createFile(file.res)
    h5createGroup(file.res, gr.accs.e)
    
    idx.trust = 0
  }
  
  # Get the corresponsing function between two references
  
  s.ref1 = paste0(gr.accs.e, '', ref1)
  f01 <- h5read(file.comb0, s.ref1)
  idx01.b = f01 != 0
  idx01 = which(idx01.b)  # idx which we trust
  f01 = f01[idx01]
  
  # Idx trust
  suppressMessages({
    idx.trust = idx.trust + idx01.b * 1
    h5write(idx.trust, file.res, v.idx.trust)
  })
  
  # Get accessions to combine
  groups0 = h5ls(file.comb0)
  groups1 = h5ls(file.comb1)
  
  accessions = intersect(groups0$name[groups0$group == gr.accs.b], 
                         groups1$name[groups1$group == gr.accs.b])  # full name of accessions
  # accessions <- sub("^.*_", "", accessions)
  
  for(acc in accessions){
    # if(acc == ref0) next
    # if(acc == ref1) next
    
    pokaz('Accession', acc, file=file.log.loop, echo=echo.loop)
    s = paste0('/',gr.accs.e, acc)
    
    # Data from the main reference
    v0 = h5read(file.comb0, s)
    # pokaz('Vector in the ref0 file', length(v0), file=file.log.loop, echo=echo.loop)
    v.final = v0
    v.final[idx01] = 0  # This is a potential problem... after all rounds of comb, we should leave only those positions, which were checked at least twice.
    v0 = v0[idx01]
    pokaz('Vector of meaningfull positions', length(v0), file=file.log.loop, echo=echo.loop)
    
    
    # Data from the second reference
    v1 = h5read(file.comb1, s)
    v.final[v.final %in% v1] = 0
    pokaz('Vector in the ref1 file', length(v1), file=file.log.loop, echo=echo.loop)
    pokaz('Length of function', length(f01), file=file.log.loop, echo=echo.loop)
    v01 = v1[abs(f01)] * sign(f01)
    
    #v01[(v0 != v01) & (v0 != 0)] = 0
    # idx.lost = (v0 != 0) & (v01 == 0) & !(v0 %in% v01) 
    # v01[idx.lost] = v0[idx.lost]
    # 
    
    v0[(v0 != v01) & (v01 != 0)] = 0
    
    pokaz('Length of resultant correspondence', length(v01), file=file.log.loop, echo=echo.loop)
    pokaz('Sum of matches', sum(v01 != 0), file=file.log.loop, echo=echo.loop)
    
    v.final[idx01] = v0
    
    dup.value = setdiff(unique(v.final[duplicated(v.final)]), 0)
    if(length(dup.value > 0)){
      v.final[v.final %in% dup.value] == 0
      pokaz('Number of duplicated', length(dup.value), file=file.log.loop, echo=echo.loop)
    }
    
    pokaz('Length of saved vector', length(v.final), file=file.log.loop, echo=echo.loop)
    
    suppressMessages({
      h5write(v.final, file.res, s)
    })
  }
  
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
  for(s.comb in pref.combinations){
    loop.function(s.comb,
                  # file.log.loop = file.log.loop, 
                  echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(s.comb = pref.combinations, 
                .packages=c('rhdf5', 'crayon'))  %dopar% { 
                  loop.function(s.comb,
                                echo.loop=echo.loop)
                }
  stopCluster(myCluster)
}

warnings()

pokaz('Done.',
      file=file.log.main, echo=echo.main)


# ***********************************************************************
# ---- Manual testing ----

if(F){
source(system.file("utils/utils.R", package = "pannagram"))
  path.cons = './'
  ref0 = '0'
  ref1 = '6046-v1.1'
  library(rhdf5)
}



