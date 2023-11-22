#' How the output files look like:
#'     group           name       otype   dclass      dim
#'         /           accs   H5I_GROUP                  
#'     /accs          acc_0 H5I_DATASET    FLOAT 28940631
#'     /accs      acc_10002 H5I_DATASET    FLOAT 28940631
#'     /accs      acc_10015 H5I_DATASET    FLOAT 28940631


suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
})

source("utils.R")
source("synteny_funcs.R")

pokazStage('Combine: alignments by chromosomes')

args = commandArgs(trailingOnly=TRUE)


option_list = list(
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to consensus directory", metavar="character"),
  make_option(c("-r", "--ref0"), type="character", default=NULL, 
              help="reference file 1", metavar="character"),
  make_option(c("-p", "--ref1"), type="character", default=NULL, 
              help="reference file 2", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# Set the number of cores for parallel processing
num.cores.max = 10
num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))
myCluster <- makeCluster(num.cores, type = "PSOCK")
registerDoParallel(myCluster)

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



# ---- Combinations of chromosomes query-base to create the alignments ----


# Find all fines with the prefix of both references and common suffixes
# path.cons = './'
# ref0 = '0'
# ref1 = '6046-v1.1'

files.pref <- lapply(c(ref0, ref1), function(ref) {
  suff <- paste0('_ref_', ref)
  files <- list.files(path = path.cons, pattern = suff, full.names = FALSE)
  sapply(files, function(s) strsplit(s, suff)[[1]][1])
})

pref.combinations <- Reduce(intersect, files.pref)
pref.combinations <- gsub("comb_", "", pref.combinations)
pokaz('Combinations', pref.combinations)



# ----  Combine correspondence  ----

pokaz('References:', ref0, ref1)


gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'
max.len.gap = 20000

#flag.for = F
#tmp = foreach(s.comb = pref.combinations, .packages=c('rhdf5', 'crayon'))  %dopar% {  # which accession to use
flag.for = T
for(s.comb in pref.combinations){

  file.comb0 = paste(path.cons, 'comb_',s.comb,'_ref_',ref0,'.h5', sep = '')
  file.comb1 = paste(path.cons, 'comb_',s.comb,'_ref_',ref1,'.h5', sep = '')
  
  # Combined file. If it exists, then use it for the growing correspondence
  file.res = paste(path.cons, 'res_',s.comb,'_ref_',ref0,'.h5', sep = '')
  if(file.exists(file.res)){
    file.comb0 = file.res
  } else {
    h5createFile(file.res)
    h5createGroup(file.res, gr.accs.e)
    h5createGroup(file.res, gr.break.e)
  }
  
  # Get the corresponsing function between two references
  s = paste(gr.accs.e, 'acc_', ref1, sep = '')
  
  pokaz('File0', file.comb0)
  print(h5ls(file.comb0))
  
  
  f01 <- h5read(file.comb0, s)
  base.len = length(f01)
  idx01 = which(f01 != 0)
  f01 = f01[idx01]
  
  groups0 = h5ls(file.comb0)
  groups1 = h5ls(file.comb1)
  
  accessions = intersect(groups0$name[groups0$group == gr.accs.b], 
                         groups1$name[groups1$group == gr.accs.b])  # full name of accessions
  # accessions <- sub("^.*_", "", accessions)
  
  for(acc in accessions){
    # if(acc == ref0) next
    # if(acc == ref1) next
    
    print(acc)
    s = paste('/',gr.accs.e, acc, sep = '')
    
    # Data from the main reference
    v0 = h5read(file.comb0, s)
    pokaz('Vector in the ref0 file', length(v0))
    v0 = v0[idx01]
    pokaz('Vector of meaningfull positions', length(v0))
    
    # Data from the second reference
    v1 = h5read(file.comb1, s)
    pokaz('Vector in the ref1 file', length(v1))
    pokaz('Length of function', length(f01))
    v01 = v1[abs(f01)] * sign(f01)
    v01[v0 != v01] = 0
    
    pokaz('Length of resultant correspondence', length(v01))
    pokaz('Sum of matches', sum(v01 != 0))
    
    # Turn into real coorsinates back
    v.final = rep(0, base.len)
    v.final[idx01] = v01
    
    pokaz('Length of saved vector', length(v.final))
    
    suppressMessages({
      h5write(v.final, file.res, s)
    })
  }
  
  H5close()
  gc()
  
}





