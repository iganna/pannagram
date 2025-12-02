# Get syntenic positions

suppressMessages({
  library(foreach)
  library(doParallel)
  library(crayon)
  library(rhdf5)
  library(pannagram)
  library(optparse)
  library(ggplot2)
})

# ***********************************************************************
# ---- Alignment types ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

# ***********************************************************************
args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.project",      type = "character", default = NULL, help = "Path to the project"),
  make_option("--path.inter.msa",    type = "character", default = NULL, help = "Path to msa dir (internal)"),
  make_option("--path.features.msa", type = "character", default = NULL, help = "Path to msa dir (features)"),
  make_option("--path.figures",      type = "character", default = "",   help = "Path to folder with figures"),
  make_option("--ref",               type = "character", default = "",   help = "Prefix of the reference file"),
  make_option("--cores",             type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--aln.type",          type = "character", default = aln.type.msa, help = "Prefix for the output file"),
  make_option("--wnd.size",          type = "integer", default = 100000, help = "Window size (plotting)")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ***********************************************************************
# ---- Variables ----

num.cores <- opt$cores
wnd.size <- opt$wnd.size
block.len.min = 20000

# ***********************************************************************
# ---- Paths ----

path.project <- opt$path.project
if(!dir.exists(path.project)) stop('Consensus folder doesn’t exist')

path.features.msa <- opt$path.features.msa
if(!dir.exists(path.features.msa)) stop('features/msa dir doesn’t exist')

path.inter.msa <- opt$path.inter.msa
if(!dir.exists(path.inter.msa)) stop('internal/msa dir doesn’t exist')

path.figures <- opt$path.figures
if(!dir.exists(path.figures)) stop('Consensus folder doesn’t exist')

# ***********************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

# Alignment prefix
if (!is.null(opt$aln.type)) {
  aln.type = opt$aln.type
} else {
  aln.type = aln.type.msa
}

# Reference genome
ref.name <- opt$ref
if(ref.name == "NULL" || is.null(ref.name)) ref.name <- ''

# Common code for aln.pref, ref.suffix and s.combinations
source(system.file("utils/chunk_combinations.R", package = "pannagram")) 
source(system.file("visualisation/panplot.R", package = "pannagram")) 

# ***********************************************************************
# ---- MAIN program body ----

file.blocks = paste0(path.inter.msa, 'syn_blocks_', aln.type, ref.suff,'.rds')
# pokaz(file.blocks)

if(!file.exists(file.blocks)){
  df.all = c()
  for(s.comb in s.combinations){
    
    pokaz('Combination', s.comb)
    # --- --- --- --- --- --- --- --- --- --- ---
    
    file.comb.in = paste0(path.features.msa, aln.pref, s.comb, ref.suff,'.h5')
    
    groups = h5ls(file.comb.in)
    accessions = groups$name[groups$group == gr.accs.b]
    
    processAcc <- function(acc) {
      pokaz('Accession', acc)
      v <- h5read(file.comb.in, paste0(gr.accs.e, acc))
      
      # save(list = ls(), file = "tmp_workspace_acc.RData")
      
      df.acc <- getBlocks(v, f.split = FALSE)
      if(nrow(df.acc) == 0) {
        return(NULL) 
      }
      df.acc$acc <- acc
      df.acc$comb <- s.comb
      
      return(df.acc)
    }
    
    if (num.cores == 1) {
      df <- do.call(rbind, lapply(accessions, processAcc))
    } else {
      myCluster <- makeCluster(num.cores, type = "PSOCK") 
      registerDoParallel(myCluster) 
      
      df <- foreach(acc = accessions, .packages = c('rhdf5', 'crayon', 'pannagram'), .combine = rbind) %dopar% {
        processAcc(acc)
      }
      stopCluster(myCluster)
    }
    
    df.all = rbind(df.all, df)
    
    gc()
    
  }
  
  idx.dir = which(df.all$own.b > df.all$own.e)
  tmp = df.all$own.b[idx.dir]
  df.all$own.b[idx.dir] = df.all$own.e[idx.dir]
  df.all$own.e[idx.dir] = tmp
  
  tmp = df.all$pan.b[idx.dir]
  df.all$pan.b[idx.dir] = df.all$pan.e[idx.dir]
  df.all$pan.e[idx.dir] = tmp
  
  saveRDS(df.all, file.blocks)  
} else {
  pokaz("File with blocks was generated in advance", file.blocks)
  df.all = readRDS(file.blocks)
}


# ***********************************************************************
# ---- Plot  ----

pokaz("Plots...")
accessions = unique(df.all$acc)

for(s.comb in unique(df.all$comb)){
  pokaz('Combination', s.comb)
  
  i.chr = parseStrings(s.comb, n = 1, split = '_', numeric = T)
  p = panplot(path.project = path.project, i.chr = i.chr, aln.type = aln.type, ref.acc = ref.name)
  
  savePDF(p, path = path.figures, name = paste0('synteny_', aln.pref, s.comb, ref.suff), width = 6, height = length(accessions) / 6 + 1)
}


warnings()

