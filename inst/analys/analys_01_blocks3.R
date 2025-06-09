# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
  library(pannagram)
})

source(system.file("utils/chunk_hdf5.R", package = "pannagram"))


args = commandArgs(trailingOnly=TRUE)

option_list <- list(
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

source(system.file("utils/chunk_logging.R", package = "pannagram"))

aln.type <- opt$aln.type
ref.name <- opt$ref
num.cores <- opt$cores
wnd.size <- opt$wnd.size


if(ref.name == "NULL") ref.name = ''
if(is.null(ref.name)) ref.name = ''


if (!is.null(opt$path.features.msa)) path.features.msa <- opt$path.features.msa
if(!dir.exists(path.features.msa)) stop('features/msa dir doesn’t exist')

pokaz(path.features.msa)

if (!is.null(opt$path.inter.msa)) path.inter.msa <- opt$path.inter.msa
if(!dir.exists(path.inter.msa)) stop('internal/msa dir doesn’t exist')

if (!is.null(opt$path.figures)) path.figures <- opt$path.figures
if(!dir.exists(path.figures)) stop('Consensus folder doesn’t exist')

# ---- Combinations of chromosomes query-base to create the alignments ----
s.pattern <- paste0("^", aln.type, ".*h5")
s.combinations <- list.files(path = path.features.msa, pattern = s.pattern, full.names = FALSE)
pokaz(s.combinations)
s.combinations = gsub(aln.type, "", s.combinations)
pokaz(s.combinations)
s.combinations = gsub(".h5", "", s.combinations)
pokaz(s.combinations)


pokaz('Reference:', ref.name)
if(ref.name != ""){
  ref.suff = paste0('_', ref.name)
  
  pokaz('Reference:', ref.name)
  s.combinations <- s.combinations[grep(ref.suff, s.combinations)]
  s.combinations = gsub(ref.suff, "", s.combinations)
  
} else {
  ref.suff = ''
}

if(length(s.combinations) == 0){
  # save(list = ls(), file = "tmp_workspace_s.RData")
  stop('No Combinations found.')

} else {
  pokaz('Combinations', s.combinations)  
}
# ***********************************************************************
# ---- MAIN program body ----

file.blocks = paste0(path.inter.msa, aln.type, 'syn_blocks', ref.suff,'.rds')

if(!file.exists(file.blocks)){
  df.all = c()
  for(s.comb in s.combinations){
    
    i.chr = as.numeric(strsplit(s.comb, '_')[[1]][1])
    pokaz('Chromosome', i.chr)
    # --- --- --- --- --- --- --- --- --- --- ---
    
    file.comb.in = paste0(path.features.msa, aln.type, s.comb, ref.suff,'.h5')
    
    groups = h5ls(file.comb.in)
    accessions = groups$name[groups$group == gr.accs.b]
    
    myCluster <- makeCluster(num.cores, type = "PSOCK") 
    registerDoParallel(myCluster) 
    
    df = foreach(acc = accessions, .packages = c('rhdf5', 'crayon', 'pannagram'), .combine = rbind) %dopar% {
      pokaz('Accession', acc)
      v = h5read(file.comb.in, paste0(gr.accs.e, acc))
      
      df.acc = getBlocks(v, f.split = F)
      df.acc$acc = acc
      df.acc$chr = i.chr
      
      df.acc
    }
    
    df.all = rbind(df.all, df)
    
    stopCluster(myCluster)
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
n.chr = max(df.all$chr)

for(i.chr in 1:n.chr){
  pokaz('Chromosome', i.chr)
  i.order = 1:length(accessions)
  
  # i.order = c(2,  6,  5,  4,  3,  7,  9, 10,  8, 11,  1, 12)
  # pokaz(i.order)
  
  df.tmp = df.all
  df.tmp$acc <- factor(df.tmp$acc, levels = accessions[i.order])
  
  # save(list = ls(), file = "tmp_workspace_blocks.RData")
  
  p = panplot(df.tmp, i.chr, accessions = accessions, i.order = i.order, wnd.size=wnd.size) 
  
  pdf(paste(path.figures, 'fig_synteny_chr',i.chr,'.pdf', sep = ''), width = 6, height = 4 / 27 * length(accessions))
  print(p)     # Plot 1 --> in the first page of PDF
  dev.off()
}


warnings()

