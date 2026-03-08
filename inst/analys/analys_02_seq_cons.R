# Get sequence alignments and consensus sequence

suppressMessages({
  library(Biostrings)
  library(rhdf5)
  library(foreach)
  library(doParallel)
  library(optparse)
  library(pannagram)
  library(crayon)
})


# ***********************************************************************
# ---- Alignment types ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) 

# ***********************************************************************
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option("--ref",               type = "character", default = "",           help = "Prefix of the reference file"),
  make_option("--path.chr",          type = "character", default = NULL,         help = "path to directory with chromosomes"),
  make_option("--path.seq",          type = "character", default = NULL,         help = "Path to seq dir"),
  make_option("--path.features.msa", type = "character", default = NULL,         help = "Path to msa dir (features)"),
  make_option("--cores",             type = "integer",   default = 1,            help = "number of cores to use for parallel processing"),
  make_option("--aln.type",          type = "character", default = aln.type.msa, help = "type of alignment ('pan', 'ref', 'extra1', etc)"),
  make_option("--log.level",        type = "character", default = NULL, help = "Level of log to be shown on the screen")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ***********************************************************************
# ---- Variables ----

num.cores <- opt$cores  # Set the number of cores for parallel processing

s.nts = c('A', 'C', 'G', 'T', '-')

# ***********************************************************************
# ---- Paths ----

path.seq <- opt$path.seq
if (!dir.exists(path.seq)) stop(paste0('No path.seq dir found!'))

path.features.msa <- opt$path.features.msa
if(!dir.exists(path.features.msa)) stop('features/msa dir does not exist')

path.chr <- opt$path.chr
if(!dir.exists(path.chr)) stop('intermediate/chromosomes dir does not exist')

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

# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(s.comb, echo = T){
# tmp = foreach(s.comb = pref.combinations, .packages=c('rhdf5', 'crayon'))  %dopar% {  # which accession to use
# # for(s.comb in pref.combinations){
  pokaz('* Combination', s.comb)
  
  # Get accessions
  file.comb = paste0(path.features.msa, aln.pref, s.comb, ref.suff, '.h5')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  if(ref.name %in% accessions){
    accessions = c(ref.name, setdiff(accessions, ref.name))
  }
  n.acc = length(accessions)
  
  # File with sequences
  file.seq = paste0(path.seq, 'seq_', s.comb, ref.suff,'.h5')
  if (file.exists(file.seq)) file.remove(file.seq)
  h5createFile(file.seq)
  h5createGroup(file.seq, gr.accs.e)
  
  mx.consensus = NULL
  # idx.negative = c()
  for(acc in accessions){
    # pokaz('Sequence of accession', acc)
    v = h5read(file.comb, paste0(gr.accs.e, acc))
    v.na = is.na(v)
    v[v.na] = 0
    if(is.null(mx.consensus)){
      mx.consensus = matrix(0, nrow = length(v), ncol = length(s.nts), dimnames = list(NULL, s.nts))
    }
    
    if(acc == ref.name){
      q.chr = strsplit(s.comb, '_')[[1]][2]
    } else {
      q.chr = strsplit(s.comb, '_')[[1]][1]  
    }
    
    pokaz('Accession', acc, 'Chromosome', q.chr)
    
    file.chr = paste0(path.chr, acc, '_chr', q.chr, '.fasta')
    if(!file.exists(file.chr)){
      stop(paste0('Chromosomal file was not found', file.chr))
    }
    genome = readFasta(file.chr)
    genome = seq2nt(genome)
    genome = toupper(genome)
    
    if(max(abs(v)) > length(genome)) stop('Length of the genome is shorter than the idex involded')
  
    s = rep('-', length(v))
    idx.plus = (v > 0)
    idx.mins = (v < 0)
    if(sum(idx.plus) > 0){
      s[idx.plus] = genome[v[idx.plus]]
    }
    if(sum(idx.mins) > 0){
      s[idx.mins] = justCompl(genome[abs(v[idx.mins])])
    }
    
    # idx.negative = c(idx.negative, which(idx.mins))
    
    for(s.nt in s.nts){
      mx.consensus[,s.nt] = mx.consensus[,s.nt] + (s == s.nt)
    }
    
    suppressMessages({
      h5write(s, file.seq, paste0(gr.accs.e, acc))
    })
    
    rm(v)
    rm(v.na)
    rm(genome)
    rm(s)
    rm(idx.plus)
    rm(idx.mins)
    gc()
    
  }
  
  suppressMessages({
    h5write(mx.consensus, file.seq, 'matrix')
  })
  
  # ---- Consensus sequence ----
  pokaz('Prepare consensus fasta-sequence')
  i.chr = comb2ref(s.comb)
  file.seq.cons = paste0(path.seq, 'seq_cons_', s.comb, ref.suff, '.fasta')
  
  n = nrow(mx.consensus)
  s.cons = rep('N', n)
  n.nt = rep(0, n)
  for(k in 1:4){
    idx.k =  mx.consensus[,k] > n.nt
    s.cons[idx.k] = s.nts[k]
    n.nt[idx.k] = mx.consensus[idx.k, k]
  }
  
  if(sum(s.cons == 'N') > 0) pokazAttention('Some nucleotides are N:', sum(s.cons == 'N'))
  s.cons = paste0(s.cons, collapse = '')
  
  pokaz('Saving consensus sequence...')
  names(s.cons) = paste0('PanGen_Chr', i.chr)
  writeFasta(s.cons, file.seq.cons)
  
  rmSafe(mx.consensus)
  rmSafe(s.cons)
  H5close()
  gc()
  
  # pokaz('Done.', file=file.log.loop, echo=echo.loop)
}


# ***********************************************************************
# ---- Loop  ----

if(num.cores == 1){
  
  for(s.comb in s.combinations){
    loop.function(s.comb)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  foreach(s.comb = s.combinations, .packages=c('rhdf5', 'crayon', 'pannagram'))  %dopar% { 
    tmp = loop.function(s.comb)
    return(tmp)
  }
  stopCluster(myCluster)
}

warnings()

