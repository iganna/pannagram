# SNP calling

suppressMessages({ library(Biostrings)
  library(rhdf5)
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("analys/analys_func.R", package = "pannagram"))

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option("--path.features.msa", type = "character", default = NULL, help = "Path to msa dir (features)"),
  make_option("--path.snp", type = "character", default = NULL, help = "Path to snp dir"),
  make_option("--path.seq", type = "character", default = NULL, help = "Path to seq dir"),
  make_option(c("-c", "--cores"), type = "integer", default = 1,       help = "number of cores to use for parallel processing"),
  make_option(c("--aln.type"),    type="character", default="default", help="type of alignment ('msa_', 'comb_', 'v_', etc)"),
  make_option(c("--ref.pref"),    type="character", default=NULL,      help="prefix of the reference file")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);


path.seq <- opt$path.seq
if (!dir.exists(path.seq)) stop('Folder `path.seq` does not exist')

path.snp <- opt$path.snp
if (!dir.exists(path.snp)){
  dir.create(path.snp)
} 
if (!dir.exists(path.snp)){
  stop(paste0('The output folder was not created'))
} 

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files


# ***********************************************************************

# Set the number of cores for parallel processing
num.cores <- opt$cores
if(num.cores > 1){
  myCluster <- makeCluster(num.cores, type = "PSOCK")
  registerDoParallel(myCluster)  
}

# Reference genome
if (is.null(opt$ref.pref) || (opt$ref.pref == "NULL")) {
  ref.pref <- ""
  ref.suff = ""
} else {
  ref.pref <- opt$ref.pref
  ref.suff <- paste0('_ref_', ref.pref)
}

# Alignment prefix
if (!is.null(opt$aln.type)) {
  aln.type = opt$aln.type
} else {
  aln.type = aln.type.msa
  pokazAttention('The defaul anighment type is used:', aln.type)
}

path.features.msa <- if (!is.null(opt$path.features.msa)) opt$path.features.msa else stop("Error: 'path.features.msa' is NULL. Please provide a valid path.")

# Paths
if(!dir.exists(path.features.msa)) stop(paste('The consensus folder does not exist:', path.features.msa))

# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type, ".*", ref.suff, "\\.h5")

# pokaz('Consensus folder', path.cons)
files <- list.files(path = path.features.msa, pattern = s.pattern, full.names = FALSE)

pref.combinations = gsub(aln.type, "", files)
pref.combinations <- sub(ref.suff, "", pref.combinations)
pref.combinations <- sub(".h5", "", pref.combinations)

if(length(pref.combinations) == 0){
  stop('No Combinations found.')
} else {
  pokaz('Combinations', pref.combinations)  
}

# ---- Main ----
flag.for = T
for(s.comb in pref.combinations){
  
  # Get Consensus
  i.chr = comb2ref(s.comb)
  file.seq.cons = paste0(path.seq, 'seq_cons_', i.chr, '.fasta')
  s.pangen = readFastaMy(file.seq.cons)
  s.pangen.name = names(s.pangen)[1]
  s.pangen = seq2nt(s.pangen)
  
  # Get accessions
  file.seq = paste0(path.seq, 'seq_', s.comb, ref.suff, '.h5')
  
  groups = h5ls(file.seq)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)

  pokaz('Round 1: get positions of differences..')
  pos.diff = c()
  for(acc in accessions){
    # pokaz('Sequence of accession', acc)
    v = h5read(file.seq, paste0(gr.accs.e, acc))
    
    pos = which((v != s.pangen) & (v != '-'))
    pos.diff = c(pos.diff, pos)
    
    rmSafe(v)
    gc()
  }
  
  # Common positions
  pokaz('Round 2: get diffs in common positions..')
  pos = sort(unique(pos.diff))
  
  snp.matrix = s.pangen[pos]
  snp.val = c()
  for(acc in accessions){
    # pokaz('Sequence of accession', acc)
    v = h5read(file.seq, paste0(gr.accs.e, acc))
    
    tmp = (v[pos] != s.pangen[pos]) * 1
    tmp[v[pos] == '-'] = -1
    snp.matrix = cbind(snp.matrix, tmp)
    
    snp.val = cbind(snp.val, v[pos])
    
    rmSafe(v)
    gc()
  }
  colnames(snp.val) = accessions
  
  # print(dim(snp.matrix))
  # print(length(c(s.pangen.name, accessions)))
  colnames(snp.matrix) = c(s.pangen.name, accessions)
  snp.matrix = cbind(pos, snp.matrix)

  pokaz('Save table...')
  file.snps = paste0(path.snp, 'snps_', s.comb, ref.suff, '_pangen.txt')
  write.table(snp.matrix, file.snps, row.names = F, col.names = T, quote = F, sep = '\t')
  
  #Save VCF-file
  pokaz('Save VCF-file...')
  file.vcf = paste0(path.snp, 'snps_', s.comb, ref.suff, '_pangen.vcf')
  saveVCF(snp.val, pos, chr.name=paste0('PanGen_Chr', i.chr), file.vcf = file.vcf)
  
  # Create the VCF-file for the first accession, the main reference.
  file.comb = paste0(path.features.msa, aln.type, s.comb, ref.suff, '.h5')
  # pokaz(file.comb)
  acc = accessions[1]
  pos.acc = h5read(file.comb, paste0(gr.accs.e, acc))
  pos.acc = pos.acc[pos]
  snp.val.acc = snp.val[pos.acc != 0,,drop=F]
  pos.acc = abs(pos.acc[pos.acc != 0])

  pokaz('Save VCF-file for the accession', acc, '...')
  file.vcf.acc = paste0(path.snp, 'snps_', s.comb, ref.suff, '_',acc,'.vcf')
  saveVCF(snp.val.acc, pos.acc, chr.name=paste0(acc,'_Chr', i.chr), file.vcf = file.vcf.acc)
  
}

if(num.cores > 1){
  stopCluster(myCluster)
}

warnings()
