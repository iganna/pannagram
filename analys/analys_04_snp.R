suppressMessages({ library(Biostrings)
  library(rhdf5)
  library('seqinr')
  library('foreach')
  library(doParallel)
  library("optparse")
})

source("utils/utils.R")
source("analys/analys_func.R")

pokazStage('SNP calling')

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--ref.pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("--path.chromosomes"), type="character", default=NULL, 
              help="path to directory with chromosomes", metavar="character"),
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to directory with the consensus", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("--aln.type"), type="character", default="default", 
              help="type of alignment ('msa_', 'comb_', 'v_', etc)", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# Set the number of cores for parallel processing
num.cores.max = 10
num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))
myCluster <- makeCluster(num.cores, type = "PSOCK")
registerDoParallel(myCluster)

# Reference genome
if (is.null(opt$ref.pref)) {
  stop("ref.pref is NULL")
} else {
  ref.pref <- opt$ref.pref
}

# Alignment prefix
if (!is.null(opt$aln.type)) {
  aln.type = opt$aln.type
} else {
  aln.type = 'msa_'
}

if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons

path.seq = paste(path.cons, 'seq/', sep = '')

path.snp = paste(path.cons, 'snps/', sep = '')
if (!dir.exists(path.snp)) dir.create(path.snp)

gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'


# ---- Combinations of chromosomes query-base to create the alignments ----


s.pattern <- paste("^", aln.type, ".*", '_ref_', ref.pref, sep = '')
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type, "", files)
pref.combinations <- sub("_ref.*$", "", pref.combinations)
pref.combinations <- pref.combinations[grep("^[0-9]+_[0-9]+$", pref.combinations)]

pokaz('Reference:', ref.pref)
if(length(pref.combinations) == 0){
  stop('No Combinations found.')
} else {
  pokaz('Combinations', pref.combinations)  
}


# ---- Main ----
flag.for = T
for(s.comb in pref.combinations){
  
  pokaz('* Combination', s.comb)
  
  # Get Consensus
  i.chr = comb2ref(s.comb)
  file.seq.cons = paste(path.seq, 'seq_cons_', i.chr, '.fasta', sep = '')
  s.pangen = readFastaMy(file.seq.cons)
  s.pangen.name = names(s.pangen)[1]
  s.pangen = seq2nt(s.pangen)
  
  # Get accessions
  file.seq = paste(path.seq, 'seq_', s.comb,'_ref_',ref.pref,'.h5', sep = '')
  
  groups = h5ls(file.seq)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)

  pos.diff = c()
  for(acc in accessions){
    pokaz('Sequence of accession', acc)
    v = h5read(file.seq, paste(gr.accs.e, acc, sep = ''))
    
    pos = which((v != s.pangen) & (v != '-'))
    pos.diff = c(pos.diff, pos)
    
    rmSafe(v)
    gc()
  }
  
  # Common positions
  pos = sort(unique(pos.diff))
  
  snp.matrix = s.pangen[pos]
  snp.val = c()
  for(acc in accessions){
    pokaz('Sequence of accession', acc)
    v = h5read(file.seq, paste(gr.accs.e, acc, sep = ''))
    
    tmp = (v[pos] != s.pangen[pos]) * 1
    tmp[v[pos] == '-'] = -1
    snp.matrix = cbind(snp.matrix, tmp)
    
    snp.val = cbind(snp.val, v[pos])
    
    rmSafe(v)
    gc()
  }
  colnames(snp.val) = accessions
  
  print(dim(snp.matrix))
  print(length(c(s.pangen.name, accessions)))
  colnames(snp.matrix) = c(s.pangen.name, accessions)
  snp.matrix = cbind(pos, snp.matrix)

  file.snps = paste(path.snp, 'snps_', s.comb,'_ref_',ref.pref,'_pangen.txt', sep = '')
  write.table(snp.matrix, file.snps, row.names = F, col.names = T, quote = F, sep = '\t')
  
  
  #Save VCF-file
  
  file.vcf = paste(path.snp, 'snps_', s.comb,'_ref_',ref.pref,'_pangen.vcf', sep = '')
  saveVCF(snp.val, pos, chr.name=paste('PanGen_Chr', i.chr, sep = ''), file.vcf = file.vcf)
  
  
}

stopCluster(myCluster)

warnings()














