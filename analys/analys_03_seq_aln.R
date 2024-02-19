suppressMessages({ library(Biostrings)
  library(rhdf5)
  library('seqinr')
  library('foreach')
  library(doParallel)
  library("optparse")
})

source("utils/utils.R")

pokazStage('Create and Alignment in Pangenome coordinates')

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--ref.pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("--path.chromosomes"), type="character", default=NULL, 
              help="path to directory with chromosomes", metavar="character"),
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to directory with the consensus", metavar="character"),
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

# Reference genome
if (is.null(opt$ref.pref)) {
  stop("ref.pref is NULL")
} else {
  ref.pref <- opt$ref.pref
}

if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons

path.seq = paste(path.cons, 'seq/', sep = '')

path.aln = paste(path.cons, 'aln_pangen/', sep = '')
if (!dir.exists(path.aln)) dir.create(path.aln)

gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'


# ---- Combinations of chromosomes query-base to create the alignments ----


s.pattern <- paste("^", 'msa_', ".*", '_ref_', ref.pref, sep = '')
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub("msa_", "", files)
pref.combinations <- sub("_ref.*$", "", pref.combinations)
pref.combinations <- pref.combinations[grep("^[0-9]+_[0-9]+$", pref.combinations)]

pokaz('Reference:', ref.pref)
pokaz('Combinations', pref.combinations)


# ---- Main ----
flag.for = T
for(s.comb in pref.combinations){
  
  pokaz('* Combination', s.comb)
  
  # Get accessions
  file.seq = paste(path.seq, 'seq_', s.comb,'_ref_',ref.pref,'.h5', sep = '')
  
  groups = h5ls(file.seq)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)
  # 
  # file.aln = paste(path.aln, 'aln_', s.comb,'_ref_',ref.pref,'_pangen.fasta', sep = '')
  # 
  # for(acc in accessions){
  #   pokaz('Sequence of accession', acc)
  #   v = h5read(file.seq, paste(gr.accs.e, acc, sep = ''))
  #   v = paste0(v, collapse='')
  #   names(v) = acc
  #   writeFastaMy(v, file.aln, append = T)
  #   
  #   rmSafe(v)
  #   gc()
  # }
  # 
  
  # ---- MAF file ----
  file.aln = paste(path.aln, 'aln_', s.comb,'_ref_',ref.pref,'_pangen.maf', sep = '')

  
  # Ddd pangen line
  i.chr = comb2ref(s.comb)
  file.seq.cons = paste(path.seq, 'seq_cons_', i.chr, '.fasta', sep = '')
  s.pangen = readFastaMy(file.seq.cons)
  len.pangen = nchar(s.pangen)
  
  
  # len.pangen = 20
  
  
  # Write initial info
  write(paste('track name=pangen',
                   '_ref_', ref.pref, 
                   '_chr_', i.chr, 
                   '\n', sep = ''), 
             file.aln)
  
  write(paste('a score=666',
              # '\n',
              sep = ''),
        file.aln, append = TRUE)
  
  
  write(paste('s ',
              names(s.pangen)[1], ' ',  # track id
                   1, ' ',
              len.pangen, ' ', 
                   '+', ' ',
              len.pangen, ' ', 
              # substr(s.pangen,1,len.pangen),
              s.pangen,
                   # '\n',
              sep = ''), 
             file.aln, append = TRUE)
  
  for(acc in accessions){
    pokaz('Sequence of accession', acc)
    v = h5read(file.seq, paste(gr.accs.e, acc, sep = ''))
    v = paste0(v, collapse='')
    
    write(paste('s ',
                acc, ' ',  # track id
                     1, ' ',
                     len.pangen, ' ', 
                     '+', ' ',
                len.pangen, ' ', 
                # substr(v,1,len.pangen),
                v,
                     # '\n',
                sep = ''), 
               file.aln, append=T)
    
    rmSafe(v)
    gc()
  }
  
}

stopCluster(myCluster)

warnings()














