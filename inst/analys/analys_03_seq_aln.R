# Create and Alignment in Pangenome coordinates
suppressMessages({ library(Biostrings)
  library(rhdf5)
  library('foreach')
  library(doParallel)
  library("optparse")
})

source(system.file("utils/utils.R", package = "pannagram"))



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

if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons

path.seq = paste0(path.cons, 'seq/')

path.aln = paste0(path.cons, 'aln_pangen/')
if (!dir.exists(path.aln)) dir.create(path.aln)

# Alignment prefix
if (!is.null(opt$aln.type)) {
  aln.type = opt$aln.type
} else {
  aln.type = 'msa_'
}

# ---- Variables ----

gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'


# ---- Combinations of chromosomes query-base to create the alignments ----


s.pattern <- paste0("^", 'msa_', ".*", '_ref_', ref.pref)
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub("msa_", "", files)
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
  
  # Get accessions
  file.seq = paste0(path.seq, 'seq_', s.comb,'_ref_',ref.pref,'.h5')
  
  groups = h5ls(file.seq)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)

  file.aln = paste0(path.aln, 'aln_', s.comb,'_ref_',ref.pref,'_pangen.fasta')

  for(acc in accessions){
    pokaz('Sequence of accession', acc)
    v = h5read(file.seq, paste0(gr.accs.e, acc))
    v = paste0(v, collapse='')
    names(v) = acc
    writeFastaMy(v, file.aln, append = T)

    rmSafe(v)
    gc()
  }

  
  # ---- MAF file ----
  file.aln = paste0(path.aln, 'aln_', s.comb,'_ref_',ref.pref,'_pangen.maf')

  # Add pangen line
  i.chr = comb2ref(s.comb)
  file.seq.cons = paste0(path.seq, 'seq_cons_', i.chr, '.fasta')
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
    v = h5read(file.seq, paste0(gr.accs.e, acc))
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














