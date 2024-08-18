# This script gets gff file and the genome and return the set of candidate sequences for merging

path.pannagram = ""
source(paste(path.pannagram, 'utils/utils.R', sep = ''))


pokazStage('Analyse counts..')

library(optparse)

option_list = list(
  make_option(c("--file.cnt"), type="character", default="", 
              help="Path to the GFF file", metavar="character"),
  make_option(c("--file.genome"), type="character", default="", 
              help="Path to the genome file", metavar="character"),
  make_option(c("--file.seqs"), type="character", default="", 
              help="Path to the sequences file", metavar="character"),
  make_option(c("--file.fix"), type="character", default="", 
              help="Path to the sequences file", metavar="character"),
  make_option(c("--copy.number"), type="integer", default=3, 
              help="Allowed minimal copu-number of genes", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# print(opt)

# ---- Parsing of parameters ----
opt = parse_args(opt_parser)

if (is.null(opt$file.cnt)) {
  stop("Error: --file.cnt is required.")
}
if (is.null(opt$file.genome)) {
  stop("Error: --file.genome is required.")
}
if (is.null(opt$file.seqs)) {
  stop("Error: --file.seqs is required.")
}
if (is.null(opt$file.fix)) {
  stop("Error: --file.fix is required.")
}

file.cnt = opt$file.cnt
file.genome = opt$file.genome
file.seqs = opt$file.seqs
file.fix = opt$file.fix
copy.number = opt$copy.number


# ---- Read the merged sequences
file.fix = '/Users/annaigolkina/Library/CloudStorage/OneDrive-Personal/vienn/pacbio/tebagra/tair12_data/simsearch_90/merged_seqs_fixed.txt'

x = read.table(file.fix, stringsAsFactors = F)

x.info = t(sapply(gsub("Chr", "", x$V2), function(s) as.numeric(strsplit(s, '\\|')[[1]][3:5])))
rownames(x.info) = NULL
colnames(x.info) = c('chr', 'beg', 'end')
x = cbind(x, x.info)

x = x[order(x$end),]
x = x[order(x$beg),]
x = x[order(x$chr),]
rownames(x) = NULL

to.remove <- c()

for (i in 2:nrow(x)) {
  for (j in 1:(i-1)) {
    if ((x$chr[i] == x$chr[j]) &&
        (x$beg[j] <= x$beg[i]) && (x$end[i] <= x$end[j])) {
      to.remove <- cbind(to.remove, c(i, j))
      break
    }
  }
}

# Read all of the sequence files and take those sequences, which match in  names

