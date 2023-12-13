# Load the necessary library
library(optparse)
source('sim/sim_func.R')
source('utils/utils.R')

# Define options
option_list = list(
  make_option(c("--in_file"), type = "character", default = NULL,
              help = "Path to the fasta file with sequences", metavar = "FILE"),
  make_option(c("--res"), type = "character", default = NULL,
              help = "Path to the BLAST results", metavar = "FILE"),
  make_option(c("--out"), type = "character", default = NULL,
              help = "Path to the output coverage file", metavar = "FILE"),
  make_option(c("--sim"), type = "numeric", default = 90,
              help = "Similarity threshold", metavar = "NUMBER")
)

# Create the option parser
opt_parser = OptionParser(option_list = option_list)

# Parse the arguments
opt = parse_args(opt_parser)

# Check for the presence of all required arguments
fasta.file <- ifelse(!is.null(opt$in_file), opt$in_file, 
                     stop("FASTA file not specified", call. = FALSE))
blast.file <- ifelse(!is.null(opt$res), opt$res, 
                     stop("BLAST file not specified", call. = FALSE))
output.file <- ifelse(!is.null(opt$out), opt$out, 
                      stop("Output file not specified", call. = FALSE))
sim.cutoff <- ifelse(!is.null(opt$sim), opt$sim, 
                     stop("Similarity threshold not specified", call. = FALSE))
sim.cutoff = as.numeric(sim.cutoff) / 100


# ---- Testing ----

# blast.file = 'tmp.txt.blast.tmp'
# source('pannagram/utils.R')
# fasta.file = 'new_genes/new_genes.fasta'
# sim.cutoff = 0.85


# ---- Main ----


v = read.table(blast.file, stringsAsFactors = F)
v = v[v$V6 >= sim.cutoff * 100,]

seqs = readFastaMy(fasta.file)
v$len1 = nchar(seqs)[v$V1]
rm(seqs)

res = findHitsInRef(v, sim.cutoff = sim.cutoff, echo = F)

# Sort V4 and V5 positions
idx.tmp = res$V4 > res$V5
# print(sum(idx.tmp))
# print(sum(res$strand == '-'))

tmp = res$V4[idx.tmp]
res$V4[idx.tmp] = res$V5[idx.tmp]
res$V5[idx.tmp] = tmp

blastres2gff(res, output.file)

write.table(res, paste(output.file, '.table', sep = ''), quote = F, row.names = F, col.names = T, sep = '\t')

pokaz('Number of hits found:', nrow(res))




