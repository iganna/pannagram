#!/usr/bin/env Rscript

# Load the necessary library
library(optparse)
source('sim_func.R')
source('utils.R')

# Define options
option_list = list(
  make_option(c("-in", "--fasta_file"), type = "character", default = NULL,
              help = "Path to the fasta file with sequences", metavar = "FILE"),
  make_option(c("-res"), type = "character", default = NULL,
              help = "Path to the BLAST results", metavar = "FILE"),
  make_option(c("-out"), type = "character", default = NULL,
              help = "Path to the output coverage file", metavar = "FILE"),
  make_option(c("-sim"), type = "numeric", default = 90,
              help = "Similarity threshold", metavar = "NUMBER")
)

# Create the option parser
opt_parser = OptionParser(option_list = option_list)

# Parse the arguments
opt = parse_args(opt_parser)
print(opt)

# Check for the presence of all required arguments
fasta.file <- ifelse(!is.null(opt$fasta_file), opt$fasta_file, 
                     stop("FASTA file not specified", call. = FALSE))
blast.file <- ifelse(!is.null(opt$res), opt$res, 
                     stop("BLAST file not specified", call. = FALSE))
output.file <- ifelse(!is.null(opt$out), opt$out, 
                      stop("Output file not specified", call. = FALSE))
sim.cutoff <- ifelse(!is.null(opt$sim), opt$sim, 
                     stop("Similarity threshold not specified", call. = FALSE))
sim.cutoff = as.numeric(sim.cutoff) / 100

pokaz(blast.file)
# v = read.table(blast.file, stringsAsFactors = F)
# 
# seqs = readFastaMy(fasta.file)
# v$len1 = nchar(seqs)
# rm(seqs)
# 
# if(sum(is.na(v$len1) > 0)) {
#   stop('Length of sequences is supposed to be')
# }
# 
# res = findHitsInRef(v, sim.cutoff = sim.cutoff, echo = F)
# 
# write.table(res, output.file, quote = F, row.names = F, col.names = T, sep = '\t')






