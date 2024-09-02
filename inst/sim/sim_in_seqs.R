# Load the necessary library
library(optparse)
source(system.file("sim/sim_func.R", package = "pannagram"))
source(system.file("utils/utils.R", package = "pannagram"))

# Define options
option_list = list(
  make_option(c("--in_file"), type = "character", default = NULL,
              help = "Path to the fasta file with sequences", metavar = "FILE"),
  make_option(c("--db_file"), type = "character", default = NULL,
              help = "Path to the database fasta file with sequences", metavar = "FILE"),
  make_option(c("--res"), type = "character", default = NULL,
              help = "Path to the BLAST results", metavar = "FILE"),
  make_option(c("--out"), type = "character", default = NULL,
              help = "Path to the output coverage file", metavar = "FILE"),
  make_option(c("--use_strand"), type = "character", default = NULL,
              help = "Use strand or not", metavar = "FILE"),
  make_option(c("--sim"), type = "numeric", default = 90,
              help = "Similarity threshold", metavar = "NUMBER"),
  make_option(c("--coverage"), type = "numeric", default = NULL,
              help = "Coverage threshold", metavar = "NUMBER")
)

# Create the option parser
opt_parser = OptionParser(option_list = option_list)

# Parse the arguments
opt = parse_args(opt_parser)

# Check for the presence of all required arguments
fasta.file <- ifelse(!is.null(opt$in_file), opt$in_file, 
                     stop("FASTA file not specified", call. = FALSE))
db.fasta.file <- ifelse(!is.null(opt$db_file), opt$db_file, 
                     stop("FASTA file not specified", call. = FALSE))
blast.file <- ifelse(!is.null(opt$res), opt$res, 
                     stop("BLAST file not specified", call. = FALSE))
output.file <- ifelse(!is.null(opt$out), opt$out, 
                      stop("Output file not specified", call. = FALSE))
sim.cutoff <- ifelse(!is.null(opt$sim), opt$sim, 
                     stop("Similarity threshold not specified", call. = FALSE))
sim.cutoff = as.numeric(sim.cutoff) / 100

use.strand <- ifelse(!is.null(opt$use_strand), as.logical(opt$use_strand), 
                     stop("Strand should be provided", call. = FALSE))

coverage <- ifelse(is.null(opt$coverage), sim.cutoff, opt$coverage/100)


# ---- Testing ----

# blast.file = 'tmp.txt.blast.tmp'
# source("pannagram/utils.R")
# fasta.file = 'new_genes/new_genes.fasta'
# sim.cutoff = 0.85


# ---- Main ----

pokaz(blast.file)

v = read.table(blast.file, stringsAsFactors = F)
v = v[v$V6 >= sim.cutoff * 100,]
v = v[v$V1 != v$V8,]

# Lengths of query sequences
seqs = readFastaMy(fasta.file)
q.len = nchar(seqs)
rm(seqs)

# Lengths of database sequences
if(db.fasta.file == fasta.file){
  pokaz('Search seqeunces on themselvs')
  db.len = q.len
} else {
  seqs = readFastaMy(db.fasta.file)
  db.len = nchar(seqs)
  rm(seqs)  
}

head(db.len)

res = findNestedness(v, use.strand=use.strand)

pokaz('Number of pairs before the sumilarity cutoff', nrow(res))

# s.touched = unique(res$V1)
res$len1 = q.len[res$V1]
res$len8 = db.len[res$V8]

res$p1 = res$C1 / res$len1
res$p8 = res$C8 / res$len8

res$cover = ((res$p1 >= sim.cutoff) | (res$p8 >= sim.cutoff)) * 1
pokaz('Number of pairs after the sumilarity cutoff', sum(res$cover == 0))
head(res)

saveRDS(res, output.file)


# # Sort V4 and V5 positions
# idx.tmp = res$V4 > res$V5
# tmp = res$V4[idx.tmp]
# res$V4[idx.tmp] = res$V5[idx.tmp]
# res$V5[idx.tmp] = tmp
# 
# blastres2gff(res, output.file)
# 
# # write.table(res, output.file, quote = F, row.names = F, col.names = T, sep = '\t')






