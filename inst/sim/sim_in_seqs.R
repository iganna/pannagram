# Load the necessary library
invisible(suppressMessages({
  library(optparse)
  library(pannagram)
}))
source(system.file("sim/sim_func.R", package = "pannagram"))


# Define options
option_list = list(
  make_option(c("--in_file"),    type = "character", default = NULL, help = "Path to the fasta file with sequences"),
  make_option(c("--db_file"),    type = "character", default = NULL, help = "Path to the database fasta file with sequences"),
  make_option(c("--res"),        type = "character", default = NULL, help = "Path to the BLAST results"),
  make_option(c("--out"),        type = "character", default = NULL, help = "Path to the output coverage file"),
  make_option(c("--use_strand"), type = "character", default = NULL, help = "Use strand or not"),
  make_option(c("--sim"),        type = "numeric",   default = NULL, help = "Similarity threshold"),
  make_option(c("--coverage"),   type = "numeric",   default = NULL, help = "Coverage threshold")
);

# Create the option parser
opt_parser = OptionParser(option_list = option_list)

# Parse the arguments
opt = parse_args(opt_parser)
# print(opt)

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
use.strand <- ifelse(!is.null(opt$use_strand), as.logical(opt$use_strand), 
                     stop("Strand should be provided", call. = FALSE))
coverage <- ifelse(!is.null(opt$coverage), opt$coverage, 
                   stop("Coverage threshold not specified", call. = FALSE))

# ---- Main ----

v = readBlast(blast.file)
v = v[v$V6 >= sim.cutoff,]
v = v[v$V1 != v$V8,]

# print(v)

if(nrow(v) == 0){
  pokazAttention('No similarity in SVs, NO SVs were genegated')
  quit(save = "no", status = 0)
}

# Lengths
uniq1 = !duplicated(v$V1)
len1 = setNames(v$V9[uniq1], v$V1[uniq1])
uniq8 = !duplicated(v$V8)
len8 = setNames(v$V10[uniq8], v$V8[uniq8])
                
# Nestedness
res = findNestedness(v, use.strand=use.strand)

pokaz('Number of pairs before the sumilarity cutoff', nrow(res))

# s.touched = unique(res$V1)
res$len1 = len1[res$V1]
res$len8 = len8[res$V8]

res$p1 = res$C1 / res$len1
res$p8 = res$C8 / res$len8

# res$cover = ((res$p1 >= sim.cutoff) | (res$p8 >= sim.cutoff)) * 1
# pokaz('Number of pairs after the sumilarity cutoff', sum(res$cover == 0))

# Change the order of columns
res <- res[,c('V1', 'V8', 'dir', 'len1', 'len8', 'C1', 'C8', 'p1', 'p8')]

saveRDS(res, output.file)

# Save the txt-file with proper column names
res <- res[,c('V1', 'V8', 'dir', 'len1', 'len8', 'p1', 'p8')]
colnames(res) <- c('name.query', 'name.target', 'strand', 'length.query', 'length.target', 'coverage.query', 'coverage.target')

# output.file.txt = sub('.rds', '.txt', output.file)
# pokaz(output.file.txt)
# write.table(res,
#             output.file.txt,
#             sep       = "\t",
#             quote     = FALSE,
#             row.names = F)


# Incorporate coverage
res = res[(res$coverage.q > coverage / 100) | (res$coverage.t > coverage / 100),]

output.file.txt = sub('\\.rds$', paste0( '_',sim.cutoff, '_', coverage, '.txt'), output.file)
write.table(res,
            output.file.txt,
            sep       = "\t",
            quote     = FALSE,
            row.names = F)

# # Sort V4 and V5 positions
# idx.tmp = res$V4 > res$V5
# tmp = res$V4[idx.tmp]
# res$V4[idx.tmp] = res$V5[idx.tmp]
# res$V5[idx.tmp] = tmp
# 
# blastres2gff(res, output.file)
# 
# # write.table(res, output.file, quote = F, row.names = F, col.names = T, sep = '\t')
