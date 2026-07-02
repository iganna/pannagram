# Load the necessary library
library(optparse)
source(system.file("sim/sim_func.R", package = "pannagram"))
source(system.file("utils/utils.R", package = "pannagram"))

# Define options
option_list = list(
  make_option(c("--in_file"), type = "character", default = NULL,
              help = "Path to the fasta file with sequences", metavar = "FILE"),
  make_option(c("--res"), type = "character", default = NULL,
              help = "Path to the BLAST results", metavar = "FILE"),
  make_option(c("--out"), type = "character", default = NULL,
              help = "Path to the output coverage file", metavar = "FILE"),
  make_option(c("--sim"), type = "numeric", default = 90,
              help = "Similarity threshold", metavar = "NUMBER"),
  make_option(c("--coverage"), type = "numeric", default = NULL,
              help = "Coverage threshold", metavar = "NUMBER"),
  make_option(c("--symmetric"), type = "character", default = "FALSE",
              help = "Require the genome footprint to also be <= 1/cov of the query length (TRUE/FALSE). FALSE keeps copies carrying insertions / nested elements. Default FALSE.", metavar = "BOOL")
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

coverage <- ifelse(is.null(opt$coverage), sim.cutoff, opt$coverage/100)

symmetric <- toupper(as.character(opt$symmetric)) %in% c("TRUE", "T", "1", "YES", "Y")

# ---- Main ----

# Rename the output file
output.file = paste(output.file, round(sim.cutoff * 100), round(coverage * 100), sep = '_')

# PATCHED (speed): data.table::fread reads the (often >100k-row) BLAST table
# ~10x faster than read.table; fall back to readBlast if data.table is absent.
if(requireNamespace("data.table", quietly = TRUE)){
  v = data.table::fread(blast.file, header = FALSE, data.table = FALSE,
                        showProgress = FALSE)
} else {
  v = readBlast(blast.file)
}
if(is.null(v) || nrow(v) == 0) quit(save = "no")
# PATCHED: do NOT pre-filter HSPs by pident here. Identity is now enforced on the
# fully assembled copy inside findHitsInRef (length-weighted), so low-identity
# internal pieces can still fill coverage gaps instead of being dropped early.

len1 = v$V9
# PATCHED: capture the optional mismatch (V11) / gaps (V12) columns (present when
# blastn is run with `... mismatch gaps` appended to outfmt 6) so findHitsInRef
# can compute RM-style substitution-only identity. Falls back gracefully if absent.
mism = if(ncol(v) >= 11) v$V11 else NULL
gaps = if(ncol(v) >= 12) v$V12 else NULL
v = v[,1:8,drop=F]
v$len1 = len1
if(!is.null(mism) && !is.null(gaps)){
  v$mism = mism
  v$gaps = gaps
}

# ---- Similarity analysis ----
res = findHitsInRef(v, sim.cutoff = sim.cutoff, coverage = coverage,
                    symmetric = symmetric, echo = F)

# Sort V4 and V5 positions
idx.tmp = res$V4 > res$V5
tmp = res$V4[idx.tmp]
res$V4[idx.tmp] = res$V5[idx.tmp]
res$V5[idx.tmp] = tmp

if(nrow(res) == 0){
  quit(save = "no")
}

# PATCHED GFF: emit the MERGED genome sub-intervals of each copy (one line per
# aligned segment, sharing the copy ID = Q<i>) instead of the min..max envelope,
# so inserted / nested bp between segments is NOT painted as covered.
gff.rows <- list()
for(i in seq_len(nrow(res))){
  segs <- strsplit(res$subiv[i], ";", fixed = TRUE)[[1]]
  cov.pct <- round(res$V7[i] / res$len1[i] * 100, 1)
  for(s in segs){
    ab <- as.numeric(strsplit(s, "-", fixed = TRUE)[[1]])
    gff.rows[[length(gff.rows) + 1]] <- data.frame(
      col1 = res$V8[i], col2 = 'blast2gff', col3 = 'query',
      col4 = min(ab), col5 = max(ab), col6 = '.', col7 = res$strand[i], col8 = '.',
      col9 = paste0('ID=Q', i, ';query=', res$V1[i], ';length=', res$len1[i],
                    ';similarity=', round(res$V6[i], 1), ';coverage=', cov.pct),
      stringsAsFactors = FALSE)
  }
}
gff.df <- do.call(rbind, gff.rows)
gff.df <- gff.df[order(gff.df$col4), ]
gff.df <- gff.df[order(gff.df$col1), ]
writeGFF(gff.df, paste0(output.file, '.gff'))

# Copy-level table/cnt: drop the helper sub-interval column.
res$subiv <- NULL
colnames(res) <- c('sequence', 'beg.seq', 'end.seq', 'beg.genome', 'end.genome', 'similarity', 'coverage', 'genome.chr', 'len.seq', 'strand')
res = res[, c('sequence', 'beg.seq', 'end.seq', 'beg.genome', 'end.genome', 'strand', 'genome.chr', 'similarity', 'coverage', 'len.seq')]
res$coverage = res$coverage / res$len.seq * 100

res$similarity = round(res$similarity, 1)
res$coverage = round(res$coverage, 1)

res = res[order(res$genome.chr),]
res = res[order(res$sequence),]
write.table(res, paste0(output.file, '.table'), quote = F, row.names = F, col.names = T, sep = '\t')

# Copy-number information
res.cnt = as.data.frame.matrix(table(res$sequence, res$genome.chr))
res.cnt$total = rowSums(res.cnt)
res.cnt = res.cnt[order(-res.cnt$total),]
write.table(res.cnt, paste0(output.file, '.cnt'), quote = F, row.names = T, col.names = NA, sep = '\t')

cnt = tapply(res$genome.chr, res$sequence, length)
pokaz('Mean, min, max number of hits per sequence:', mean(cnt),  min(cnt),  max(cnt))
pokaz('Number of hits found:', nrow(res))
