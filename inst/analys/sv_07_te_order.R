# Structural order of Mobile Element Families (MEF) found in SVs:
# LTR (terminal direct repeats), TIR (terminal inverted repeats)
# or a one-sided poly-A tail (LINE-like).
# Everything is inferred from the SV sequences themselves.

suppressMessages({
  library(optparse)
  library(pannagram)
})

args = commandArgs(trailingOnly = TRUE)

option_list = list(
  make_option("--path.sv",        type = "character", default = NULL, help = "Path to sv dir"),
  make_option("--file.families",  type = "character", default = NULL, help = "File with SV families (SV <tab> family)"),
  make_option("--file.seqs",      type = "character", default = NULL, help = "FASTA with SV sequences"),
  make_option("--file.out",       type = "character", default = NULL, help = "Output table: one row per family"),
  make_option("--file.out.repr",  type = "character", default = NULL, help = "Output table: one row per representative"),
  make_option("--min.ltr.len",    type = "integer",   default = 100,  help = "Minimum length of LTR"),
  make_option("--min.tir.len",    type = "integer",   default = 10,   help = "Minimum length of TIR"),
  make_option("--min.polya.len",  type = "integer",   default = 10,   help = "Minimum length of the poly-A tail"),
  make_option("--term.motifs",    type = "character", default = "CACTA,CACTG",
              help = "Known short terminal signatures, comma-separated; empty string switches the check off"),
  make_option("--min.ident",      type = "double",    default = 0.8,  help = "Minimum identity of terminal repeats"),
  make_option("--hel.stem",       type = "integer",   default = 8,
              help = "Minimum length of the stem of the Helitron hairpin"),
  make_option("--hel.dist",       type = "integer",   default = 40,
              help = "Maximum distance between the Helitron hairpin and the 3'-terminus"),
  make_option("--ltr.offset",     type = "integer",   default = 20,
              help = "Maximum distance between the LTR and the terminus of the SV"),
  make_option("--tir.offset",     type = "integer",   default = 10,
              help = "Maximum distance between the TIR and the terminus of the SV"),
  make_option("--tir.ext.offset", type = "integer",   default = 0,
              help = "Shifts of the termini to try during the ungapped extension; every pair of shifts is one more test, so the rate of random hits grows quickly"),
  make_option("--n.repr",         type = "integer",   default = 10,   help = "Maximum number of representatives per family"),
  make_option("--len.frac",       type = "double",    default = 0.9,  help = "Minimum length of a representative relative to the longest member"),
  make_option("--min.support",    type = "double",    default = 1/3,  help = "Minimum fraction of representatives supporting the class")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, args = args)

# ***********************************************************************
# ---- Paths and files ----

path.sv = opt$path.sv
if(is.null(path.sv)) stop('Path to the SV dir is not provided')
if(!dir.exists(path.sv)) stop(paste0('No SV dir! ', path.sv))

file.seqs = if(is.null(opt$file.seqs)) paste0(path.sv, 'seq_sv_large.fasta') else opt$file.seqs
if(!file.exists(file.seqs)){
  pokazAttention('No SVs found:', file.seqs)
  quit(save = "no")
}

file.families = opt$file.families
if(is.null(file.families)) stop('File with families is not provided')
if(!file.exists(file.families)){
  pokazAttention('No families found:', file.families)
  pokazAttention('Please run the step -sv_families first')
  quit(save = "no")
}

file.out = if(is.null(opt$file.out)) paste0(path.sv, 'sv_te_order.txt') else opt$file.out
file.out.repr = if(is.null(opt$file.out.repr)) sub('\\.txt$', '_repr.txt', file.out) else opt$file.out.repr

# ***********************************************************************
# ---- Reading the data ----

sv.seqs = readFasta(file.seqs)
pokaz('Number of SV sequences:', length(sv.seqs))

families = read.table(file.families, sep = '\t', header = FALSE,
                      stringsAsFactors = FALSE, quote = '', comment.char = '')
if(ncol(families) < 2) stop('The file with families should have two columns')
pokaz('Number of SVs in families:', nrow(families))
pokaz('Number of families:', length(unique(families[[2]])))

# ***********************************************************************
# ---- Structural order ----

term.motifs = trimws(strsplit(opt$term.motifs, ',')[[1]])
term.motifs = term.motifs[nchar(term.motifs) > 0]
if(length(term.motifs) == 0){
  term.motifs = character(0)
  pokaz('Known terminal signatures are not checked')
} else {
  pokaz('Known terminal signatures:', paste0(term.motifs, collapse = ', '))
}

res = svTeOrderFamilies(seqs = sv.seqs,
                        families = families,
                        n.repr = opt$n.repr,
                        len.frac = opt$len.frac,
                        min.support = opt$min.support,
                        ltr.par   = list(min.len = opt$min.ltr.len, min.ident = opt$min.ident,
                                         max.offset = opt$ltr.offset),
                        tir.par   = list(min.len = opt$min.tir.len, min.ident = opt$min.ident,
                                         max.offset = opt$tir.offset, ext.offset = opt$tir.ext.offset),
                        polya.par = list(min.len = opt$min.polya.len, min.ident = opt$min.ident),
                        motif.par = list(motifs = term.motifs),
                        hel.par   = list(min.stem = opt$hel.stem, max.dist = opt$hel.dist))

# ***********************************************************************
# ---- Saving ----

write.table(res$families, file.out, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(res$repr, file.out.repr, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

pokaz('Table of families:', file.out)
pokaz('Table of representatives:', file.out.repr)

cnt = table(res$families$te.class)
for(s.class in names(cnt)){
  pokaz('Families of the order', s.class, ':', cnt[[s.class]])
}
