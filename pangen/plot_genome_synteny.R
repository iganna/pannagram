
### imports
source("visualisation/visualisation.R")
source("pangen/plot_genome_synteny_funcs.R")
source("utils/utils.R")

suppressMessages({
    require("optparse")
})

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
    make_option(c("--ref"),
        type = "character",
        default = NULL,
        help = "Name af reference sequence",
        metavar = "character"
    ),
    make_option(c("--path_ref"),
        type = "character",
        default = NULL,
        help = "Path to reference sequence",
        metavar = "character"
    ),
    make_option(c("--path_in"),
        type = "character",
        default = NULL,
        help = "Path to genome fasta files",
        metavar = "character"
    ),
    make_option(c("--path_out"),
        type = "character",
        default = NULL,
        help = "Path to the output directory",
        metavar = "character"
    ),
    make_option(c("--algn_path"),
        type = "character",
        default = NULL,
        help = "Path to the alignment directory",
        metavar = "character")    
)

opt <- parse_args(OptionParser(option_list=option_list))

if (
    is.null(opt$ref) ||
    is.null(opt$path_ref) ||
    is.null(opt$path_in) ||
    is.null(opt$path_out) ||
    is.null(opt$algn_path)
) stop("All mandatory arguments must be provided.")

ref <- opt$ref
path.ref <- opt$path_ref
path.in <- opt$path_in
path.out <- opt$path_out
algn.path <- opt$algn_path

query.ids <- getPrefixes(algn.path)
pokaz("Number of pdf files:", length(query.ids))

ref.name <- normalizePath(file.path(path.ref, paste0(ref, ".fna")), mustWork = FALSE)
pdf.path <- normalizePath(file.path(path.out, paste0("plots_", ref)), mustWork = FALSE)
dir.create(pdf.path, showWarnings = FALSE, recursive = TRUE)

for (id in query.ids){
    query.name <- normalizePath(file.path(path.in, paste0(id, ".fna")), mustWork = FALSE)
    pdf.name <- paste0(ref, "-", id)
    p <- plotGenomeAgainstRef(algn.path, query.name, ref.name, seq.order="alphanum")
    savePDF(p, path = pdf.path, name = pdf.name)
}
