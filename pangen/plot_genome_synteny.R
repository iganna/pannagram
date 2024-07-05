source("visualisation/visualisation.R")
source("pangen/plot_genome_synteny_funcs.R")
source("utils/utils.R")

suppressMessages({
    require("optparse")
    library("parallel")
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
        metavar = "character"
    ),
    make_option(c("--cores"),
        type = "integer",
        default = 1,
        help = "Number of cores to use",
        metavar = "character"
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

if (
    is.null(opt$ref)      ||
    is.null(opt$path_ref) ||
    is.null(opt$path_in)  ||
    is.null(opt$path_out) ||
    is.null(opt$algn_path)||
    is.null(opt$cores)
) stop("All mandatory arguments must be provided.")

ref <- opt$ref
path.ref <- opt$path_ref
path.in <- opt$path_in
path.out <- opt$path_out
algn.path <- opt$algn_path
cores <- opt$cores

query.ids <- getPrefixes(algn.path)
pokaz("Number of pdf files:", length(query.ids))

ref.name <- findFastaFile(path.ref, ref)
if (is.null(ref.name)) stop(paste("Reference file not found for", ref))

pdf.path <- normalizePath(file.path(path.out, paste0("plots_", ref)), mustWork = FALSE)
dir.create(pdf.path, showWarnings = FALSE, recursive = TRUE)


cores <- min(cores, length(query.ids))
cl <- makeCluster(cores)
clusterExport(cl, ls(envir = .GlobalEnv), envir = .GlobalEnv)
invisible(clusterEvalQ(cl, {
    source("visualisation/visualisation.R")
    source("pangen/plot_genome_synteny_funcs.R")
    source("utils/utils.R")
    NULL
}))


results <- parLapply(cl, query.ids, function(id) {
  query.name <- findFastaFile(path.in, id)
  if (is.null(query.name)) {
    return(paste("Query file not found for", id))
  }
  
  pdf.name <- paste0(ref, "-", id)
  p <- plotGenomeAgainstRef(algn.path, query.name, ref.name)
  savePDF(p, path = pdf.path, name = pdf.name)
  
  return(paste("Successfully processed", id))
})

stopCluster(cl)

for (result in results) {
  if (grepl("Query file not found", result)) {
    warning(result)
  } else {
    # message(result)
    NULL
  }
}
