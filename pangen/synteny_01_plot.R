
### imports
source("visualisation/visualisation.R")
source("pangen/synteny_func_plot.R")
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
        metavar = "character"),
    make_option(c("--path.log"), 
        type = "character", 
        default = NULL,
        help = "Path for log files",
        metavar = "character"),
    make_option(c("--log.level"), 
        type = "character", 
        default = NULL,
        help = "Level of log to be shown on the screen", 
        metavar = "character")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Logging
source('utils/chunk_logging.R') # a common code for all R logging

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
path.aln<- opt$algn_path

# Extracting only ids from path
pattern <- ".*_[0-9]+_[0-9]+_maj\\.rds$"
files.aln <- list.files(path = path.aln, pattern = pattern, full.names = F)
query.ids <- unique(sapply(files.aln, function(s) sub("^(.*?)_\\d+_\\d+_maj\\.rds$", "\\1", s)))
pokaz('Genomes analysed:', query.ids, file=file.log.main, echo=echo.main)

# Find the file with the reference genome
ext <- c('fasta', 'fna', 'fa', 'fas')
ref.name <- findGenomeFile(genome.pref = ref, 
                           path.genome = path.ref,
                           ext = ext)
if (is.null(ref.name)) stop('No reference genome files found in the specified folder')

# Output folder with plots
pdf.path <- normalizePath(file.path(path.out, paste0("plots_", ref)), mustWork = FALSE)
dir.create(pdf.path, showWarnings = FALSE, recursive = TRUE)

for (id in query.ids){
  pokaz('Accession', id, file=file.log.main, echo=echo.main)
  query.name <- findGenomeFile(genome.pref = id, 
                               path.genome = path.in,
                               ext = ext)
  if (is.null(query.name)) stop('No target genome files found in the specified folder')
  
  pdf.name <- paste0(ref, "-", id)
  p <- plotGenomeAgainstRef(path.aln, query.name, ref.name, seq.order="alphanum")
  savePDF(p, path = pdf.path, name = pdf.name)
}

pokaz('Done.',
      file=file.log.main, echo=echo.main)

