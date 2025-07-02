#' Extract a subregion from an alignment matrix
#'
#' This function extracts a region from a multiple sequence alignment (MSA)
#' or sequence coordinate matrix for a specific chromosome and accession.
#'
#' @param acc Accession name (or 'pangenome', 'pannagram' etc. for the pangenome coordinate).
#' @param i.chr Chromosome number.
#' @param p.beg Start position in the Accession.
#' @param p.end End position in the Accession.
#' @param path.proj Path to project folder.
#' @param mode Mode of extraction: either `"seq"` (sequence coordinates) or `"pos"` (alignment coordinates).
#' @param aln.type Prefix for alignment files (default `"msa_"`).
#' @param ref.acc Reference accession ID if you built in the reference-based alignment.
#' @param echoWhether to print verbose messages (default `FALSE`).
#'
#' @return A matrix where each row corresponds to an accession and columns represent aligned bases in the specified region.
#'
#' @export
cutAln <- function(acc, i.chr, p.beg, p.end,
                   path.proj = NULL,
                   mode = 'seq',
                   aln.type = "msa_", 
                   ref.acc = '',
                   echo = FALSE, ...) {
  
  # --- Variables ---
  s.pangenome <- c("pangen", "pannagram", "pangenome")
  gr.accs.e <- "accs/"
  gr.accs.b <- "/accs"
  
  # --- Check input arguments ---
  if (p.beg > p.end) stop("Argument p.beg must be less than or equal to p.end")
  if (p.beg <= 0) stop("Argument p.beg must be positive")
  if (!mode %in% c("seq", "pos")) stop("Argument mode must be either 'seq' or 'pos'")
  
  # --- Determine Pannagram version and corresponding paths ---
  dot.args <- list(...)
  if ("path.msa" %in% names(dot.args)) {
    pokazAttention("You are working with output of Pannagram v1.1")
    if (!is.null(path.proj)) stop("Path to the project folder should be provided only for Pannagram v2.X")
    path.msa <- dot.args$path.msa
    path.seq <- file.path(path.msa, "seq")
    pannagram.version <- 1
  } else {
    pannagram.version <- 2
    if (is.null(path.proj)) stop("Path to the project folder must be provided!")
    path.msa <- file.path(path.proj, "features", "msa")
    path.seq <- file.path(path.proj, "features", "seq")
  }
  
  # --- Construct file suffix and MSA file path ---
  ref.suff <- if (ref.acc == '') '' else paste0('_', ref.acc)
  file.msa <- file.path(path.msa, paste0(aln.type, i.chr, '_', i.chr, ref.suff, '.h5'))
  
  if (!file.exists(file.msa)) stop(paste("File", file.msa, "does not exist"))
  
  # --- Extract accessions from MSA file ---
  groups <- rhdf5::h5ls(file.msa)
  accessions <- groups$name[groups$group == gr.accs.b]
  
  # --- Determine file path depending on mode ---
  if (mode == "seq") {
    if (!dir.exists(path.seq)) stop("Please run script 'features' with flag -seq.")
    file.mode <- file.path(path.seq, paste0("seq_", i.chr, "_", i.chr, ref.suff, ".h5"))
  } else if (mode == "pos") {
    file.mode <- file.msa
  }
  
  if (!file.exists(file.mode)) stop(paste("File", file.mode, "does not exist"))
  
  # --- Map genomic positions to alignment indices ---
  if (echo) pokaz("Define new pos based on the accession", acc)
  
  if (tolower(acc) %in% s.pangenome) {
    info <- rhdf5::h5ls(file.msa)
    info <- info[info$group == gr.accs.b, ]
    v <- 1:as.numeric(info$dim[1])
  } else {
    if (!(acc %in% accessions)) stop("Provided acc is not in the alignment")
    v <- rhdf5::h5read(file.msa, paste0(gr.accs.e, acc))
  }
  
  p.beg.acc <- which(v == p.beg)
  p.end.acc <- which(v == p.end)
  
  if (length(p.beg.acc) == 0) stop(paste("Position", p.beg, "is not found in the alignment of the accession", acc))
  if (length(p.end.acc) == 0) stop(paste("Position", p.end, "is not found in the alignment of the accession", acc))
  
  v <- v[p.beg.acc:p.end.acc]
  v <- v[v != 0]
  if (is.unsorted(v)) stop("The region is not in one synteny block")
  
  p.beg <- p.beg.acc
  p.end <- p.end.acc
  
  # --- Read aligned sequences for all accessions ---
  groups <- rhdf5::h5ls(file.mode)
  accessions <- groups$name[groups$group == gr.accs.b]
  
  aln.mx = c()
  for(acc in accessions){
    if(echo) pokaz('Sequence of accession', acc)
    v = rhdf5::h5read(file.mode, paste0(gr.accs.e, acc))
    aln.mx = rbind(aln.mx, v[p.beg:p.end])
  }
  rownames(aln.mx) <- accessions
  return(aln.mx)
}
