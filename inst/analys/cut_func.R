#' Determine Paths to MSA and Sequence Folders
#'
#' This function determines the paths to the MSA and sequence directories
#' based on whether the input is from Pannagram v1.1 or v2.X.
#'
#' @param path.proj A character string specifying the path to the project folder (for Pannagram v2.X).
#'        Required unless using output from Pannagram v1.1.
#' @param dot.args A list of additional arguments. If `path.msa` is present,
#'        it indicates usage of output from Pannagram v1.1.
#'
#' @return A list with two character elements: `path.msa` and `path.seq`, indicating the
#'         paths to the MSA and sequence directories, respectively.
getPannagramPaths <- function(path.proj = NULL, dot.args = list()) {
  if ("path.msa" %in% names(dot.args)) {
    pokazAttention("You are working with output of Pannagram v1.1")
    if (!is.null(path.proj)) stop("Path to the project folder should be provided only for Pannagram v2.X")
    path.msa <- dot.args$path.msa
    path.seq <- file.path(path.msa, "seq")
  } else if ("path.cons" %in% names(dot.args)) {
    pokazAttention("You are working with output of Pannagram v1.1")
    if (!is.null(path.proj)) stop("Path to the project folder should be provided only for Pannagram v2.X")
    path.msa <- dot.args$path.cons
    path.seq <- file.path(path.msa, "seq")
  } else {
    if (is.null(path.proj)) stop("Path to the project folder must be provided!")
    path.msa <- file.path(path.proj, "features", "alignments")
    path.seq <- file.path(path.proj, "features", "consensus")
  }
  path.msa = paste0(path.msa, '/')
  path.seq = paste0(path.seq, '/')
  return(list(path.msa = path.msa, path.seq = path.seq))
}



#' Extract a subregion from an alignment matrix
#'
#' This function extracts a region from a Pangenome alignment
#' or sequence coordinate matrix for a specific chromosome and accession.
#'
#' @param i.chr Chromosome number.
#' @param acc Accession name (or 'pangenome', 'pannagram' etc. for the pangenome coordinate).
#' @param p.beg Start position in the Accession.
#' @param p.end End position in the Accession.
#' @param path.proj Path to project folder.
#' @param mode Mode of extraction: either `"seq"` (sequence coordinates) or `"pos"` (alignment coordinates).
#' @param aln.type Prefix for alignment files (default `"pan"`).
#' @param ref.acc Reference accession ID if you built in the reference-based alignment.
#' @param echoWhether to print verbose messages (default `FALSE`).
#'
#' @return A matrix where each row corresponds to an accession and columns represent aligned bases in the specified region.
#'
#' @export
getRegion <- function(i.chr, acc, p.beg, p.end,
                      path.proj = NULL,
                      mode = 'seq',
                      aln.type = "pan", 
                      ref.acc = '',
                      echo = FALSE, 
                      acc.aln = NULL,...) {
  
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
  pannagram.paths = getPannagramPaths(path.proj, dot.args)
  path.msa = pannagram.paths$path.msa
  path.seq = pannagram.paths$path.seq
  
  # --- Construct file suffix and MSA file path ---
  ref.suff <- if (ref.acc == '') '' else paste0('_', ref.acc)
  if(ref.suff != '') aln.type='ref'
  file.msa <- file.path(path.msa, paste0(aln.type, '_', i.chr, '_', i.chr, ref.suff, '.h5'))  
  
  
  if (!file.exists(file.msa)) stop(paste("File", file.msa, "does not exist"))
  
  # --- Extract accessions from MSA file ---
  if(echo) pokaz("Extract accessions from MSA file")
  groups <- rhdf5::h5ls(file.msa)
  accessions <- groups$name[groups$group == gr.accs.b]
  
  if(!is.null(acc.aln)){
    accessions = intersect(accessions, acc.aln)
    if(length(accessions) == 0) stop('Please provide relevant accession names')
    if(length(accessions) == length(acc.aln)) pokazAttention('Some accession names are not relevant:', 
                                                             setdiff(acc.aln, accessions))
    pokaz('Generating alignment for', length(accessions), 'accessions')
  }
  
  # --- Determine file path depending on mode ---
  if(echo) pokaz("Determine file path depending on mode")
  if (mode == "seq") {
    if (!dir.exists(path.seq)) stop("Please run script 'features' with flag -seq.")
    file.mode <- file.path(path.seq, paste0("seq_", i.chr, "_", i.chr, ref.suff, "_chunked.h5"))
    if(!file.exists(file.mode)){
      file.mode <- file.path(path.seq, paste0("seq_", i.chr, "_", i.chr, ref.suff, ".h5"))
    }
  } else if (mode == "pos") {
    
    file.mode <- file.path(path.msa, paste0(aln.type, '_', i.chr, '_', i.chr, ref.suff, '_chunked.h5'))
    if(!file.exists(file.mode)){
      file.mode <- file.path(path.msa, paste0(aln.type, '_', i.chr, '_', i.chr, ref.suff, '.h5'))  
    }
  }
  
  if (!file.exists(file.mode)) stop(paste("File", file.mode, "does not exist"))
  
  pokaz(file.msa)
  pokaz(file.mode)
  
  # --- Map genomic positions to alignment indices ---
  if (echo) pokaz("Define new pos based on the accession", acc)
  
  if (!(tolower(acc) %in% s.pangenome)) {
    
    if (!(acc %in% accessions)) stop("Provided acc is not in the alignment")
    v <- rhdf5::h5read(file.msa, paste0(gr.accs.e, acc))
    
    p.beg <- which(v == p.beg)
    p.end <- which(v == p.end)
    
    if (length(p.beg) == 0){
      pokazAttention("Position", p.beg, "is not found in the alignment of the accession", acc,
                     '\nReturn empty matrix')
      return(NULL)
    } 
    if (length(p.end) == 0){
      pokazAttention("Position", p.end, "is not found in the alignment of the accession", acc,
                     '\nReturn empty matrix')
      return(NULL)
    }
    
    v <- v[p.beg:p.end]
    v <- v[v != 0]
    if (is.unsorted(v)) stop("The region is not in one synteny block")
    
  }
  
  # --- Read aligned sequences for all accessions ---
  
  n_acc <- length(accessions)
  w <- p.end - p.beg + 1
  aln.mx <- matrix(NA, nrow = n_acc, ncol = w,
                   dimnames = list(accessions, NULL))
  
  for (i in seq_along(accessions)) {
    acc_i <- accessions[i]
    if (echo) pokaz("Sequence of accession", acc_i)
    aln.mx[i, ] <- rhdf5::h5read(
      file.mode,
      paste0(gr.accs.e, acc_i),
      index = list(p.beg:p.end)
    )
  }
  return(aln.mx)
}


cutAln <- function(...) {
  pokazAttention("Function 'cutAln()' is deprecated. Please use 'getRegion()' instead.")
  getRegion(...)
}

