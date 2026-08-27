# Interval-format readers (_n). The shared helper `getPannagramPaths()` is
# defined in analys/cut_func.R and reused here.

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
#' @note INTERVAL FORMAT (_n). Reads alignments written by `pannagram_n`, where
#'   /accs/<acc> is an interval table (see utils/interval_func.R). In `mode = "pos"`
#'   only the requested window is materialised -- the rest of the alignment is never
#'   touched -- and the accession position is mapped by binary search instead of the
#'   former `which(v == pos)` scan over the whole pangenome.
#' @export
getRegion_n <- function(i.chr, acc, p.beg, p.end,
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
    # Interval tables are chunked by construction, so there is no separate
    # "_chunked" variant to prefer any more.
    file.mode <- file.msa
  }
  
  if (!file.exists(file.mode)) stop(paste("File", file.mode, "does not exist"))
  
  pokaz(file.msa)
  pokaz(file.mode)
  
  # --- Map genomic positions to alignment indices ---
  if (echo) pokaz("Define new pos based on the accession", acc)
  
  if (!(tolower(acc) %in% s.pangenome)) {
    
    if (!(acc %in% accessions)) stop("Provided acc is not in the alignment")
    ivl.ref <- h5IvlRead(file.msa, acc)
    len.pan <- h5AccLen(file.msa, acc)

    # Binary search instead of `which(v == pos)`. The legacy code compared the SIGNED
    # stored value, i.e. it only ever matched plus-strand positions -- keep that.
    mapPos <- function(q){
      i <- ivlPanAt(ivl.ref, q)
      i <- i[i > 0]
      if (length(i) == 0) return(integer(0))
      i[ivlToVecRange(ivl.ref, i[1], i[1]) == q]
    }
    p.beg <- mapPos(p.beg)
    p.end <- mapPos(p.end)
    
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
    
    v <- ivlToVecRange(ivl.ref, p.beg, p.end)
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
    if (mode == "pos") {
      # Only the window is expanded; nothing else of the accession is read.
      aln.mx[i, ] <- ivlToVecRange(h5IvlRead(file.mode, acc_i), p.beg, p.end)
    } else {
      aln.mx[i, ] <- rhdf5::h5read(
        file.mode,
        paste0(gr.accs.e, acc_i),
        index = list(p.beg:p.end)
      )
    }
  }
  return(aln.mx)
}


cutAln_n <- function(...) {
  pokazAttention("Function 'cutAln_n()' is deprecated. Please use 'getRegion_n()' instead.")
  getRegion_n(...)
}

