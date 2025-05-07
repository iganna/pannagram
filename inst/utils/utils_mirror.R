#' Validate and Normalize a Replacement Mask
#'
#' This function checks whether a replacement mask is valid.
#' It ensures that all names and values are uppercase and that the values are unique.
#'
#' @param mask A named character vector where each name is a symbol to be replaced,
#'   and the corresponding value is the symbol it should be replaced with.
#'
#' @return A normalized named character vector with uppercase names and values.
#'
# @export
checkMask <- function(mask){
  names(mask) <- toupper(names(mask))
  mask <- toupper(mask)
  if (sum(duplicated(mask)) != 0) stop("The mask values must be unique")
  return(mask)
}


#' Replace a Sequence of Letters According to a Mask
#'
#' Applies a replacement mask to a sequence of characters.
#'
#' @param seq A character vector representing the input sequence.
#' @param mask A named character vector where each name is a symbol to be replaced,
#'   and the corresponding value is the symbol it should be replaced with.
#'
#' @return A character vector with characters replaced according to the mask.
#'
# @export
maskVect <- function(seq, mask) {
  seq <- toupper(seq)
  mask <- checkMask(mask)
  
  invalid <- setdiff(seq, names(mask))
  if (length(invalid) > 0) {
    stop(paste("Invalid symbols found in input sequence:", paste(invalid, collapse = ", ")))
  }
  
  seq.msk <- mask[seq]
  return(seq.msk)
}


#' Invert a Replacement Mask
#'
#' Given a valid replacement mask, returns its inverse,
#' where the values become the new keys and the keys become the new values.
#'
#' @param mask A named character vector where each name is a symbol to be replaced,
#'   and the corresponding value is the replacement symbol.
#'
#' @return A named character vector where names and values of the original mask are swapped.
#'
# @export
invertMask <- function(mask) {
  mask <- checkMask(mask)
  inv <- names(mask)
  names(inv) <- mask
  return(inv)
}


#' Apply a Replacement Mask to a Set of Sequences
#'
#' @param seqs A named character vector of sequences.
#' @param mask A named character vector specifying the replacement mask. Each name is a symbol to be replaced,
#'   and each value is the symbol to replace it with.
#'
#' @return A named character vector of masked sequences.
#'
#' @examples
#' maskSeqs(c(seq1 = "ATGC"), c(A = "X", T = "Y", G = "Z", C = "W"))
#'
maskSeqs <- function(seqs, mask) {
  seqs.masked <- sapply(seqs, function(s) {
    s <- seq2nt(s)
    s <- maskVect(s, mask)
    nt2seq(s)
  }, USE.NAMES = TRUE)
  
  return(seqs.masked)
}


genMask <- function(bases = c("A", "C", "G", "T"), echo=T){
  permute <- function(x) {
    if (length(x) == 1) return(list(x))
    out <- list()
    for (i in seq_along(x)) {
      sub <- x[-i]
      for (p in permute(sub)) {
        out <- append(out, list(c(x[i], p)))
      }
    }
    return(out)
  }
  
  # Generate permutations
  perms_list <- permute(bases)
  
  # To matrix
  mx <- do.call(rbind, perms_list)
  
  colnames(mx) = bases
  
  # Print
  if(echo){
    for(irow in 1:nrow(mx)){
      pokaz(sprintf("%2d:", irow-1), mx[irow,])
    }  
  }
  
  # Remove the identical mask
  mx = mx[-1,]
  
  return(mx)
}


#' Get Mask 
#'
#' @param id Integer. The row index to extract from the mask matrix.
#' @param bases Character vector. Nucleotide bases to include in the mask.
#'
#' @return A named vector representing the selected mask row.
#'
#' @examples
#' getMask()
#' getMask(id = 2, bases = c("A", "T"))
#'
# @export
getMask <- function(id = 9, bases = c("A", "C", "G", "T")){
  mx = genMask(bases, echo=F)
  mask = setNames(mx[id,], colnames(mx))
  return(mask)
}






