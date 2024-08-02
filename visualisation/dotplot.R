# Dotplot

#' Create a Dotplot for Two Nucleotide Sequences
#'
#' @description
#' `dotplot` generates a dotplot to visualize the relationship between two nucleotide sequences.
#' It compares two nucleotide sequences, represented as character vectors, using
#' a specified window size and a minimum number of matches. The function produces a graphical
#' representation with colored tiles, showing areas of similarity and differences
#' between the two sequences, including both forward (dark) and reverse complement (pink) alignments.
#' The color intensity of each tile reflects the number of matches within the specified window size.
#'
#' @param seq1 A character vector representing the first nucleotide sequence.
#' @param seq2 A character vector representing the second nucleotide sequence.
#' @param wsize Window size for the comparison (an integer). This specifies the length of
#' the sequence segment to be compared at each position.
#' @param nmatch Minimum number of matches required within the window (an integer). This
#' defines the threshold for considering a match significant.
#'
#' @return A ggplot object representing the dotplot. 
#'
#' @examples
#' # Example usage:
#' seq1 <- c("A", "C", "G", "T", "A", "C", "G", "T", "A", "C", "G", "T")
#' seq2 <- c("T", "G", "C", "A", "T", "G", "C", "A", "T", "G", "C", "A")
#' p <- dotplot(seq1, seq2, wsize = 3, nmatch = 2)
#' p
#'
#' @export
#'
dotplot <- function(seq1, seq2, wsize, nmatch) {
  
  if(wsize < nmatch) stop('wsize must be larger than nmatch')
  seq2.rc = revCompl(seq2)
  
  mx1 = toupper(seq2mx(seq1, wsize))
  mx2 = toupper(seq2mx(seq2, wsize))
  
  result = mxComp(mx1, mx2, wsize, nmatch)
  
  mx2.rc = toupper(seq2mx(seq2.rc, wsize))
  result.rc = mxComp(mx1, mx2.rc, wsize, nmatch)
  result.rc$values = -result.rc$values
  result.rc$col = length(seq2) - result.rc$col - wsize + 2
  result = rbind(result.rc, result)
  result = rbind(result, data.frame(row=1, col=length(seq2)-1, values=0))
  
  len1 = length(seq1)
  len2 = length(seq2)
  
  p = invisible(
        ggplot(result, aes(x = row, y = col, fill = values, color = values)) +
        geom_tile(width = 1, height = 1, linewidth = 0.5) +
        xlab(NULL) + ylab(NULL) +
        # xlim(c(0, len1)) +
        # ylim(c(0, len2)) +
        theme_minimal() + coord_fixed() +
        scale_x_continuous(expand = c(0, 0), limits = c(0, length(seq1))) + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, length(seq2))) +
        scale_fill_gradient2(low = "#CE1F6A", mid = "white", high = "#27374D",
                             breaks = c(-wsize, 0, wsize)) +
        scale_color_gradient2(low = "#CE1F6A", mid = "white", high = "#27374D",
                              breaks = c(-wsize, 0, wsize)) +
        theme(panel.border = element_rect(colour = "grey", fill = NA, size = 1)) +
        guides(fill = FALSE, color = FALSE) + 
        annotate('text', x = len1/2, y = len2, 
                 label = paste('(',wsize,',',nmatch, ')', sep = ''), 
                 vjust = 1.2, hjust = 0.5)
    )
    
  # p
  return(p )
}


#'  Create a Dotplot for Two Nucleotide Sequences
#'
#' @description
#' The same as `dotplot` but the sequences can be provided as strings
#' 
dotplot.s <- function(seq1, seq2, wsize, nmatch) {
  return(dotplot(seq2nt(seq1), seq2nt(seq2), wsize, nmatch))
}

#' Dotplot for Sequences and Their Reverse Complements
#'
#' @description
#' This function generates two dotplots: one comparing a nucleotide sequence to itself,
#' and another comparing the sequence to its reverse. 
#' 
#' @param seq A character vector representing the nucleotide sequence.
#' @param wsize Window size for the comparison (an integer).
#' @param nmatch Minimum number of matches within the window (an integer).
#'
#' @return A combined ggplot object displaying both dotplots side by side.
#'
#' @examples
#' # Example usage:
#' seq <- c("A", "C", "G", "T", "A", "C", "G", "T", "A", "C", "G", "T")
#' dotspiegel(seq, wsize = 3, nmatch = 2)
#'
#' @export
dotspiegel <- function(seq, wsize, nmatch) {
  
  p.fw = dotplot(seq, seq, wsize, nmatch) + ggtitle('seq vs seq')
  p.rev = dotplot(seq, rev(seq), wsize, nmatch) + ggtitle('seq vs rev(seq)')
  
  pp = invisible(grid.arrange(p.fw, p.rev, nrow = 1))
  
  return(pp)
}

#' Dotplot for Sequences and Their Reverse Complements
#'
#' @description
#' The same as `dotspiegel` but the sequence can be provided as a string
#'
#' @param seq Sequence as a string.
#' @param wsize Comparison window size.
#' @param nmatch Minimum matches in window.
#'
#' @return ggplot of sequence self and reverse comparison.
#' 
#' @export
dotspiegel.s <- function(seq, wsize, nmatch) {
  return(dotspiegel(seq2nt(seq), wsize, nmatch))
}

#' Several dotplots with varying nmatch
#'
#' @param seq1 The first sequence to be compared.
#' @param seq2 The second sequence to be compared.
#' @param wsize The window size to be used in the dot plot.
#' @param nmatch.beg Optional; the beginning number of matches required. 
#'                   Defaults to the window size if not specified.
#' @param nmatch.end Optional; the ending number of matches required.
#'                   Defaults to the window size minus 5 if not specified.
#' @param n.row Optional; the number of rows to arrange the plots in the grid.
#'              Defaults to 2.
#' @return A graphical object containing the grid of dotplots.
#' 
#' @export
dotfacet <- function(seq1, seq2, wsize, nmatch.beg=NULL, nmatch.end=NULL, n.row = 2){
  if(is.null(nmatch.beg)){
    nmatch.beg = wsize
    nmatch.end = NULL
  }
  
  if(is.null(nmatch.end)){
    nmatch.end = wsize - 5
  }
  
  nmatch.range = nmatch.beg:nmatch.end
  # pokaz('Range:', nmatch.range)
  
  p.list = list()
  for(nmatch in nmatch.range){
    p.list[[length(p.list) + 1]] = dotplot(seq1, seq2, wsize, nmatch)
  }
  
  pp = cowplot::plot_grid(plotlist = p.list, nrow=n.row)
  return(pp)
  
}



dotfacet.s <- function(seq1, seq2, ...) {
  return(dotfacet(seq2nt(seq1), seq2nt(seq2),...))
}


#' Compare Two Matrices for Dotplot Generation
#'
#' @description
#' `mxComp` compares two matrices and identifies matches based on a specified window size 
#' and minimum number of matches. It is used for generating a dotplot by comparing 
#' nucleotide sequences represented in matrix form.
#'
#' @param mx1 The first matrix representing a nucleotide sequence.
#' @param mx2 The second matrix representing a nucleotide sequence.
#' @param wsize Window size for the comparison. This parameter is not directly used in the function, 
#' but it's kept for consistency with the overall dotplot generation process.
#' @param nmatch Minimum number of matches required for a significant comparison result.
#'
#' @return A data frame with indices of matching positions and the number of matches at these positions.
#' The data frame has three columns: row index, column index, and the value (number of matches).
#'
mxComp <- function(mx1, mx2, wsize, nmatch){
  mx.res = 0
  for(s in c('A', 'C', 'G', 'T')){
    mx.res = mx.res + (mx1 == s) %*% t(mx2 == s)
  }
  # mx.res = (mx.res >= nmatch) * 1
  mx.res[mx.res < nmatch] = 0
  
  indices <- which(mx.res != 0, arr.ind = TRUE)
  values <- mx.res[indices]
  result <- cbind(indices, values)
  result = as.data.frame(result)
  return(result)
}


#' Calculate Sequence Complexity
#'
#' @description
#' `seqComplexity` calculates the complexity of a nucleotide sequence. The complexity is 
#' determined based on the number of matches found in a dotplot comparison of the sequence 
#' against itself and its reverse complement. The method used for comparison can be specified.
#'
#' @param seq1 A character vector representing the nucleotide sequence.
#' @param method The method used for complexity calculation, default is 'dotplot'.
#' @param wsize Window size for the comparison, default is 10.
#' @param nmatch Minimum number of matches required for considering a segment as complex, default is 9.
#'
#' @return A numeric value representing the complexity of the sequence. This value is 
#' calculated as the total number of matching segments normalized by the length of the sequence.
#'
seqComplexity <- function(seq1, method='dotplot', wsize=10, nmatch=9) {
  
  mx1 = toupper(seq2mx(seq1, wsize))
  result = mxComp(mx1, mx1, wsize, nmatch)
  
  seq1.rc = revCompl(seq1)
  mx1.rc = toupper(seq2mx(seq1.rc, wsize))
  result.rc = mxComp(mx1, mx1.rc, wsize, nmatch)
  
  n.match = (nrow(result) + nrow(result.rc)) / length(seq1) 
  
  return(n.match)
}


dotScore <- function(seq1, seq2, method='dotplot', wsize=10, nmatch=9) {
  
  mx1 = toupper(seq2mx(seq1, wsize))
  mx2 = toupper(seq2mx(seq2, wsize))
  result = mxComp(mx1, mx2, wsize, nmatch)
  n.forward = nrow(result)
  
  seq2.rc = revCompl(seq2)
  mx2.rc = toupper(seq2mx(seq2.rc, wsize))
  
  result.rc = mxComp(mx1, mx2.rc, wsize, nmatch)
  n.backward = nrow(result.rc)
  
  n = (length(seq1) + length(seq2)) / 2
  return(c(n.forward, n.backward) / n)
}

