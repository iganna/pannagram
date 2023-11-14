library(ggplot2)

#' ----------------------------------------------------------------------
#' Plot Synteny Between Two Genomic Regions
#'
#' This function generates a ggplot2-based visualization of synteny in the alignment.
#' Synteny blocks are displayed as segments.
#'
#' @param x A data frame representing the alignment. Must contain columns `V2`, `V3`, `V4`, and `V5`.
#' (`V2`, `V3`) - begin-end positions in a query sequence
#' (`V4`, `V5`) - begin-end positions in a reference sequence
#' 
#' @param base.len Numeric length of the reference sequence. 
#' If provided, the reference coordinated in `x` will be processed with the `getBase` function.
#' 
#' @param hlines Numeric vector. Y-coordinates at which horizontal lines should be drawn. Default is `NULL`.
#' @param vlines Numeric vector. X-coordinates at which vertical lines should be drawn. Default is `NULL`.
#' @param col.fw Color for forward-oriented synteny blocks. Default is '#27374D'.
#' @param col.rc Color for reverse-complement synteny blocks. Default is '#CE1F6A'.
#' @param col.line Color for the optional horizontal and vertical lines. Default is '#362FD9'.
#'
#' @return A `ggplot2` plot object.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(V2 = c(1, 4), V3 = c(3, 5), V4 = c(1, 2), V5 = c(3, 4))
#' plotSynteny(data)
#' }
#'
#' @import ggplot2
#' @export
plotSynteny <- function(x, base.len = NULL, hlines=NULL, vlines=NULL,
                        col.fw = '#27374D',
                        col.rc = '#CE1F6A',
                        col.line = '#362FD9'){
  
  if(!is.null(base.len)){
    x = getBase(x, base.len)
  }
  p <- ggplot(x, aes(x = V2, y=V4, xend = V3, yend = V5, color = as.factor(V4 < V5))) + 
    # geom_point(show.legend = FALSE, size = 0.1) + 
    geom_segment(show.legend = FALSE) + 
    theme_bw() + 
    xlab('query') + 
    ylab('base') +
    scale_color_manual(values = c("FALSE" = col.rc, "TRUE" = col.fw)) +
    coord_fixed(ratio = 1) 
  
  if(!is.null(hlines)){
    p <- p + geom_hline(yintercept=hlines, color= col.line)
  }
  
  if(!is.null(vlines)){
    p <- p + geom_vline(xintercept=vlines, color = col.line)
  }
  return(p)
}


plotDot <- function(..., alpha = 1) {
  p = plotSynteny(...) + geom_point(show.legend = FALSE, size = 1, alpha = alpha)
  return(p)
}

#' ----------------------------------------------------------------------
#' Plot Synteny Blocks Between Two Genomic Regions (Alias for \code{\link{plotSynteny}})
#'
#' This function is an alias for the \code{\link{plotSynteny}} function and 
#' generates a ggplot2-based visualization of synteny in the alignment. 
#' 
plotSyntenyBlocks <- function(...) {
  plotSynteny(...)
}
