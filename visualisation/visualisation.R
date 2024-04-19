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
#' @param show.point Flag to show the plot thicker with the help of dots in the beginning of synteny segments.
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
                        col.line = '#362FD9',
                        show.point = F){
  
  if(!is.null(base.len)){
    x = getBase(x, base.len)
  }
  
  range.x = max(c(x$V2,x$V3)) - min(c(x$V2,x$V3))
  range.y = max(c(x$V4,x$V5)) - min(c(x$V4,x$V5))
  range.m = max(range.x, range.y)
  max.lim = max(c(x$V4,x$V5,x$V2,x$V3))
  
  n.scale.step = 5
  for(i.scale in 6:2){
    n.scale = 10^i.scale 
    if((range.m > n.scale* n.scale.step) || (i.scale == 2)){
      seq.lab.n = round(max.lim / n.scale) + 2
      seq.lab =  seq(0, seq.lab.n, by = n.scale.step)
      v.ticks = seq.lab * n.scale
      s.ticks = paste(seq.lab, 'e', i.scale ,sep = '' )
      break
    }
  }
  # pokaz(v.ticks)
  # pokaz(s.ticks)
  
  p <- ggplot(x, aes(x = V2, y=V4, xend = V3, yend = V5, color = as.factor(V4 < V5))) + 
    # geom_point(show.legend = FALSE, size = 0.1) + 
    geom_segment(show.legend = FALSE) + 
    theme_bw() + 
    xlab('query') + 
    ylab('base') +
    scale_color_manual(values = c("FALSE" = col.rc, "TRUE" = col.fw)) +
    coord_fixed(ratio = 1) +
    scale_y_continuous(breaks = v.ticks, labels = s.ticks) +
    scale_x_continuous(breaks = v.ticks, labels = s.ticks)
  
  
  
  if(!is.null(hlines)){
    p <- p + geom_hline(yintercept=hlines, color= col.line)
  }
  
  if(!is.null(vlines)){
    p <- p + geom_vline(xintercept=vlines, color = col.line)
  }
  
  if(show.point){
    p <- p + geom_point(show.legend = FALSE, size = 1, alpha = alpha)
  }
  
  return(p)
}


ploDot <- function(..., alpha = 1) {
  pokazAttention('Please use plotSynteny(..., show.dot=T)')
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

#' ----------------------------------------------------------------------
plotPanAcc <- function(file.msa, acc){
  
  # Setup
  gr.accs.e <- "accs/"
  axis.breaks = seq(0, 30, by = 5)
  idx.step = 10000
  idx = seq(1, pan.len, idx.step)
  c.red = '#872341'
  c.grey = '#EBEBEB'
  
  # Read the correpondence for one accession
  v.acc = h5read(file.msa, paste(gr.accs.e, acc, sep = ''))
  v.acc = data.frame(pan = idx, acc = v.acc[idx])
  v.acc = v.acc[v.acc$acc != 0,]
  
  # Generate a ggplot for each accession
  p = ggplot(v.acc, aes(x = pan, y = abs(acc), color=as.factor(sign(acc)))) + 
    geom_abline(slope = 1, intercept = 0, color = c.grey) +
    geom_point(size = 0.01) + 
    theme_minimal() +
    labs(x = "Pangenome coord, Mbp", 
         y = paste(acc, ', Mbp', sep = ''), 
         title = NULL) +
    scale_color_manual(values = c("-1" = c.red, "1" = "black")) +
    theme(legend.position = "none") + 
    scale_y_continuous(breaks = axis.breaks * 10^6, labels = axis.breaks) + 
    scale_x_continuous(breaks = axis.breaks * 10^6, labels = axis.breaks)
  
  return(p)
  
}


