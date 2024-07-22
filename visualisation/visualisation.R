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
#' @param point.alpha Float value for points alpha channel if `show.point` is `TRUE`. Default is `1.0`.
#' @param query.label String to label x axis of the plot
#' @param ref.label String to label y axis of the plot
#' @param expand.axis A vector of range expansion constants used to add axis padding. Defaults to `waiver()`.
#' To eliminate axis offsets set to `c(0,0)`. See docs for `scale_x_continuous` or `scale_y_continuous` for more details.
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
                        show.point = F,
                        point.alpha = 1.0,
                        query.label = NULL,
                        ref.label = NULL,
                        expand.axis = waiver()
){
  if(!is.null(base.len)) x = getBase(x, base.len)
  if (is.null(query.label)) query.label = 'query'
  if (is.null(ref.label)) ref.label = 'base'

  x.limits = c(x$V2,x$V3)
  y.limits = c(x$V4,x$V5)
  
  range.x = max(x.limits) - min(x.limits)
  range.y = max(y.limits) - min(y.limits)
  range.m = max(range.x, range.y)
  max.lim = max(x.limits, y.limits)

  if (range.m <1.05e6){
    n.scale.step <- 0.1
  } else if (range.m <5.1e6){
    n.scale.step <- 0.5
  } else if (range.m <10.2e6){
    n.scale.step <- 1
  } else {
    n.scale.step <- 5
  }
  # n.scale.step = 0.5
  for(i.scale in 6:2){
    n.scale = 10^i.scale
    if((range.m > n.scale*n.scale.step) || (i.scale == 2)){
      seq.lab.n = round(max.lim / n.scale) + 2
      seq.lab =  seq(0, seq.lab.n, by = n.scale.step)
      v.ticks = seq.lab * n.scale
      s.ticks = paste(seq.lab, 'e', i.scale ,sep = '' )
      break
    }
  }
  
  # Replase 0eX to 0
  s.ticks <- ifelse(grepl("^0e", s.ticks), "0", s.ticks)
  
  p <- ggplot(x, aes(x=V2, y=V4, xend=V3, yend=V5, color=as.factor(V4<V5))) +
    geom_segment(show.legend = FALSE) +
    theme_bw() +
    xlab(query.label) +
    ylab(ref.label) +
    scale_color_manual(values = c("FALSE"=col.rc, "TRUE"=col.fw)) +
    coord_fixed(ratio=1) +
    scale_y_continuous(breaks = v.ticks, labels = s.ticks, expand=expand.axis) +
    scale_x_continuous(breaks = v.ticks, labels = s.ticks, expand=expand.axis)

  if(!is.null(hlines)) p <- p + geom_hline(yintercept=hlines, color= col.line)
  if(!is.null(vlines)) p <- p + geom_vline(xintercept=vlines, color= col.line)
  if(show.point)       p <- p + geom_point(show.legend = FALSE, size = 0.8, alpha = point.alpha)
  
  return(p)
}


plotSynDot <- function(..., alpha = 1) {
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

#' ----------------------------------------------------------------------
#' Plot Synteny Between Two Genomes
#'
#' This function wraps `plotSynteny` to visualize synteny regions between two full genomes, reading the files, produced by the pangen pipeline
#'
#' @param alignments.path Char with path to **directory**, where `.rds` files with alignments are located. By default, pangen pipeline writes such files into dir with **alignments_** prefix
#' @param query.name Char with path to **FASTA file** with genome, used as query (horizontal axis)
#' @param ref.name Char with path to **FASTA file** with genome, used as reference (vertical axis)
#' @param seq.order Sequence reordering type (both axis): 'default': fasta file order, 'descending': sorted from longest to shortest sequence, 'alphanum': sorted by sequence labels (excluding accession number)
#' @param query.label Char for x axis label (query). Default `NULL` will automatically deduce label based on given `query.name`
#' @param ref.label Char for y axis label (reference). Default `NULL` will automatically deduce label based on given `ref.name`
#'
#' @return A `ggplot2` plot object.
#'
#' @import ggplot2
#' @export
plotGenomeAgainstRef <- function(alignments.path, acc.name, ref.name,
                                 seq.order = "default",
                                 acc.label = NULL,
                                 ref.label = NULL) {
  # if none given, deducing axis labels straight from filenames
  if (is.null(acc.label)) {
    acc.label <- tools::file_path_sans_ext(basename(acc.name))
  }
  if (is.null(ref.label)) {
    ref.label <- tools::file_path_sans_ext(basename(ref.name))
  }

  # Read FASTA files
  acc.genome <- readFastaMy(acc.name)
  ref.genome <- readFastaMy(ref.name)
  
  # Filter out accession numbers, that are typically present before the space character
  acc.labels <- sub("^[^ ]+ ", "", names(acc.genome))
  ref.labels <- sub("^[^ ]+ ", "", names(ref.genome))
  
  # === === === === Reordering === === === ===
  if (seq.order=="descending") {
    order.acc <- order(-nchar(acc.genome))
    order.ref <- order(-nchar(ref.genome))
  } else if (seq.order == "alphanum"){
    order.acc <- order(acc.labels)
    order.ref <- order(ref.labels)
  } else if (seq.order == "default"){
    order.acc <- seq(1, length(nchar(acc.genome)))
    order.ref <- seq(1, length(nchar(ref.genome)))
  } else {
    stop("Unknown keyword for `seq.order`. Options: 'default', 'alphanum', 'descending'")
  }
  acc.genome <- acc.genome[order.acc]
  acc.labels <- acc.labels[order.acc]
  
  ref.genome <- ref.genome[order.ref]
  ref.labels <- ref.labels[order.ref]
  
  # Cummulative lengths
  len.acc <- nchar(acc.genome)
  len.ref <- nchar(ref.genome)
  cum.acc <- c(0, cumsum(len.acc))
  cum.ref <- c(0, cumsum(len.ref))
  
  # === === === Filling the new data.frame === === ===
  df <- data.frame()

  acc_prefix <- tools::file_path_sans_ext(basename(acc.name))
  max.ref.chr = 0
  max.acc.chr = 0
  for (i in seq_along(acc.genome)) {
    for (j in seq_along(ref.genome)) {
      file_name = paste0(acc_prefix, "_", order.acc[i], "_", order.ref[j], "_maj.rds")
      file_path = file.path(alignments.path, file_name)

      if (file.exists(file_path)) {
        # Max chr numbers
        max.acc.chr = i
        max.ref.chr = j
        
        # Synteny
        data.ij <- readRDS(file_path)
        data.ij[, c(2, 3)] = data.ij[, c(2, 3)] + cum.acc[i]
        data.ij[, c(4, 5)] = data.ij[, c(4, 5)] + cum.ref[j]
        df <- rbind(df, data.ij)
      }
    }
  }
  
  v.sep.acc = cum.acc[-length(cum.acc)]
  h.sep.ref = cum.ref[-length(cum.ref)]
  
  # Annotations
  annot.acc = data.frame(x = cum.acc[-1], 
                           y = rep(0, length(cum.acc) - 1),
                           label = acc.labels)
  annot.ref = data.frame(x = rep(0, length(cum.ref) - 1), 
                           y = cum.ref[-1],
                           label = ref.labels)
  
  # Remain only the existing chromosomes
  v.sep.acc = v.sep.acc[1:max.acc.chr]
  h.sep.ref = h.sep.ref[1:max.ref.chr]
  
  annot.acc = annot.acc[1:max.acc.chr,]
  annot.ref = annot.ref[1:max.ref.chr,]
  
  
  # Plot
  synteny.plot = plotSynteny(df,
                   acc.label = acc.label,
                   ref.label = ref.label,
                   hlines = h.sep.ref,
                   vlines = v.sep.acc,
                   col.line = "#3530D966",
                   show.point = TRUE,
                   expand = c(0, 0)
  ) +
    annotate("text", x = cum.acc[-1], y = rep(0, length(cum.acc) - 1),
             label = acc.labels,
             vjust = -0.3, hjust = -0.02,
             size = 2.7, angle=90,
             color = "#7D7D7D") +
    annotate("text", x = rep(0, length(cum.ref) - 1), y = cum.ref[-1],
             label = ref.labels,
             vjust = 1.2, hjust = -0.02,
             size = 2.7,
             color = "#7D7D7D")
  return(synteny.plot)
}
