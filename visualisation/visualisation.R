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
#' @param axis.ticks Numeric vector for axis ticks to denote Mbases of sequence lengths. Default is `seq(0, 10, by = 5)`
#' @param expand.axis A vector of range expansion constants used to add axis padding. Defaults to `waiver()`. See docs for `scale_x_continuous` or `scale_y_continuous` for more details
#'
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
                        axis.ticks = seq(0, 10, by = 5),
                        expand.axis = waiver()
){
  if(!is.null(base.len)) x = getBase(x, base.len)
  if (is.null(query.label)) query.label = 'query'
  if (is.null(ref.label)) ref.label = 'base'
  
  p <- ggplot(x, aes(x=V2, y=V4, xend=V3, yend=V5, color=as.factor(V4 < V5))) + 
    geom_segment(show.legend = FALSE) + 
    theme_bw() + 
    xlab(query.label) + 
    ylab(ref.label) +
    scale_color_manual(values = c("FALSE"=col.rc, "TRUE"=col.fw)) +
    coord_fixed(ratio=1) +
    scale_y_continuous(breaks=axis.ticks*1e6, labels=axis.ticks, expand=expand.axis) +
    scale_x_continuous(breaks=axis.ticks*1e6, labels=axis.ticks, expand=expand.axis)
  
  if(!is.null(hlines)) p <- p + geom_hline(yintercept=hlines, color= col.line)
  if(!is.null(vlines)) p <- p + geom_vline(xintercept=vlines, color= col.line)
  if(show.point)       p <- p + geom_point(show.legend = FALSE, size = 0.8, alpha = point.alpha)
  
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

#' ----------------------------------------------------------------------
#' Plot Synteny Between Two Genomes
#'
#' This function wraps `plotSynteny` to visualize synteny regions between two full genomes, reading the files, produced by the pangen pipeline
#'
#' @param alignments.path Char with path to **directory**, where `.rds` files with alignments are located. By default, pangen pipeline writes such files into dir with **alignments_** prefix
#' @param query.name Char with path to **FASTA file** with genome, used as query (horizontal axis)
#' @param ref.name Char with path to **FASTA file** with genome, used as reference (vertical axis)
#' @param sort.descending Flag to sort sequences (chromosome, plasmids, fragments, etc.) of a genome by their lengths in descending order. Default: `TRUE`
#' @param query.label Char for x axis label (query). Default `NULL` will automatically deduce label based on given `query.name`
#' @param ref.label Char for y axis label (reference). Default `NULL` will automatically deduce label based on given `ref.name`
#'
#' @return A `ggplot2` plot object.
#'
#' @import ggplot2
#' @export
plotGenomeAgainstRef <- function(alignments.path, query.name, ref.name,
                                 sort.descending=T,
                                 query.label = NULL,
                                 ref.label = NULL) {
  # if none given, deducing axis labels straight from filenames
  if (is.null(query.label)) query.label = tools::file_path_sans_ext(basename(query.name))
  if (is.null(ref.label)) ref.label = tools::file_path_sans_ext(basename(ref.name))
  
  # Read FASTA files
  fasta_query <- readFastaMy(query.name)
  fasta_ref <- readFastaMy(ref.name)
  
  # =============== reordering ===================
  if (sort.descending){
    order_query <- order(-nchar(fasta_query))
    order_ref <- order(-nchar(fasta_ref))
    fasta_query <- fasta_query[order_query]
    fasta_ref <- fasta_ref[order_ref]
  } else {
    order_query <- seq(1, length(nchar(fasta_query)))
    order_ref <- seq(1, length(nchar(fasta_ref)))
  }
  
  len_query = nchar(fasta_query)
  len_ref = nchar(fasta_ref)
  cum_query = c(0, cumsum(len_query))
  cum_ref = c(0, cumsum(len_ref))
  
  
  #=========== Filling the new dataframe =====================
  df <- data.frame()
  query_prefix <- tools::file_path_sans_ext(basename(query.name))
  for (i in seq_along(fasta_query)) {
    for (j in seq_along(fasta_ref)) {
      
      file_name = paste(query_prefix, "_", order_query[i], "_", order_ref[j], "_maj.rds", sep='')
      file_path = file.path(alignments.path, file_name)
      if (file.exists(file_path)) {
        data_ij <- readRDS(file_path)
        
        data_ij[, c(2, 3)] = data_ij[, c(2, 3)] + cum_query[i]
        data_ij[, c(4, 5)] = data_ij[, c(4, 5)] + cum_ref[j]
        df = rbind(df, data_ij)
      }
    }
  }
  
  v.plasmid.divisors = cum_query[-length(cum_query)]
  h.plasmid.divisors = cum_ref[-length(cum_ref)]
  
  query_labels = sub("^[^ ]+ ", "", names(fasta_query))
  ref_labels = sub("^[^ ]+ ", "", names(fasta_ref))
  
  pS = plotSynteny(df,
                   query.label=query.label,
                   ref.label=ref.label,
                   axis.ticks=seq(0, 10, by = 2),
                   hlines=h.plasmid.divisors,
                   vlines=v.plasmid.divisors,
                   col.line = "#3530D966",
                   show.point = T,
                   expand = c(0,0)
  ) +
    annotate("text", x = cum_query[-1], y = rep(0, length(cum_query) - 1),
             label = query_labels,
             vjust = -0.3, hjust = -0.02,
             size = 2.7, angle=90,
             color="#7D7D7D") +
    annotate("text", x = rep(0, length(cum_ref) - 1), y = cum_ref[-1], 
             label = ref_labels,
             vjust = 1.2, hjust = -0.02,
             size = 2.7,
             color="#7D7D7D")
  return(pS)
}
