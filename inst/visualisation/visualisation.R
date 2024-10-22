library(ggplot2)


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


#' @export
plotSynDot <- function(..., alpha = 1) {
  pokazAttention('Please use plotSynteny(..., show.dot=T)')
  p = plotSynteny(...) + geom_point(show.legend = FALSE, size = 1, alpha = alpha)
  return(p)
}


#' Plot Synteny Blocks Between Two Genomic Regions (Alias for \code{\link{plotSynteny}})
#'
#' This function is an alias for the \code{\link{plotSynteny}} function and 
#' generates a ggplot2-based visualization of synteny in the alignment. 
#' 
#' @export
plotSyntenyBlocks <- function(...) {
  plotSynteny(...)
}

#' @export
plotPanAcc <- function(file.msa, acc){
  
  # Setup
  gr.accs.e <- "accs/"
  axis.breaks = seq(0, 30, by = 5)
  idx.step = 10000
  idx = seq(1, pan.len, idx.step)
  c.red = '#872341'
  c.grey = '#EBEBEB'
  
  # Read the correpondence for one accession
  v.acc = h5read(file.msa, paste0(gr.accs.e, acc))
  v.acc = data.frame(pan = idx, acc = v.acc[idx])
  v.acc = v.acc[v.acc$acc != 0,]
  
  # Generate a ggplot for each accession
  p = ggplot(v.acc, aes(x = pan, y = abs(acc), color=as.factor(sign(acc)))) + 
    geom_abline(slope = 1, intercept = 0, color = c.grey) +
    geom_point(size = 0.01) + 
    theme_minimal() +
    labs(x = "Pangenome coord, Mbp", 
         y = paste0(acc, ', Mbp'), 
         title = NULL) +
    scale_color_manual(values = c("-1" = c.red, "1" = "black")) +
    theme(legend.position = "none") + 
    scale_y_continuous(breaks = axis.breaks * 10^6, labels = axis.breaks) + 
    scale_x_continuous(breaks = axis.breaks * 10^6, labels = axis.breaks)
  
  return(p)
  
}


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
plotSynAllChr <- function(path.aln, 
                               acc, 
                               ref,
                               chr.len,
                               seq.order = "default"
                               ) {
  
  # # Get alignments and numbers of chromosomes
  # pattern <- ".*_[0-9]+_[0-9]+_maj\\.rds$"
  # files.aln <- list.files(path = path.aln, pattern = pattern, full.names = F)
  # acc.ids <- unique(sapply(files.aln, function(s) sub("^(.*?)_\\d+_\\d+_maj\\.rds$", "\\1", s)))
  # acc.chrs <- unique(sapply(files.aln, function(s) sub(".*_([0-9]+)_[0-9]+_maj\\.rds$", "\\1", s)))
  # ref.chrs <- unique(sapply(files.aln, function(s) sub(".*_[0-9]+_([0-9]+)_maj\\.rds$", "\\1", s)))
  
  if(!(ref %in% chr.len$acc)) stop('Number of chromosomes in the reference genomes in not defined')
  if(!(acc %in% chr.len$acc)) stop('Number of chromosomes in the accession genomes in not defined')
  
  # Number of chromosomes
  n.chr.ref = max(chr.len$chr[chr.len$acc == ref])
  n.chr.acc = max(chr.len$chr[chr.len$acc == acc])
  
  # Lengths of chromosomes
  chr.len.acc = chr.len[chr.len$acc == acc,]
  chr.len.ref = chr.len[chr.len$acc == ref,]
  chr.len.acc = chr.len.acc[order(chr.len.acc$chr),]
  chr.len.ref = chr.len.ref[order(chr.len.ref$chr),]
  
  # Check that all choromosomes are presented
  if(nrow(chr.len.acc) != n.chr.acc) stop('Wrong number of chromosomes in accession')
  if(nrow(chr.len.ref) != n.chr.ref) stop('Wrong number of chromosomes in reference')
  
  
  # === === === === Reordering === === === ===
  if (seq.order=="descending") {
    order.acc <- order(-n.chr.acc)
    order.ref <- order(-n.chr.ref)
  } else if (seq.order == "alphanum"){
    stop('Wrong name')
    # order.acc <- order(acc.labels)
    # order.ref <- order(ref.labels)
  } else if (seq.order == "default"){
    order.acc <- 1:n.chr.acc
    order.ref <- 1:n.chr.ref
  } else {
    stop("Unknown keyword for `seq.order`. Options: 'default', 'descending'")
  }
  
  # === === === === Main === === === ===
  
  # Cummulative lengths
  cum.acc <- c(0, cumsum(chr.len.acc$len[order.acc]))
  cum.ref <- c(0, cumsum(chr.len.ref$len[order.ref]))
  
  # pokaz('Number of chromosomes ref and acc:', n.chr.ref, n.chr.acc)
  
  # Read the alignments
  df <- data.frame()
  for (i.acc in order.acc) {
    for (i.ref in order.ref) {
      file.aln = paste0(path.aln, acc, "_", i.acc, "_", i.ref, "_maj.rds")
      pokaz(file.aln)
      if (file.exists(file.aln)) {
        pokaz('have it')
        # Synteny
        data.ij <- readRDS(file.aln)
        data.ij[, c(2, 3)] = data.ij[, c(2, 3)] + cum.acc[i.acc]
        data.ij[, c(4, 5)] = data.ij[, c(4, 5)] + cum.ref[i.ref]
        df <- rbind(df, data.ij)
      }
    }
  }
  
  v.sep.acc = cum.acc[-length(cum.acc)]
  h.sep.ref = cum.ref[-length(cum.ref)]
  
  # Annotations
  
  if("name" %in% colnames(chr.len)){
    acc.labels = chr.len.acc$name[order.acc]
    ref.labels = chr.len.ref$name[order.ref]
  } else {
    acc.labels = paste('Chr', order.acc)
    ref.labels = paste('Chr', order.ref)  
  }
  
  
  annot.acc = data.frame(x = cum.acc[-1], 
                           y = rep(0, length(cum.acc) - 1),
                           label = acc.labels)
  annot.ref = data.frame(x = rep(0, length(cum.ref) - 1), 
                           y = cum.ref[-1],
                           label = ref.labels)
  
  # # Remain only the existing chromosomes
  # v.sep.acc = v.sep.acc[1:max.acc.chr]
  # h.sep.ref = h.sep.ref[1:max.ref.chr]
  # 
  # annot.acc = annot.acc[1:max.acc.chr,]
  # annot.ref = annot.ref[1:max.ref.chr,]
  # 
  # 
  
  # Plot
  synteny.plot = plotSynteny(df,
                             query.label = acc,
                             ref.label = ref,
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


