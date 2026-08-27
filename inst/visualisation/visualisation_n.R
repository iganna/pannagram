# ============================================================================
#  Interval-format (_n) counterpart of the visualisation.R reader.
#  Only pangrowth() reads /accs; the rest of visualisation.R is format-independent.
#  See utils/interval_func.R.
# ============================================================================

#' Plot Pangenome Coordinate Projection for One Accession
#'
#' This function visualizes how the coordinates of a single accession map
#' onto a pangenome coordinate for a given chromosome. 
#' 
#' @param path.project Path to the Pannagram project directory.
#' @param acc Accession identifier whose coordinates will
#'   be plotted on the y-axis.
#' @param i.chr Chromosome index.
#'
#' @import ggplot2
#' @export
pangrowth_n <- function(path.project, acc, i.chr, aln.type='pan', ref.acc='', size = 0.5){
  
  ref.suff <- if (ref.acc == '') '' else paste0('_', ref.acc)
  if(ref.suff != '') aln.type='ref'
  
  path.msa = paste0(path.project, '/', 'features/alignments/')
  file.msa = paste0(path.msa, aln.type, '_', i.chr, '_', i.chr, ref.suff, '.h5')
  
  # Setup
  gr.accs.e <- "accs/"
  axis.breaks = seq(0, 30, by = 5)
  idx.step = 10000
  
  c.red = '#CE1F6A'
  c.black = '#27374D'
  c.grey = '#EBEBEB'
  
  # Read the correspondence for one accession. The plot subsamples every
  # idx.step-th position, so only those positions are looked up in the interval
  # table -- the alignment is never expanded.
  ivl.acc = h5IvlRead(file.msa, acc)
  pan.len = h5AccLen(file.msa, acc)
  idx = seq(1, pan.len, idx.step)
  v.acc = data.frame(pan = idx, acc = ivlAccAt(ivl.acc, idx))
  v.acc = v.acc[v.acc$acc != 0,]
  
  pokaz(sum(v.acc < 0))
  
  # Generate a ggplot for the accession
  p = ggplot(v.acc, aes(x = pan, y = abs(acc), color=as.factor(sign(acc)))) + 
    # geom_abline(slope = 1, intercept = 0, color = c.grey) +
    geom_point(size = size) + 
    theme_bw() +
    labs(x = "Pangenome, Mbp", 
         y = paste0(acc, ', Mbp'), 
         title = NULL) +
    scale_color_manual(values = c("-1" = c.red, "1" = c.black)) +
    theme(legend.position = "none") +
    scale_x_continuous(limits = c(0, max(v.acc$pan)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, max(abs(v.acc$acc))), expand = c(0, 0)) +
    coord_fixed(ratio=1) 
  
  return(p)
  
}
