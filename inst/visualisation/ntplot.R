#' Plot nucleotide frequencies along a sequence
#'
#' This function calculates and visualizes the frequencies of nucleotides (A, C, G, T) along a given sequence
#' using a sliding window approach.
#'
#' @param sequence A character vector representing the nucleotide sequence.
#' @param wnd An integer specifying the sliding window size (default: 40).
#' @param nt.separate A logical indicating whether to plot each nucleotide separately using facets (default: FALSE).
#'
#' @return A ggplot object showing the nucleotide frequency distribution along the sequence.
#'
#' @examples
#' seq_example <- sample(c("A", "C", "G", "T"), 1000, replace = TRUE)
#' ntplot(seq_example, wnd = 50, nt.separate = TRUE)
#'
#' @import ggplot2
#' @import reshape2
#' @export
ntplot <- function(sequence, wnd=100, nt.separate = F) {
  
  msa.cols = c("A" = "#8ACD9D", 
               "C" = "#EE7571", 
               "G" = "#7E9CC8", 
               "T" = "#FFD97C", 
               '-' = '#EEEDEF',
               'N' = '#31363F')
  
  sequence = toupper(sequence)
  
  positions <- seq(1, length(sequence) - wnd + 1)
  result <- data.frame(Position = positions)
  
  for (s.nt in c("A", "C", "G", "T")) {
    result[[s.nt]] <- sapply(positions, function(pos) {
      sum(sequence[pos:(pos + wnd - 1)] == s.nt) / wnd
    })
  }
  
  df.freq <- reshape2::melt(result, id.vars = "Position", variable.name = "Base", value.name = "Frequency")
  df.freq$Base = factor(df.freq$Base, levels = c('A', 'T', 'G', 'C'))
  
  p = ggplot(df.freq, aes(x = Position, y = Frequency, fill = Base)) +
    geom_bar(stat = "identity", width = 1) +
    theme_minimal()+ scale_fill_manual(values = msa.cols) +
    theme(legend.position = "bottom", legend.direction = "horizontal")
  
  if(nt.separate){
    p <- p + facet_wrap(~Base, scales = "free_y", ncol = 1) + theme(legend.position = "none") 
  }
  p
  
  return(p)
}

#' @export
ntplot.s <- function(s,  wnd = 100, nt.separate = F, ...) {
  return(ntplot(seq2nt(s),  wnd=wnd, nt.separate = nt.separate, ...))
}

