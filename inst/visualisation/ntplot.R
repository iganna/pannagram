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
  
  sequence <- prepareNtSeq(sequence)
  
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
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "bottom", legend.direction = "horizontal")
  
  if(nt.separate){
    p <- p + facet_wrap(~Base, scales = "free_y", ncol = 1) + theme(legend.position = "none") 
  }
  p
  
  return(p)
}

#' Calculate GC Content in DNA Sequences
#'
#' This function calculates the GC content for a given set of DNA sequences, using a specified window size.
#'
#' @param sequences A character vector of DNA sequences. Each element of the vector should be a DNA sequence in the form of a string.
#' @param wnd An integer specifying the window size for calculating the GC content. Default is 100.
#' @return A list where each element corresponds to a sequence from the input `s`. Each element is a numeric vector
#' representing the GC content for each window in the sequence.
#' @examples
#' # Example usage
#' sequences <- c("ATGCGATCGATCG", "GCGCGCGCGCGCG", "ATATATATATAT")
#' gcContent(sequences, wnd = 5)
#'
#' @export
gcContent <- function(sequences, wnd = NULL){

  if(is.null(wnd)){
    gc.list = c()
  } else {
    gc.list = list()
  }
  
  sequences = toupper(sequences)
  
  for(i.s in seq_along(sequences)){
    s = sequences[i.s]
    
    if(is.null(wnd)){
      s = seq2nt(s)
      gc.tmp = sum((s == 'G') | (s == 'C')) / length(s)
      gc.list[i.s] = gc.tmp
    } else if(nchar(s[i.seq]) <= wnd){
      s = seq2nt(s)
      gc.tmp = sum((s == 'G') | (s == 'C')) / length(s)
      gc.list[[i.s]] = gc.tmp
    } else {
      m = splitSeq(s[i.seq], n = wnd, merge = F)
      gc.tmp = rowSums(m == 'G') + rowSums(m == 'C')
      gc.tmp = gc.tmp / wnd  
      gc.list[[i.s]] = gc.tmp
    }
  }
  names(gc.list) = names(sequences)
  return(gc.list)
}



