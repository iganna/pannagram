#' Plot ORFs with Horizontal Arrows
#'
#' This function takes a data frame containing information about Open Reading Frames (ORFs) 
#' and plots them as horizontal arrows. Each arrow represents an ORF, 
#' with the direction indicated by the `strand` column. 
#' The function automatically adjusts the vertical positioning of the arrows to prevent overlap.
#'
#' @param df A data frame with columns `beg`, `end`, and `strand`, 
#'          where `beg` and `end` are the start and end positions of the ORFs, 
#'          and `strand` indicates the direction (`+` for forward, `-` for reverse).
#' @return A ggplot object representing the ORFs as horizontal arrows. 
#' 
#' @examples
#' df <- data.frame(
#'   beg = c(100, 500, 900),
#'   end = c(400, 800, 1200),
#'   strand = c('+', '-', '+')
#' )
#' orfplot(df)
#' @export
orfplot <- function(df, optimal = F, 
                    s.color = 'strand', 
                    y = NULL,
                    show.legend = F,
                    arrow.size = 0.05){
  
  if(optimal){
    # Initialize row.number for vertical positioning and a position array to track filled positions
    df$row.number = 0
    pos.array = matrix(0, 1, max(c(df$beg, df$end)))
    seq.dist = 50  # Buffer distance between ORFs
    
    # Loop through each ORF to assign a row without overlap
    for(ipos in 1:nrow(df)){
      flag.row = F  # Flag to track if a row has been found
      
      # Determine start and end positions with buffer
      p1 = max(1, (min(df$beg[ipos], df$end[ipos]) - seq.dist))
      p2 = min(ncol(pos.array), (max(df$beg[ipos], df$end[ipos]) + seq.dist))
      
      idx.pos = p1:p2
      # Check each row for available space
      for(irow in 1:nrow(pos.array)){
        if(sum(pos.array[irow, idx.pos]) == 0){
          df$row.number[ipos] = irow
          pos.array[irow, idx.pos] = 1
          flag.row = T
          break
        } 
      }
      
      # If no row found, add a new row to pos.array
      if(flag.row == F){
        pos.array = rbind(pos.array, 0)
        irow = nrow(pos.array)
        df$row.number[ipos] = irow
        pos.array[irow, idx.pos] = 1
      }
    }
  } else {
    df$row.number = 1:nrow(df)
  }
  
  if(!is.null(y)){
    df$row.number = y
  }
  
  # Plot ORFs 
  if(s.color %in% colnames(df)){
    p.orf <- ggplot(df, aes(x = beg, xend = end, y = row.number, yend = row.number, colour = !!sym(s.color))) +
      geom_segment(arrow = arrow(length = unit(arrow.size, "inches")), size = 1) + 
      theme_minimal() + xlab(NULL) + ylab(NULL) 
    
    if(!show.legend){
      p.orf = p.orf + theme(legend.position = "none", 
                    axis.text.y = element_blank(),  
                    axis.ticks.y = element_blank()) 
    }
    
    if(sum(unique(df[,s.color]) %in% c('+', '-')) == 2){
      p.orf = p.orf + scale_colour_manual(values = c('-' = '#40679E', '+' = '#FF407D')) 
    }
    
  } else {
    p.orf <- ggplot(df, aes(x = beg, xend = end, y = row.number, yend = row.number, color = as.factor(row.number))) +
      geom_segment(arrow = arrow(length = unit(arrow.size, "inches")), size = 1) + 
      theme_minimal() + xlab(NULL) + ylab(NULL) +  
      theme(legend.position = "none", 
            axis.text.y = element_blank(),  
            axis.ticks.y = element_blank())
  }
  
  
  return(p.orf) 
}

#' Find the Longest ORFs in DNA Sequences
#'
#' @param seqs A list or vector of DNA sequences.
#' @param n.best An integer specifying the number of longest ORFs to select for each sequence (default is 1).
#' @return A named list of the longest ORFs found, where names indicate the original sequence and ORF identifier.
#' 
#' @export
orfBest <- function(seqs, n.best = 1) {
  seqs.names <- names(seqs)
  if (is.null(seqs.names)) {
    seqs.names <- paste0('s_', seq_along(seqs))
  }
  
  orf.best <- lapply(seq_along(seqs), function(i) {
    res <- orfFinder(seqs[[i]])
    if (is.null(res$pos)) return(NULL)
    
    orf.tmp <- head(res$orf, n.best)
    names(orf.tmp) <- paste(seqs.names[i], names(orf.tmp), sep = '|')
    return(orf.tmp)
  })
  
  orf.best = unlist(orf.best, recursive = FALSE)
  
  if(is.null(orf.best)) return(NULL)
  
  orf.best = orf.best[order(-nchar(orf.best))]
  
  return(orf.best)
}
