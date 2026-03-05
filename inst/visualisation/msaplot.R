#' Plot Multiple Sequence Alignment (MSA)
#'
#' This function creates a visual representation of a multiple sequence alignment (MSA) using ggplot2.
#' It supports both nucleotide (`nt`) and amino acid (`aa`) sequences, and colors the alignment according 
#' to the specified color scheme.
#'
#' @param aln A matrix or vector representing the sequences to be plotted. Each row corresponds to 
#' a sequence, and each column corresponds to a position in the alignment.
#' @param msa.cols A named vector of colors to be used for each nucleotide or amino acid. If `NULL`, 
#' a default color scheme will be used based on the sequence type (`nt` or `aa`).
#' @param seq.type A character string indicating the type of sequences provided: either `'nt'` 
#' for nucleotide sequences or `'aa'` for amino acid sequences. Default is `'nt'`.
#'
#' @return A ggplot object representing the MSA plot.
#' 
#' @details The function first checks the `seq.type` to determine whether the sequences are nucleotides 
#' or amino acids, and sets a default color scheme if none is provided. It then processes the sequences, 
#' converting them into a format suitable for plotting with ggplot2. The plot is generated using `geom_tile` 
#' to represent each position in the alignment, with colors corresponding to the nucleotide or amino acid 
#' at each position.
#'
#' @export
msaplot <- function(aln, seq.type='nt', msa.cols = NULL, show.legend=F){
  
  # Input handling
  if(is.vector(aln) && is.character(aln)){
    aln <- aln2mx(aln)
  }
  
  # If input is matrix, ensure correct type
  if(!is.matrix(aln)){
    stop("Input must be either a character vector of aligned sequences or a matrix.")
  }
  
  # Colors
  if(is.null(msa.cols)){
    if(seq.type == 'nt'){
      msa.cols = c("A" = "#8ACD9D", 
                   "C" = "#EE7571", 
                   "G" = "#7E9CC8", 
                   "T" = "#FFD97C", 
                   '-' = '#EEEDEF',
                   'N' = '#31363F')
    } else if(seq.type == 'aa'){
      msa.cols <- c(
        "A" = "#8EDB8E",  # Light green
        "G" = "#8EDB8E",  # Light green
        "C" = "#ACEC78",  # Green
        "D" = "#70B94A",  # Dark green
        "E" = "#70B94A",  # Dark green
        "N" = "#70B94A",  # Dark green
        "Q" = "#70B94A",  # Dark green
        "I" = "#7BB9F9",  # Blue
        "L" = "#7BB9F9",  # Blue
        "M" = "#7BB9F9",  # Blue
        "V" = "#7BB9F9",  # Blue
        "F" = "#9999F8",  # Lilac
        "W" = "#9999F8",  # Lilac
        "Y" = "#9999F8",  # Lilac
        "H" = "#5555F6",  # Dark blue
        "K" = "#F7CE83",  # Orange
        "R" = "#F7CE83",  # Orange
        "P" = "#E4ADAA",  # Pink
        "S" = "#EC545B",  # Red
        "T" = "#EC545B",  # Red
        "-" = "#EEEDEF",  # Gap
        "X" = "#31363F"   # Unknown
      )
    } else {
      stop('Wrong type of sequences. Only `nt` or `aa` are allowed')
    }
  }
  
  # Row names
  if(is.null(row.names(aln))){
    pokazAttention('Names of sequences are not provided. They will be .')
    row.names(aln) = paste0('s.', 1:nrow(aln))
  }
  
  if(length(unique(row.names(aln))) != nrow(aln)){
    pokazAttention('Names of sequences are not unique. They were modified.')
    row.names(aln) = paste0('s.', 1:nrow(aln))
  }
  
  if (is.vector(aln)){
    aln = t(matrix(aln))
  }
  
  if(is.null(rownames(aln))){
    rownames(aln) = paste('s', 1:nrow(aln), sep = '_')
  }
  
  aln = toupper(aln)
  df <- reshape2::melt(aln)
  df$Var1 = factor(df$Var1, levels = rev(rownames(aln)))
  df$Var2 = as.numeric(df$Var2)
  
  names(msa.cols) = toupper(names(msa.cols))
  
  g.msa = ggplot(df, aes(x = Var2, y = Var1, fill = value, color = value)) + 
    geom_tile() +
    scale_fill_manual(values = msa.cols) +
    scale_color_manual(values = msa.cols) +
    theme_bw() + 
    scale_x_continuous(limits = c(0, ncol(aln)+1), expand = c(0, 0)) +
    theme(panel.grid = element_blank(),
          panel.border = element_blank()) + ylab('') + 
    xlab(NULL) +
    theme(legend.position = "none")
  
  if(show.legend){
    g.msa = g.msa + theme(legend.position = "bottom") + 
      guides(color = guide_legend(title = NULL),
             fill = guide_legend(title = NULL, 
                                 override.aes = list(colour = "black"))) 
  }
  
  return(g.msa)
  
}

#' Highlight Differences in Multiple Sequence Alignment (MSA)
#'
#' This function identifies and highlights differences in a multiple sequence alignment (MSA) relative 
#' to a reference sequence. It marks positions that are the same as, different from, or gaps compared 
#' to the reference sequence.
#'
#' @param aln A matrix representing the sequences to be analyzed. Each row corresponds to a sequence, 
#' and each column corresponds to a position in the alignment.
#' @param i.ref An integer indicating the index of the reference sequence within `aln`. The default is 1.
#'
#' @return A ggplot object representing the MSA plot, with differences, similarities, and gaps highlighted 
#' in different colors.
#' 
#' @details The function compares each sequence in the alignment to a reference sequence. Positions that 
#' differ from the reference are marked as "diff", positions that are the same are marked as "same", and 
#' gaps are marked as "gap". These differences are then visualized using the `msaplot` function.
#'
#' @export
msadiff <- function(aln, i.ref=1, show.legend=F){
  
  # Input handling
  if(is.vector(aln) && is.character(aln)){
    aln <- aln2mx(aln)
  }
  
  # If input is matrix, ensure correct type
  if(!is.matrix(aln)){
    stop("Input must be either a character vector of aligned sequences or a matrix.")
  }
  
  aln = toupper(aln)
  bin.mx <- t(apply(aln, 1, function(row) as.integer(row != aln[i.ref, ]))) + 2
  bin.mx[,aln[i.ref,] == '-'] = 2
  bin.mx[aln == '-'] = 1
  values = c('gap', 'same', 'diff')
  bin.mx <- t(apply(bin.mx, 1, function(row) values[row]))
  
  b.msa = msaplot(bin.mx,
                  msa.cols = c("same" = "grey80", "diff" = "grey20", "gap" = "white"),
                  show.legend = show.legend)
  
  return(b.msa)
  
}


