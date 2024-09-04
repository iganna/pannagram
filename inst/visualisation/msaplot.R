#' Plot Multiple Sequence Alignment (MSA)
#'
#' This function creates a visual representation of a multiple sequence alignment (MSA) using ggplot2.
#' It supports both nucleotide (`nt`) and amino acid (`aa`) sequences, and colors the alignment according 
#' to the specified color scheme.
#'
#' @param seqs.mx A matrix or vector representing the sequences to be plotted. Each row corresponds to 
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
#' @examples
#' # Example nucleotide sequences
#' nt_sequences <- matrix(c("A", "T", "G", "C", "A", "T", "-", "G", "C", "A"), nrow = 2, byrow = TRUE)
#' msaplot(nt_sequences, seq.type = 'nt')
#'
#' # Example amino acid sequences
#' aa_sequences <- matrix(c("A", "G", "C", "D", "E", "F", "G", "H", "I", "K"), nrow = 2, byrow = TRUE)
#' msaplot(aa_sequences, seq.type = 'aa')
#'
#' @export
msaplot <- function(seqs.mx, msa.cols = NULL, seq.type='nt'){
  
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
  

  if (is.vector(seqs.mx)){
    seqs.mx = t(matrix(seqs.mx))
  }
  
  if(is.null(rownames(seqs.mx))){
    rownames(seqs.mx) = paste('s', 1:nrow(seqs.mx), sep = '_')
  }
  
  seqs.mx = toupper(seqs.mx)
  df <- reshape2::melt(seqs.mx)
  df$Var1 = factor(df$Var1, levels = rev(rownames(seqs.mx)))
  df$Var2 = as.numeric(df$Var2)
  
  names(msa.cols) = toupper(names(msa.cols))
  
  g.msa = ggplot(df, aes(x = Var2, y = Var1, fill = value, color = value)) + 
    geom_tile() +
    scale_fill_manual(values = msa.cols) +
    scale_color_manual(values = msa.cols) +
    theme_bw() + 
    scale_x_continuous(limits = c(0, ncol(seqs.mx)+1), expand = c(0, 0)) +
    theme(panel.grid = element_blank(),
          panel.border = element_blank()) + ylab('') + 
    xlab(NULL) +
    theme(legend.position = "none")
  
  return(g.msa)
  
}

#' Highlight Differences in Multiple Sequence Alignment (MSA)
#'
#' This function identifies and highlights differences in a multiple sequence alignment (MSA) relative 
#' to a reference sequence. It marks positions that are the same as, different from, or gaps compared 
#' to the reference sequence.
#'
#' @param seqs.mx A matrix representing the sequences to be analyzed. Each row corresponds to a sequence, 
#' and each column corresponds to a position in the alignment.
#' @param i.ref An integer indicating the index of the reference sequence within `seqs.mx`. The default is 1.
#'
#' @return A ggplot object representing the MSA plot, with differences, similarities, and gaps highlighted 
#' in different colors.
#' 
#' @details The function compares each sequence in the alignment to a reference sequence. Positions that 
#' differ from the reference are marked as "diff", positions that are the same are marked as "same", and 
#' gaps are marked as "gap". These differences are then visualized using the `msaplot` function.
#'
#' @examples
#' # Example sequences
#' sequences <- matrix(c("A", "T", "G", "C", "A", "T", "-", "G", "C", "A"), nrow = 2, byrow = TRUE)
#' msadiff(sequences, i.ref = 1)
#'
#' @export
msadiff <- function(seqs.mx, i.ref=1){
  
  seqs.mx = toupper(seqs.mx)
  bin.mx <- t(apply(seqs.mx, 1, function(row) as.integer(row != seqs.mx[i.ref, ]))) + 2
  bin.mx[,seqs.mx[i.ref,] == '-'] = 2
  bin.mx[seqs.mx == '-'] = 1
  values = c('gap', 'same', 'diff')
  bin.mx <- t(apply(bin.mx, 1, function(row) values[row]))
  
  b.msa = msaplot(bin.mx,
                  msa.cols = c("same" = "grey80", "diff" = "grey20", "gap" = "white"))
  
  return(b.msa)
  
}


