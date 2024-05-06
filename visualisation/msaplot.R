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


msadiff <- function(seqs.mx, i.ref=1){
  
  bin.mx <- t(apply(seqs.mx, 1, function(row) as.integer(row != seqs.mx[i.ref, ]))) + 2
  bin.mx[,seqs.mx[i.ref,] == '-'] = 2
  bin.mx[seqs.mx == '-'] = 1
  values = c('gap', 'same', 'diff')
  bin.mx <- t(apply(bin.mx, 1, function(row) values[row]))
  
  b.msa = msaplot(bin.mx,
                  msa.cols = c("same" = "grey80", "diff" = "grey20", "gap" = "white"))
  
  return(b.msa)
  
}


