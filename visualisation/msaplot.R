msaplot <- function(seqs.mx, 
                    msa.cols = c("A" = "#8ACD9D", "C" = "#EE7571", "G" = "#7E9CC8", "T" = 
"#FFD97C", '-'='#EEEDEF')){
  
  if (is.vector(seqs.mx)){
    seqs.mx = t(matrix(seqs.mx))
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


