
# Dotplots

seq2mx <- function(seq, wsize){
  
  num_rows <- length(seq) - wsize + 1
  matrix_seq <- matrix(nrow = num_rows, ncol = wsize)
  for (i in 1:num_rows) {
    matrix_seq[i, ] <- seq[i:(i + wsize - 1)]
  }
  
  return(matrix_seq)
}

mxComp <- function(mx1, mx2, wsize, nmatch){
  mx.res = 0
  for(s in c('A', 'C', 'G', 'T')){
    mx.res = mx.res + (mx1 == s) %*% t(mx2 == s)
  }
  # mx.res = (mx.res >= nmatch) * 1
  mx.res[mx.res < nmatch] = 0
  
  indices <- which(mx.res != 0, arr.ind = TRUE)
  values <- mx.res[indices]
  result <- cbind(indices, values)
  result = as.data.frame(result)
  return(result)
}

dotplot <- function(seq1, seq2, wsize, nmatch) {
  seq2.rc = rev(seqinr::comp(seq2))
  
  mx1 = toupper(seq2mx(seq1, wsize))
  mx2 = toupper(seq2mx(seq2, wsize))
  
  result = mxComp(mx1, mx2, wsize, nmatch)
  
  mx2.rc = toupper(seq2mx(seq2.rc, wsize))
  result.rc = mxComp(mx1, mx2.rc, wsize, nmatch)
  result.rc$values = -result.rc$values
  result.rc$col = length(seq2) - result.rc$col - wsize + 2
  result = rbind(result.rc, result)
  
  
  p = ggplot(result, aes(x = row, y = col, fill = values, color=values)) +
    geom_tile(width = 1, height = 1, linewidth = 0.5) +
    # xlab(name1) + ylab(name2) +
    # xlab(paste0(strsplit(name1, '\\|')[[1]][7:9], collapse = '|')) + 
    # ylab(paste0(strsplit(name2, '\\|')[[1]][7:9], collapse = '|')) +
    xlab(NULL) + ylab(NULL) +
    theme_minimal() + coord_fixed() +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradient2(low = "#CE1F6A", mid = "white", high = "#27374D", midpoint = 0) +
    scale_color_gradient2(low = "#CE1F6A", mid = "white", high = "#27374D", midpoint = 0) +
    theme(panel.border = element_rect(colour = "grey", fill=NA, size=1)) +
    guides(fill = FALSE, color = NULL) 
  p 
  return(p + theme(legend.position = "none"))
}

# Sequence sould be just a vector of letters
# > head(seq1, 10)
# [1] "a" "a" "t" "t" "a" "t" "a" "t" "c" "t"
seqComplexity <- function(seq1, method='dotplot', wsize=10, nmatch=9) {
  
  mx1 = toupper(seq2mx(seq1, wsize))
  result = mxComp(mx1, mx1, wsize, nmatch)
  
  seq1.rc = rev(seqinr::comp(seq1))
  mx1.rc = toupper(seq2mx(seq1.rc, wsize))
  result.rc = mxComp(mx1, mx1.rc, wsize, nmatch)
  
  n.match = (nrow(result) + nrow(result.rc)) / length(seq1) 
  
  return(n.match)
}
