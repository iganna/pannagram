#' Compute Distance Matrices for Blast hits
#'
#'
computeDistBwBlastHits <- function(b, p = 0.5) {
  # Order by V2
  b$idx.order = 1:nrow(b)
  b <- b[order(b$V2),]
  n <- nrow(b)
  dx <- matrix(0, nrow = n, ncol = n)
  dy <- matrix(0, nrow = n, ncol = n)
  
  # Compute rectangle centers and dimensions
  b$x <- (b$V2 + b$V3) / 2
  b$y <- (b$V4 + b$V5) / 2
  b$xlen <- abs(b$V3 - b$V2)
  b$ylen <- abs(b$V5 - b$V4)
  b$xd <- b$xlen / 2
  b$yd <- b$ylen / 2
  
  # Compute pairwise combinations
  combinations <- t(combn(1:n, 2))
  for (k in 1:nrow(combinations)) {
    j <- combinations[k, 1]
    i <- combinations[k, 2]
    
    dx.tmp <- abs(b$x[i] - b$x[j]) - b$xd[i] - b$xd[j]
    dy.tmp <- abs(b$y[i] - b$y[j]) - b$yd[i] - b$yd[j]
    
    # Apply overlap correction
    if (dx.tmp < 0) {
      if (abs(dx.tmp) > p * b$xlen[i]) dx.tmp <- Inf
      if (abs(dx.tmp) > p * b$xlen[j]) dx.tmp <- Inf
    }
    
    if (dy.tmp < 0) {
      if (abs(dy.tmp) > p * b$ylen[i]) dy.tmp <- Inf
      if (abs(dy.tmp) > p * b$ylen[j]) dy.tmp <- Inf
    }
    
    dx[i, j] <- dx[j, i] <- dx.tmp
    dy[i, j] <- dy[j, i] <- dy.tmp
    
    # Increasing only
    if(b$y[i] > b$y[j]){
      dy[i, j] <- Inf
    }
    
    if(b$y[j] > b$y[i]){
      dy[j, i] <- Inf
    }
    
    if(b$x[i] > b$x[j]){
      dx[i, j] <- Inf
    }
    
    if(b$x[j] > b$x[i]){
      dx[j, i] <- Inf
    }
    
  }
  
  # Final distance matrix
  d <- dx + dy
  d[d < 0] <- 0
  d[d == 0] <- 1
  
  # Return distance matrix
  return(list(dist = d, order = b$idx.order))
}
