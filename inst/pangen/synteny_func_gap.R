transform_positions <- function(df) {
  df.gap <- df[, c('V2', 'V3', 'V4', 'V5', 'idx', 'dir')]
  tmp <- df.gap$V4[df.gap$dir == 1]
  df.gap$V4[df.gap$dir == 1] <- df.gap$V5[df.gap$dir == 1]
  df.gap$V5[df.gap$dir == 1] <- tmp
  df.gap$dir <- 0
  
  x.min <- min(df.gap$V2)
  y.min <- min(df.gap$V4)
  df.gap[, c('V2', 'V3')] <- df.gap[, c('V2', 'V3')] - x.min + 1
  df.gap[, c('V4', 'V5')] <- df.gap[, c('V4', 'V5')] - y.min + 1
  df.gap$len.y <- df.gap$V5 - df.gap$V4 + 1
  
  return(df.gap)
}

greedy_loop <- function(df.gap, overlap.cutoff = 0.2) {
  y.ticks <- sort(unique(c(df.gap$V4, df.gap$V5 + 1)))
  cell.y.pos <- data.frame(beg = y.ticks[-length(y.ticks)], end = y.ticks[-1] - 1)
  cell.y.pos$len <- cell.y.pos$end - cell.y.pos$beg + 1
  pos.y.attr <- rep(0, max(cell.y.pos$end))
  
  df.gap$X4 <- 0
  df.gap$X5 <- 0
  y.list <- list()
  
  for (irow in 1:nrow(df.gap)) {
    pos.tmp <- which((cell.y.pos$beg >= df.gap$V4[irow]) & (cell.y.pos$end <= df.gap$V5[irow]))
    df.gap$X4[irow] <- min(pos.tmp)
    df.gap$X5[irow] <- max(pos.tmp)
  }
  
  df.gap$len.x <- df.gap$V3 - df.gap$V2 + 1
  df.gap$len.y <- abs(df.gap$V4 - df.gap$V5) + 1
  
  pos.y.occ <- rep(0, nrow(cell.y.pos))
  pos.x <- 0
  idx.added <- c()
  
  df.gap <- df.gap[order(df.gap$V2), ]
  
  while (TRUE) {
    d.x <- df.gap$V3 - pos.x
    idx.next <- which((d.x > 0) & (d.x / df.gap$len.x > (1 - overlap.cutoff)))
    
    if (length(idx.next) == 0) break
    
    i.next.found <- 0
    for (i.next in idx.next) {
      pos.y.i.next <- df.gap$X4[i.next]:df.gap$X5[i.next]
      pos.y.i.next.overlap <- pos.y.i.next[pos.y.occ[pos.y.i.next] > 0]
      
      y.overlap <- sum(cell.y.pos$len[pos.y.i.next.overlap]) / df.gap$len.y[i.next]
      if (y.overlap >= overlap.cutoff) next
      
      i.next.found <- i.next
      break
    }
    
    if (i.next.found == 0) break
    
    idx.added <- c(idx.added, i.next.found)
    pos.y.occ[pos.y.i.next] <- pos.y.occ[pos.y.i.next] + 1
    pos.x <- df.gap$V3[i.next]
  }
  
  return(idx.added)
}


posShift <- function(x.gap){
  
  # Set correct position
  x.gap$beg.q = as.numeric(sapply(strsplit(x.gap[,1], "\\|"), "[", 2)) - 1
  x.gap$end.q = as.numeric(sapply(strsplit(x.gap[,1], "\\|"), "[", 3)) - 1
  x.gap$len.q = x.gap$end.q - x.gap$beg.q + 1
  x.gap$beg.r = as.numeric(sapply(strsplit(x.gap[,10], "\\|"), "[", 2)) - 1
  x.gap$end.r = as.numeric(sapply(strsplit(x.gap[,10], "\\|"), "[", 3)) - 1
  x.gap$len.r = x.gap$end.r - x.gap$beg.r + 1
  
  # Modify Q-positions
  pos.q = unique(x.gap[, c('V1','beg.q', 'end.q', 'len.q')])
  pos.q = pos.q[order(pos.q$beg.q),]
  if(nrow(pos.q) > 1){
    pos.q$shift = c(0, cumsum(pos.q$len.q[-nrow(pos.q)]))  
  } else {
    pos.q$shift = 0
  }
  rownames(pos.q) = pos.q$V1
  
  # Modify R-positions
  pos.r = unique(x.gap[, c('V10','beg.r', 'end.r', 'len.r')])
  pos.r = pos.r[order(pos.r$beg.r),]
  if(nrow(pos.r) > 1){
    pos.r$shift = c(0, cumsum(pos.r$len.r[-nrow(pos.r)]))  
  } else {
    pos.r$shift = 0
  }
  rownames(pos.r) = pos.r$V10
  
  return(list(q = pos.q, r = pos.r))
}
