mergeOverlaps <- function(breaks){
  
  # Input validation: check if the input is a data frame
  if (!is.data.frame(breaks)) {
    stop("Input must be a data frame.")
  }
  
  # Input validation: check if the data frame is empty
  if (nrow(breaks) == 0) {
    warning("Input data frame is empty.")
    return(breaks)
  }
  
  # Input validation: check if the required columns are present
  required_cols <- c('idx.beg', 'idx.end')
  if (!all(required_cols %in% colnames(breaks))) {
    stop("The input data frame must contain the following columns: 'idx.beg', 'idx.end'.")
  }
  
  # Input validation: check if 'idx.beg' and 'idx.end' are numeric
  if (!is.numeric(breaks$idx.beg) || !is.numeric(breaks$idx.end)) {
    stop("'idx.beg' and 'idx.end' must be numeric.")
  }
  
  # ---- Merge coverages ----
  n.init = nrow(breaks)
  breaks = breaks[,c('idx.beg', 'idx.end')]
  
  breaks = breaks[order(-breaks$idx.end),]
  breaks = breaks[order(breaks$idx.beg),]
  
  breaks$id = 1:nrow(breaks)
  
  breaks = breaks[!duplicated(breaks[, c('idx.beg', 'idx.end')]),]
  breaks$cnt = c(breaks$id[-1], n.init+1) - breaks$id
  if(!(sum(breaks$cnt) == n.init)) stop('Checkpoint1')
  if (any(breaks$cnt < 0)) stop('Checkpoint2')
  
  breaks <- breaks[order(-breaks$idx.end), ]
  breaks <- breaks[order(breaks$idx.beg), ]
  
  
  # Merge overlaps
  n = 0
  while (n != nrow(breaks)) {
    n = nrow(breaks)
    
    idx_full_cover = which(breaks$idx.beg[-1] <= breaks$idx.end[-nrow(breaks)])
    idx_full_cover = setdiff(idx_full_cover, idx_full_cover + 1)
    
    if (length(idx_full_cover) == 0) break
    
    breaks$cnt[idx_full_cover] = breaks$cnt[idx_full_cover] + breaks$cnt[idx_full_cover + 1]
    breaks$idx.end[idx_full_cover] = pmax(breaks$idx.end[idx_full_cover], breaks$idx.end[idx_full_cover + 1])
    
    breaks = breaks[-(idx_full_cover + 1), ]
  }
  
  if(is.unsorted(breaks$idx.beg)) stop('Chrckpoint sorted 1')
  if(is.unsorted(breaks$idx.end)) stop('Chrckpoint sorted 2')
  
  breaks$len = breaks$idx.end - breaks$idx.beg + 1
  
  return(breaks)
}



# mergeOverlaps2 <- function(breaks){
#   
#   # Input validation: check if the input is a data frame
#   if (!is.data.frame(breaks)) {
#     stop("Input must be a data frame.")
#   }
#   
#   # Input validation: check if the data frame is empty
#   if (nrow(breaks) == 0) {
#     warning("Input data frame is empty.")
#     return(breaks)
#   }
#   
#   # Input validation: check if the required columns are present
#   required_cols <- c('idx.beg', 'idx.end')
#   if (!all(required_cols %in% colnames(breaks))) {
#     stop("The input data frame must contain the following columns: 'idx.beg', 'idx.end'.")
#   }
#   
#   # Input validation: check if 'idx.beg' and 'idx.end' are numeric
#   if (!is.numeric(breaks$idx.beg) || !is.numeric(breaks$idx.end)) {
#     stop("'idx.beg' and 'idx.end' must be numeric.")
#   }
#   
#   # ---- Merge coverages ----
#   n.init = nrow(breaks)
#   
#   breaks <- breaks[order(-breaks$idx.end), ]
#   breaks <- breaks[order(breaks$idx.beg), ]
#   breaks$id = 1:nrow(breaks)
#   breaks$cnt = 1
#   
#   breaks$len.add = breaks$len.acc - breaks$len.comb
#   
#   breaks = breaks[,c('idx.beg', 'idx.end', 'len.add')]
#   breaks = breaks[!duplicated(breaks[, c('idx.beg', 'idx.end')]),]
#   
#   breaks <- aggregate(len.add ~ idx.beg + idx.end, data = breaks, FUN = mean)
#   
#   breaks.len = breaks$len.add + breaks$idx.end - breaks$idx.beg + 1
#   
#   # Merge overlaps
#   n = 0
#   while (n != nrow(breaks)) {
#     n = nrow(breaks)
#     
#     idx_full_cover = which(breaks$idx.beg[-1] <= breaks$idx.end[-nrow(breaks)])
#     # 
#     new.lengths = pmax(breaks$idx.end[idx_full_cover], breaks$idx.end[idx_full_cover + 1]) - 
#       breaks$idx.beg[idx_full_cover] + (breaks$len.plus[idx_full_cover] + breaks$
#     
#     if(any(is.na(new.lengths))) stop("NA are found")
#     
#     idx.large = new.lengths > len.large
#     if(sum(idx.large) > 0){
#       idx.rm = idx_full_cover[idx.large]  # Deside which to remove
#       
#       # Do not consider both in this iteration
#       idx_full_cover = setdiff(idx_full_cover, idx.rm)
#       idx_full_cover = setdiff(idx_full_cover, idx.rm + 1)
#       
#       idx.rm = ifelse(breaks$cnt[idx.rm] < breaks$cnt[idx.rm+1], idx.rm, idx.rm+1)
#     } else {
#       idx.rm = c()
#     }
#     
#     idx_full_cover = setdiff(idx_full_cover, idx_full_cover + 1)
#     
#     if (length(idx_full_cover) == 0) break
#     
#     breaks$cnt[idx_full_cover] = breaks$cnt[idx_full_cover] + breaks$cnt[idx_full_cover + 1]
#     breaks$idx.end[idx_full_cover] = pmax(breaks$idx.end[idx_full_cover], breaks$idx.end[idx_full_cover + 1])
#     
#     b.len[breaks$id[idx_full_cover]] = b.len[breaks$id[idx_full_cover]] + b.len[breaks$id[idx_full_cover+1]]
#     
#     if(length(idx.rm) == 0){
#       b.len[breaks$id[idx_full_cover+1],] = NA
#       breaks = breaks[-(idx_full_cover + 1), ]
#     } else {
#       stop('test')
#       idx.rm = c(idx.rm, idx_full_cover+1)
#       if(length(idx.rm) != length(unique(idx.rm))) stop('Problem with idx.rm')
#       b.len[breaks$id[idx.rm],] = NA
#       breaks = breaks[-idx.rm, ]
#     }
#   }
#   
#   
#   
#   
#   if(is.unsorted(breaks$idx.beg)) stop('Chrckpoint sorted 1')
#   if(is.unsorted(breaks$idx.end)) stop('Chrckpoint sorted 2')
#   
#   breaks$len = breaks$idx.end - breaks$idx.beg + 1
#   
#   return(breaks)
# }



solveLong <- function(breaks, breaks.init, len.large) {
  # Identify breaks where the length exceeds the len.large
  idx.to.solve = which(breaks$len > len.large)
  
  breaks.init$id <- 1:nrow(breaks.init)
  
  idx.rem.init <- c()
  for (irow in idx.to.solve) {
    pos.b <- breaks$idx.beg[irow]
    pos.e <- breaks$idx.end[irow]
    n <- pos.e - pos.b + 1
    
    # Select initial breaks that fit within the current merged break
    breaks.tmp <- breaks.init[(breaks.init$idx.beg >= pos.b) & (breaks.init$idx.end <= pos.e),]
    
    # Sort by length (starting from the largest)
    breaks.tmp <- breaks.tmp[order(-breaks.tmp$len.comb),]
    
    # Determine the number of breaks which should be removed
    for (n.rem in 1:nrow(breaks.tmp)) {
      breaks.tmp.m <- mergeOverlaps(breaks.tmp[-(1:n.rem),])  # Merge overlaps for remaining breaks
      if (max(breaks.tmp.m$len) < len.large) break
    }
    
    # pokaz('Potential removes:', n.rem)
    
    # Maybe put some breaks back
    idx.rem <- c()
    for (i in 1:n.rem) {
      # Check for overlaps between the remaining breaks
      tmp <- (breaks.tmp.m$idx.beg <= breaks.tmp$idx.end[i]) & 
        (breaks.tmp.m$idx.end >= breaks.tmp$idx.beg[i])
      
      if (sum(tmp) == 0) {  # If no overlap
        tmp <- breaks.tmp[i,]
        
        # Assign a dummy variables
        tmp$len <- tmp$len.comb
        tmp$id <- -1  
        tmp$cnt <- 1
        
        # Add the break into merged breaks
        breaks.tmp.m <- rbind(breaks.tmp.m, tmp[, colnames(breaks.tmp.m)]) 
        next
      }
      
      if (sum(tmp) != 1) {  # If more than one overlap
        idx.rem <- c(idx.rem, i)  # Remove it
        next
      }
      
      k <- which(tmp)  # Get the overlap index
      
      # Adjust boundaries 
      p.b <- min(breaks.tmp.m$idx.beg[k], breaks.tmp$idx.beg[i])
      p.e <- max(breaks.tmp.m$idx.end[k], breaks.tmp$idx.end[i])
      if ((p.e - p.b + 1) < len.large) {
        breaks.tmp.m$idx.beg[k] <- p.b
        breaks.tmp.m$idx.end[k] <- p.e
      } else {
        idx.rem <- c(idx.rem, i)  # Remove it if too large
      }
    }
    
    # pokaz('Actually removed:', length(idx.rem))
    
    # Save
    idx.rem.init <- c(idx.rem.init, breaks.tmp[idx.rem,]$id)
  }
  
  return(idx.rem.init)
}


#' Find Breaks in a Vector
#'
#' This function identifies breaks in the input vector `v`, detects blocks based 
#' on ranks, and creates a data frame of the break points. It filters blocks based 
#' on specific conditions and adds attributes like accuracy and lengths to the result.
#'
#' @param v A numeric vector to analyze.
#' 
findBreaks <- function(v) {
  v[is.na(v)] = 0
  # Remove zeros and retain indices
  v.idx <- 1:length(v)
  i.rm = v != 0
  v.idx <- v.idx[i.rm]
  v <- v[i.rm]
  
  # Rank the values and adjust for negative ranks
  v.r <- rank(abs(v))
  v.r[v < 0] <- v.r[v < 0] * (-1)
  
  # Find continuous blocks in ranked data
  v.b <- findRuns(v.r)
  
  # Assign start and end values and their indices
  v.b$v.beg <- v[v.b$beg]
  v.b$v.end <- v[v.b$end]
  v.b$i.beg <- v.idx[v.b$beg]
  v.b$i.end <- v.idx[v.b$end]
  
  # Initialize an array for block accumulation
  blocks.acc <- rep(0, max(v))
  for(irow in 1:nrow(v.b)) {
    blocks.acc[abs(v.b$v.beg[irow]):abs(v.b$v.end[irow])] <- irow
  }
  
  # Identify breaks where consecutive values differ by more than 1
  i.br.acc <- which(abs(diff(v)) != 1)
  df <- data.frame(
    val.beg = v[i.br.acc],
    val.end = v[i.br.acc + 1],
    idx.beg = v.idx[i.br.acc],
    idx.end = v.idx[i.br.acc + 1]
  )
  
  # Filter based on blocks
  df <- df[blocks.acc[abs(df$val.beg)] == blocks.acc[abs(df$val.end)],]
  
  # Add attributes: accuracy and lengths
  df$acc <- acc
  df$len.acc <- abs(df$val.end - df$val.beg) - 1
  df$len.comb <- abs(df$idx.end - df$idx.beg) - 1
  
  return(df)
}
