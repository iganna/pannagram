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
  
  # ---- Merge coverage ----
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


mergeOverlapsTolerance <- function(breaks, 
                                   len.tol = 0.95,
                                   dist.tol = 0.5,
                                   member.tol = 0.5){
  
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
  
  # ---- Merge coverage ----
  n.init = nrow(breaks)
  # Keep the untouched per-accession rows BEFORE the columns are dropped: after the
  # overlap merge below, len.acc holds only the MAX over the accessions of a merged
  # interval. That max is the right size to RESERVE, but a poor statistic to COMPARE
  # two breaks by -- one accession whose insertion is a few bp longer than its
  # neighbours' shifts the ratio and can push a genuine pair below len.tol. We
  # recompute a robust representative length from these rows once intervals are final,
  # and use the accession column to tell a real second event from alignment noise.
  has.acc = 'acc' %in% colnames(breaks)
  breaks.rows = breaks[, c('idx.beg', 'idx.end', 'len.acc',
                           if (has.acc) 'acc'), drop = FALSE]
  breaks = breaks[,c('idx.beg', 'idx.end', 'len.acc')]

  breaks <- breaks[order(breaks$idx.beg,
                         -breaks$idx.end,
                         -breaks$len.acc), ]
  
  breaks$id = 1:nrow(breaks)
  
  # Remove duplicates
  idx.keep = c(T, diff(breaks$idx.beg) != 0 | diff(breaks$idx.end) != 0)
  breaks = breaks[idx.keep,]
  
  breaks$cnt = c(breaks$id[-1], n.init+1) - breaks$id
  if(!(sum(breaks$cnt) == n.init)) stop('Checkpoint1')
  if (any(breaks$cnt < 0)) stop('Checkpoint2')
  
  
  # Merge overlaps
  n = 0
  while (n != nrow(breaks)) {
    n = nrow(breaks)
    
    idx_full_cover = which(breaks$idx.beg[-1] <= breaks$idx.end[-nrow(breaks)])
    idx_full_cover = setdiff(idx_full_cover, idx_full_cover + 1)
    
    if (length(idx_full_cover) == 0) break
    
    breaks$cnt[idx_full_cover] = breaks$cnt[idx_full_cover] + breaks$cnt[idx_full_cover + 1]
    
    breaks$idx.end[idx_full_cover] = pmax(breaks$idx.end[idx_full_cover], 
                                          breaks$idx.end[idx_full_cover + 1])
    
    breaks$len.acc[idx_full_cover] = pmax(breaks$len.acc[idx_full_cover],
                                          breaks$len.acc[idx_full_cover + 1])
    
    breaks = breaks[-(idx_full_cover + 1), ]
  }
  
  # ---- Representative length per (merged) interval ----
  # Assign every original row to the interval that swallowed it, drop the small
  # indels that share the window with a big insertion (< half of the max), and
  # take the median of what is left. On a clean interval this equals the max; on
  # a ragged one it ignores the single longest outlier.
  k = findInterval(breaks.rows$idx.beg, breaks$idx.beg)
  k[k < 1] = 1
  grp = split(breaks.rows$len.acc, k)
  rep.len = vapply(grp, function(v){
    v = v[v >= 0.5 * max(v)]
    stats::median(v)
  }, numeric(1))
  breaks$len.rep = breaks$len.acc
  pos = as.integer(names(rep.len))
  ok = pos >= 1 & pos <= nrow(breaks)
  breaks$len.rep[pos[ok]] = rep.len[ok]
  
  # ---- Accessions that substantially belong to each interval ----
  # A break may list an accession that contributes only a couple of nucleotides
  # there while all of its actual sequence sits in the neighbouring break. That is
  # alignment noise, not evidence of two independent events, so membership is
  # counted only above member.tol of the interval's representative length.
  acc.sets = vector('list', nrow(breaks))
  if (has.acc) {
    rows.by = split(seq_len(nrow(breaks.rows)), k)
    for (nm in names(rows.by)) {
      pp = as.integer(nm)
      if (pp < 1 || pp > nrow(breaks)) next
      rr = rows.by[[nm]]
      keep = breaks.rows$len.acc[rr] >= member.tol * breaks$len.rep[pp]
      acc.sets[[pp]] = unique(breaks.rows$acc[rr][keep])
    }
  }
  
  len.max.sv = 15000
  gap.max = 5000
  
  breaks$merge = 0
  if(nrow(breaks) > 1){
    for(i in 1:(nrow(breaks)-1)){
      if(breaks$len.acc[i] > len.max.sv)  next
      
      for(j in (i+1):nrow(breaks)){
        len.sim = min(breaks$len.rep[i], breaks$len.rep[j]) / 
          max(breaks$len.rep[i], breaks$len.rep[j])
        
        gap =  (breaks$idx.beg[j] - breaks$idx.end[i] + 1)
        if(gap > gap.max) break
        
        dist.sim = gap / max(breaks$len.acc[i], breaks$len.acc[j])
        if(dist.sim > dist.tol) break
        
        if(len.sim >= len.tol){
          # Refuse to merge when one accession substantially occupies more than one
          # of the intervals being joined: that is two real events, not one split.
          if (!has.acc || !any(duplicated(unlist(acc.sets[i:j])))) {
            breaks$merge[i:j] = 1
          }
        }
      }
    }  
  }
  
  
  if(sum(breaks$merge) > 0){
    df.merges = findOnes((breaks$merge > 0) * 1)
    breaks.add = data.frame(matrix(0, nrow = nrow(df.merges), 
                                   ncol = ncol(breaks),
                                   dimnames = list(NULL, colnames(breaks))))
    for(irow in 1:nrow(df.merges)){
      breaks.add$idx.beg[irow] = breaks$idx.beg[df.merges$beg[irow]]
      breaks.add$idx.end[irow] = breaks$idx.end[df.merges$end[irow]]
      breaks.add$len.acc[irow] = max(breaks$len.acc[df.merges$beg[irow]:
                                                      df.merges$end[irow]])
      breaks.add$len.rep[irow] = max(breaks$len.rep[df.merges$beg[irow]:
                                                      df.merges$end[irow]])
      breaks.add$cnt[irow] = sum(breaks$cnt[df.merges$beg[irow]:
                                                  df.merges$end[irow]])
      breaks.add$id = breaks.add$idx.beg
      
      breaks$merge[df.merges$beg[irow]:df.merges$end[irow]] = irow
    }
    breaks = rbind(breaks[breaks$merge == 0,], breaks.add)
    
    breaks <- breaks[order(breaks$idx.beg,
                           -breaks$idx.end,
                           -breaks$len.acc), ]
  }
  
  
  if(is.unsorted(breaks$idx.beg)) stop('Chrckpoint sorted 1')
  if(is.unsorted(breaks$idx.end)) stop('Chrckpoint sorted 2')
  
  breaks$len = breaks$idx.end - breaks$idx.beg + 1
  
  return(breaks)
}


solveLong <- function(breaks, breaks.init, len.max) {
  # Identify breaks where the length exceeds the len.max
  idx.to.solve = which(breaks$len > len.max)
  
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
      if (max(breaks.tmp.m$len) < len.max) break
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
      if ((p.e - p.b + 1) < len.max) {
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



solveLong2 <- function(breaks, breaks.init, len.max) {
  # Identify breaks where the length exceeds the len.max
  idx.to.solve = which(breaks$len > len.max)
  
  breaks.init$id <- 1:nrow(breaks.init)
  
  idx.rem.init <- c()
  for (irow in idx.to.solve) {
    pokaz(irow)
    pos.b <- breaks$idx.beg[irow]
    pos.e <- breaks$idx.end[irow]
    n <- pos.e - pos.b + 1
    
    # Select initial breaks that fit within the current merged break
    breaks.tmp <- breaks.init[(breaks.init$idx.beg >= pos.b) & (breaks.init$idx.end <= pos.e),]
    
    # Remove some breaks pased on idx
    breaks.tmp.len = pos.e - pos.b + 1
    pos.cov = rep(0, breaks.tmp.len)
    for(irow in 1:nrow(breaks.tmp)){
      irow.idx = (breaks.tmp$idx.beg[irow]:breaks.tmp$idx.end[irow]) - pos.b + 1
      pos.cov[irow.idx] = pos.cov[irow.idx] + 1
    }
    
    removed = removePositionsOverlap(pos.cov, len.max)
    idx.rem = c()
    for(i.pos in removed){
      pos.remove = i.pos - 1 + pos.b  
      idx.rem = c(idx.rem, which((breaks.tmp$idx.beg <= pos.remove) & (breaks.tmp$idx.end >= pos.remove)))
    }
    idx.rem = unique(idx.rem)
    idx.rem.init = c(idx.rem.init, breaks.tmp[idx.rem,]$id)
    
  }
  
  return(idx.rem.init)
}

removePositionsOverlap  <- function(pos.cov, len.max) {
  # Remove minimum of sequences to split into streaches of 25000
  n <- length(pos.cov)
  
  # a[1] = dummy 0 on the left
  # a[2..n+1] = original array
  # a[n+2] = dummy 0 on the right
  a <- c(0, pos.cov, 0)
  
  dp <- rep(0, n + 2)
  prev_pos <- rep(-1, n + 2)
  
  # deque implemented with arrays head..tail, storing indices of positions in dp
  q <- integer(n + 2)
  head <- 1
  tail <- 1
  q[1] <- 1   # index of the left dummy position
  
  # i = 2..n+2 corresponds to original positions 1..n and the right dummy position
  for (i in 2:(n + 2)) {
    # remove positions that are too far:
    # (i-1) - (j-1) > len.max  <=>  j < i - len.max
    while (head <= tail && q[head] < i - len.max) {
      head <- head + 1
    }
    
    best_j <- q[head]
    dp[i] <- dp[best_j] + a[i]
    prev_pos[i] <- best_j
    
    while (head <= tail && dp[q[tail]] >= dp[i]) {
      tail <- tail - 1
    }
    tail <- tail + 1
    q[tail] <- i
  }
  
  # reconstruct removed positions
  removed <- integer()
  cur <- n + 2
  while (cur != 1) {
    if (cur != n + 2) {
      # convert index from R-model back to array position: cur -> cur - 1
      removed <- c(cur - 1, removed)
    }
    cur <- prev_pos[cur]
  }

  return(removed)
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

  # A fully-absent accession (all zeros) has no breaks.
  if(length(v) == 0){
    return(data.frame(
      val.beg = numeric(0), val.end = numeric(0),
      idx.beg = numeric(0), idx.end = numeric(0),
      acc     = character(0), len.acc = numeric(0), len.comb = numeric(0)
    ))
  }

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
  
  # Initialize an array for block accumulation (indexed by abs() below, so it
  # must be sized by the max ABSOLUTE value, not the signed max).
  blocks.acc <- rep(0, max(abs(v)))
  for(irow in 1:nrow(v.b)) {
    blocks.acc[abs(v.b$v.beg[irow]):abs(v.b$v.end[irow])] <- irow
  }
  
  # Identify breaks where consecutive values differ by more than 1
  i.br.acc <- which(abs(diff(v)) != 1)
  if(length(i.br.acc) == 0){
    df <- data.frame(
      val.beg  = numeric(0),
      val.end  = numeric(0),
      idx.beg  = numeric(0),
      idx.end  = numeric(0),
      acc      = character(0),
      len.acc  = numeric(0),
      len.comb = numeric(0)
    )
    
  } else {
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
  }
 
  return(df)
}
