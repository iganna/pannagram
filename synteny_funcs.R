#' Determine Small Overlaps
#'
#' This function identifies small overlaps between rows in a data frame.
#' It compares adjacent rows based on their overlap and 
#' removes the smaller overlap based on column 'V7' (length).
#' The function adds a new column 'rm.len' which indicates the length of the overlap removed.
#'
#' @param x.df A data frame with alignment.
#'               Expected to be sorted based on 'p.beg'.
#'
#' @return A data frame similar to the input but with column 'rm.len' added.
#
defineOverlapps <- function(x.df){
  
  # sort.flaf = isSorted(x.df$p.beg)
  # new.col.name = 'tmp'
  # if(!sort.flaf){
  #   x.df[,new.col.name] = 1:nrow(x.df)
  #   x.df[order(x.df$p.beg),]
  # }
  
  if(sum(c('p.beg', 'p.end') %in% colnames(x.df)) != 2){
    x.df$p.beg <- ifelse(x.df$V4 < x.df$V5, x.df$V4, x.df$V5)
    x.df$p.end <- ifelse(x.df$V4 < x.df$V5, x.df$V5, x.df$V4)
  }
  
  
  # if not sorted - SORT!
  if(is.unsorted(x.df$p.beg)){
    # Sort by the position in the reference
    x.df = x.df[order(-x.df$V7),]
    x.df = x.df[order(x.df$p.beg),]
  }
  
  
  
  x.df$rm.len = 0
  idx.overlap = which(x.df$p.beg[-1] <= x.df$p.end[-nrow(x.df)])
  
  # x.df[c(irow, irow+1),]
  
  for(irow in idx.overlap){
    # We cut either irow or [irow+1]
    
    # Which row to cut
    icut = ifelse(x.df$V7[irow] > x.df$V7[irow + 1], irow + 1, irow)
    ibig = ifelse(x.df$V7[irow] > x.df$V7[irow + 1], irow, irow + 1)
    
    # Find left(value = 1) or right(value = -1) tails 
    # of irow(i.e., tail.irow) or [irow+1](i.e., tail.next), 
    # which are involved in the overlap
    tail.iсut = ifelse((x.df$V4[icut] >= x.df$p.beg[ibig]) & 
                         (x.df$V4[icut] <= x.df$p.end[ibig]), 
                       1, -1)
    
    # how much to cut
    ncut = length(intersect(x.df$V4[irow]:x.df$V5[irow],
                            x.df$V4[irow+1]:x.df$V5[irow+1]))
    # Remember the cut
    x.df$rm.len[icut] = ncut * tail.iсut
  }
  
  # if(!sort.flaf){
  #   x.df[order(x.df[,new.col.name]),]
  #   x.df = x.df[,colnames(x.df) != new.col.name]
  # }

  
  return(x.df)
}


cleanBigOverlaps <- function(x.df, rm.threshold = 0.5){
  
  while(T){
    n.row = nrow(x.df)
    x.df = defineOverlapps(x.df)
    x.df =  x.df[((abs(x.df$rm.len) / x.df$V7) <= rm.threshold) | (x.df$rm.len == 0),]
    if(n.row == nrow(x.df)) break
    pokaz(n.row)
  }
  return(x.df)
}

#' Remove Small Overlaps
#'
#' The function identifies and processes peviously defined overlaps
#'
#' @param x.sk A data frame with the alignment
#' @param rm.threshold A threshold value to decide which rows to process. 
#' Rows which should loose more than rm.threshold of their lengths - are removed
#'
#' @return A modified data frame based on the overlap criteria.
#' 
cutSmallOverlaps <- function(x.sk){
  
  # if(!isSorted(x.sk$p.beg)) pokazAttention('2!!')
  
  idx.cut = which(x.sk$rm.len != 0)
  for(irow in idx.cut){
    
    # if(x.sk$dir[irow] == 1) stop('dir')
    
    adjustment <- x.sk$rm.len[irow]  # sing says about "begin" or "end"
    adjustment.dir <- ifelse(x.sk$dir[irow] == 0, adjustment, (-1) * adjustment)
    n.symbols = abs(adjustment)
    
    if(adjustment > 0){
      # stop('3')
      # From the befinning!
      
      # Adjust strings in the alignment
      seq = seq2nt(x.sk$V9[irow])
      aln.adjust = which(seq != '-')[n.symbols]
      x.sk$V8[irow] = substr(x.sk$V8[irow], (aln.adjust+1), nchar(x.sk$V8[irow]))
      x.sk$V9[irow] = substr(x.sk$V9[irow], (aln.adjust+1), nchar(x.sk$V9[irow]))
      
      # Adjust positions
      s.q.cut = seq2nt(x.sk$V8[irow])
      adjustment.q = sum(s.q.cut != '-') - 1
      x.sk$V2[irow] = x.sk$V3[irow] - adjustment.q
      
      # x.sk$V2[irow] = x.sk$V2[irow] + adjustment.q * sign(adjustment)
      x.sk$V4[irow] = x.sk$V4[irow] + adjustment.dir
      
    }
    if(adjustment < 0){
      # stop('4')
      # From the ending!
      
      # Adjust strings in the alignment
      seq = seq2nt(x.sk$V9[irow])
      aln.adjust = tail(which(seq != '-'), n=n.symbols)[1]
      x.sk$V8[irow] = substr(x.sk$V8[irow], 1, (aln.adjust-1))
      x.sk$V9[irow] = substr(x.sk$V9[irow], 1, (aln.adjust-1))
      
      # Adjust positions
     
      s.q.cut = seq2nt(x.sk$V8[irow])
      adjustment.q = sum(s.q.cut != '-') - 1
      x.sk$V3[irow] = x.sk$V2[irow] + adjustment.q
      
      # x.sk$V3[irow] = x.sk$V3[irow] + adjustment.q * sign(adjustment)
      x.sk$V5[irow] = x.sk$V5[irow] + adjustment.dir
    }
    
    x.sk$V7[irow] = nchar(x.sk$V8[irow])
    
    # print(x.sk[irow, -c(8,9)])
    # if(x.sk$V1[irow] == 'acc_10001|chr_1|part_157|780001') stop()
  }
  
  
  # Undate begin-eng positions
  x.sk$p.beg <- ifelse(x.sk$V4 < x.sk$V5, x.sk$V4, x.sk$V5)
  x.sk$p.end <- ifelse(x.sk$V4 < x.sk$V5, x.sk$V5, x.sk$V4)
  
  return(x.sk)
}



checkCorrespToGenome <- function(t, base.fas.fw, base.fas.bw, query.fas, k = 10){
  for(irow in 1:nrow(t)) {
    s1 = toupper(remainLastN(t[irow, 'V8'], k))
    s1 = gsub("\\-","",s1)
    
    if(nchar(s1) == 0){
      s2 = s1
    } else {
      s2 = toupper(paste0(query.fas[(-(nchar(s1)-1):0) + t[irow, 'V3']], collapse = ''))
    }
    
    
    if(s1 != s2 ){
      pokaz('Row', irow)
      pokaz(s1)
      pokaz(s2)
      stop('Checking the query not passed')
    }
    
    s1 = toupper(remainLastN(t[irow, 'V9'], k))
    s1 = gsub("\\-","",s1)
    if(t$dir[irow] == 0) {
      base.fas = base.fas.fw
    } else {
      base.fas = base.fas.bw
    }
    
    if(nchar(s1) == 0){
      s2 = s1
    } else {
      s2 = toupper(paste0(base.fas[(-(nchar(s1)-1):0) + t[irow, 'V5']], collapse = ''))
    }
    
    
    if(s1 != s2 ){
      pokaz('Row', irow)
      pokaz(s1)
      pokaz(s2)
      stop('Checking the base was not passed')
    }
  }
}


#' Traverse of blast results, to find the optimal path
#'
#' This function traverses a graph allowing the search in a last skipped window
#'
#' @param x.tmp A data frame containing graph information.
#' @param irow The current row index.
#' @param x.top The top x-coordinate (query)
#' @param y.top The top y-coordinate (base)
#' @param w.beg The beginning of the window (base)
#' @param w.end The end of the window (base)
#' @param vist.info A list containing traversal information.
#'
#' @return A list containing updated traversal information.
#'
graphTraverseWnd <- function(x.tmp, irow, x.top, y.top, w.beg, w.end, vist.info,
                             forward = F) {
  
  if(irow == nrow(x.tmp)) return(vist.info)
  
  # Find neighbours
  i.rem = -(1:irow)
  jrow = which(((x.tmp$V3[i.rem] - x.top) > 0) &
                 ((x.tmp$V5[i.rem] - y.top) > 0)) + irow
  
  x.over = x.top - x.tmp$V2[jrow]
  y.over = y.top - x.tmp$V4[jrow]
  max.over = pmax(ifelse(x.over > 0, x.over, 0), ifelse(y.over > 0, y.over, 0))
  w.pure = x.tmp$w[jrow] - 2 * max.over
  delta = (-x.over + max.over) + (-y.over + max.over)
  
  d = vist.info$d.to[irow] + delta - w.pure + max.over*2
  
  
  if((w.end -  w.beg) > 0){
    jrow.add =  which((x.tmp$V3[i.rem] > x.top) & 
                        (x.tmp$V5[i.rem] > w.beg) & 
                        (x.tmp$V4[i.rem] < w.end)) + irow 
    
    y.over = ifelse(x.tmp$V5[jrow.add] > w.end, 
                    x.tmp$V5[jrow.add] - w.end, 0) +  # up
      ifelse(x.tmp$V4[jrow.add] < w.beg, 
             w.beg - x.tmp$V4[jrow.add], 0)  # low
    x.over = x.top - x.tmp$V2[jrow.add]
    max.over = pmax(ifelse(x.over > 0, x.over, 0), y.over)
    w.pure = x.tmp$w[jrow.add] - 2 * max.over

    d.add = vist.info$d.to[irow] + (-x.over + max.over) - w.pure + max.over*2

    
    jrow = c(jrow, jrow.add)
    d = c(d, d.add)
  }
  
  if(length(jrow) == 0) return(vist.info)
  
  # Order neighbours
  order.neighbours <- order(d)
  d     <- d[order.neighbours]
  jrow  <- jrow[order.neighbours]
  
  for (j in 1:length(jrow)) {
    if(d[j] >= vist.info$d.to[jrow[j]]) next
    
    vist.info$d.to[jrow[j]] = d[j]
    vist.info$v.prev[jrow[j]] = irow
    
    if(x.tmp$V5[jrow[j]] > y.top){  # Продвижение вперед
      # New window
      w.beg.next = y.top
      w.end.next = x.tmp$V4[jrow[j]]
    } else {
      # Shrink previous window
      w.beg.next = w.beg
      w.end.next = x.tmp$V4[jrow[j]]
    }
    
    # Remove the window, if it's covered now
    if((w.end.next - w.beg.next < 250)){
      w.beg.next = 0
      w.end.next = 0
    }
    
    if(forward){
      w.beg.next = 0
      w.end.next = 0
    }
    
    vist.info = graphTraverseWnd(x.tmp, jrow[j], 
                                 x.tmp$V3[jrow[j]], 
                                 max(x.tmp$V5[jrow[j]], y.top), 
                                 w.beg.next, 
                                 w.end.next, 
                                 vist.info)
  }
  return(vist.info)
}


reconstructTraverse <- function(predecessor) {
  start = 1
  end = length(predecessor)
  path <- integer(0)
  while (!is.na(end)) {
    path <- c(end, path)
    if (end == start) break
    end <- predecessor[end]
  }
  return(path)
}


#' Count consecutive occurrences of a specified value in a vector.
#'
#' This function starts counting from the next position to a specified one (`irow`) in the vector (`vec`)
#' and continues until a different value is encountered or the end of the vector is reached.
#' By default, the function looks for consecutive occurrences of -1.
#'
#' @param vec A numeric or character vector.
#' @param irow Starting position in `vec` to begin counting from.
#' @param value The value to count consecutive occurrences of. Default is -1.
#' 
#' @return An integer representing the number of consecutive occurrences of `value`
#' starting from position `irow`.
#' @examples
#' vec <- c(1, 2, -1, -1, -1, 3, 4, -1)
#' countValueStretch(vec, 2) # 3
#' countValueStretch(vec, 1) # 0
#' countValueStretch(vec, 7) # 1
#' 
countValueStretch <- function(vec, irow, value = -1) {
  if(irow > length(vec)) stop('Wong position in vector, irow should be <= the length.')
  
  count <- 0
  if(irow == length(vec)) return(0)
  for (i in (irow + 1):length(vec)) {
    if (vec[i] == value) {
      count <- count + 1
    } else {
      break
    }
  }
  return(count)
}


#' Count consecutive occurrences of a specified value in a vector, moving backward.
#'
#' This function starts counting from the previous position to a specified one (`irow`) in the vector (`vec`)
#' and continues backward until a different value is encountered or the beginning of the vector is reached.
#' By default, the function looks for consecutive occurrences of -1.
#'
#' @param vec A numeric or character vector.
#' @param irow Starting position in `vec` to begin counting from.
#' @param value The value to count consecutive occurrences of. Default is -1.
#' 
#' @return An integer representing the number of consecutive occurrences of `value`
#' starting from position `irow` and moving backward.
#' @examples
#' vec <- c(1, 2, -1, -1, -1, 3, 4, -1)
#' countValueStretchBW(vec, 6) # 3
#' countValueStretchBW(vec, 8) # 0
#' countValueStretchBW(vec, 4) # 1
#' 
countValueStretchBW <- function(vec, irow, value = -1) {
  if(irow < 1) stop('Wong position in vector, irow should be >= 1')
  
  count <- 0
  if(irow == 1) return(0)
  for (i in (irow - 1):1) {
    if (vec[i] == value) {
      count <- count + 1
    } else {
      break
    }
  }
  return(count)
}


#' Initialize Visit Information for Graph Traversal in Bellman-Ford Algorithm
#'
#' This function prepares the initial structure required for the Bellman-Ford algorithm.
#' It initializes data structure with `v.prev` indicating the predecessor of each vertex,
#' and `d.to` indicating the distance to each vertex. At the start, `v.prev` is set to `NA`
#' for all vertices, indicating no predecessors, and `d.to` is set to `Inf` for all vertices,
#' indicating an infinite distance, except for the first vertex, which is set to 0.
#'
#' @param n An integer representing the number of vertices in the graph.
#' 
#' @return A list with two components:
#' \itemize{
#'   \item \code{v.prev}: A numeric vector of length `n`
#'                        representing the predecessor of each vertex.
#'   \item \code{d.to}: A numeric vector of length `n` representing 
#'                        distances to each vertex.
#' }
#' 
initVisitInfo <- function(n){
  if(n <= 1) stop('Number of visits should be positive')
  vist.info = list(v.prev = rep(NA, n),
                   d.to = rep(Inf, n))
  vist.info$d.to[1] = 0
  return(vist.info)
}

#' Remove Specified Columns from a Data Frame or Matrix
#'
#' @param x A data frame or matrix.
#' @param columnt.to.remove A numeric vector of column indices to remove. Default is c(8,9).
#' 
#' @return Modified data frame or matrix without the specified columns.
#' 
info <- function(x, columnt.to.remove = c(8,9)){
  return(x[,-columnt.to.remove])
}



pathDownPlus <- function(x.tmp){
  
  # Flip base
  b.max = max(c(x.tmp$V4, x.tmp$V5)) + 1
  x.tmp$V4 = b.max - x.tmp$V4
  x.tmp$V5 = b.max - x.tmp$V5
  
  # Flip query
  q.max = max(c(x.tmp$V2, x.tmp$V3)) + 1
  x.tmp$V2 = q.max - x.tmp$V2
  x.tmp$V3 = q.max - x.tmp$V3
  
  # DON'T SORT
  # PLEASE NO SORTING AND NO CHANGES OF V4-V5
  
  # Get the best path
  idx.visit = pathUpPlus(x.tmp)
  
  return(idx.visit)
}

pathDownMinus <- function(x.tmp){
  
  # Flip query
  q.max = max(c(x.tmp$V2, x.tmp$V3)) + 1
  x.tmp$V2 = q.max - x.tmp$V2
  x.tmp$V3 = q.max - x.tmp$V3
  
  # DON'T SORT
  # PLEASE NO SORTING AND NO CHANGES OF V4-V5
  
  # Get the best path
  idx.visit = pathUpPlus(x.tmp)
  
  return(idx.visit)
}


pathUpMinus <- function(x.tmp, irow, visit.info){
  
  # Flip base
  b.max = max(c(x.tmp$V4, x.tmp$V5)) + 1
  x.tmp$V4 = b.max - x.tmp$V4
  x.tmp$V5 = b.max - x.tmp$V5
  
  # DON'T SORT
  # PLEASE NO SORTING AND NO CHANGES OF V4-V5
  
  idx.visit = pathUpPlus(x.tmp)
  return(idx.visit)
}


pathUpPlus <- function(x.tmp){
  
  # shave the initial order
  x.tmp$ord = 1:nrow(x.tmp)
  
  # Flip begin-end points
  idx.flip = x.tmp$V2 > x.tmp$V3
  tmp = x.tmp$V2[idx.flip]
  x.tmp$V2[idx.flip] = x.tmp$V3[idx.flip]
  x.tmp$V3[idx.flip] = tmp
  
  tmp = x.tmp$V4[idx.flip]
  x.tmp$V4[idx.flip] = x.tmp$V5[idx.flip]
  x.tmp$V5[idx.flip] = tmp
  
  x.tmp = x.tmp[order(x.tmp$V2),]
  
  # Add final destination
  x.tmp = rbind(x.tmp, 0)
  n.tmp = nrow(x.tmp)
  
  x.max = max(c(x.tmp$V2, x.tmp$V3)) + 1
  y.max = max(c(x.tmp$V4, x.tmp$V5)) + 1
  
  x.tmp$V2[n.tmp] = x.max
  x.tmp$V3[n.tmp] = x.max
  x.tmp$V4[n.tmp] = y.max
  x.tmp$V5[n.tmp] = y.max
  
  # Set up direction
  idx.dir = which(x.tmp$V5 < x.tmp$V4)
  tmp = x.tmp$V4[idx.dir]
  x.tmp$V4[idx.dir] = x.tmp$V5[idx.dir]
  x.tmp$V5[idx.dir] = tmp
  x.tmp$dir = 0
  x.tmp$dir[idx.dir] = 1
   
  # Run the traverse
  visit.info = initVisitInfo(nrow(x.tmp))
  visit.info = graphTraverseWnd(x.tmp, 1, x.tmp$V3[1], x.tmp$V5[1], 0, 0, visit.info)
  idx.visit = reconstructTraverse(visit.info$v.prev)  # Reconstruct the optimal path
  idx.visit = idx.visit[-c(1, length(idx.visit))]
  # idx.visit = idx.visit[-c(length(idx.visit))]  # Remain the anchor

  
  # Return previous order
  idx.visit = x.tmp$ord[idx.visit]
  
  return(idx.visit)
}


# traverseUpPlus <- function(x.tmp, irow, visit.info) {
#   
#   if(irow == nrow(x.tmp)) return(visit.info)
#   # Find neighbours
#   jrow = which((x.tmp$V3[-(1:irow)] > x.tmp$V3[irow]) & 
#                  (x.tmp$V5[-(1:irow)] > x.tmp$V5[irow])) + irow 
#   if(length(jrow) == 0) return(visit.info)
#   
#   d = abs(x.tmp$V2[jrow] -  x.tmp$V3[irow]) + abs(x.tmp$V4[jrow] -  x.tmp$V5[irow])
#   # d = pmin(abs(x.tmp$V2[jrow] -  x.tmp$V3[irow]), abs(x.tmp$V4[jrow] -  x.tmp$V5[irow])) * 2
#   d = d + visit.info$d.to[irow] - x.tmp$w[jrow]
#   
#   # Order neighbours
#   order.neighbours <- order(d)
#   d     <- d[order.neighbours]
#   jrow  <- jrow[order.neighbours]
#   
#   # # Remove those, which are globally not advantageous
#   # idx.stretch = d < 0
#   # if(sum(idx.stretch) == 0) return(visit.info)
#   # jrow <- jrow[idx.stretch]
#   # d <- d[idx.stretch]
#   
#   for (i in 1:length(jrow)) {
#     if(d[i] < visit.info$d.to[jrow[i]]){
#       visit.info$d.to[jrow[i]] = d[i]
#       visit.info$v.prev[jrow[i]] = irow
#       visit.info = traverseUpPlus(x.tmp, jrow[i], visit.info)
#     }
#   }
#   return(visit.info)
# }

# reconstructUpPlus <- function(predecessor) {
#   start = 1
#   end = length(predecessor)
#   while(is.na(predecessor[end])){
#     if (end == start) return(NULL)
#     end = end - 1
#   }
#   path <- integer(0)
#   while (!is.na(end)) {
#     path <- c(end, path)
#     if (end == start) break
#     end <- predecessor[end]
#   }
#   return(path)
# }


# ----  OLD FUNCTIONS  ----

getCorresponding2BasePositions <- function(t, base.len){
  t.base = getBase(t, base.len)
  pos.corresp = rep(0, base.len)
  for(irow in 1:nrow(t)){
    # if(echo) print(irow)
    positions.q = rep(0, t$V7[irow])
    positions.q[strsplit(t$V8[irow],'')[[1]] != '-'] = t$V2[irow]:t$V3[irow]
    positions.b = rep(0, t$V7[irow])
    if(t$dir[irow] == 0){
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V4[irow]:t.base$V5[irow]  
    } else {
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V5[irow]:t.base$V4[irow]
    }
    
    
    positions.q = positions.q[positions.b != 0]
    positions.b = positions.b[positions.b != 0]
    pos.corresp[positions.b] = positions.q
    # 
    # positions.q.all = c(positions.q.all, positions.q)
    # positions.b.all = c(positions.b.all, positions.b)
  }
  return(pos.corresp)
}

getCorresponding2BasePositionsSign <- function(t, base.len){
  t.base = getBase(t, base.len)
  pos.corresp = rep(0, base.len)
  for(irow in 1:nrow(t)){
    # if(echo) print(irow)
    positions.q = rep(0, t$V7[irow])
    positions.q[strsplit(t$V8[irow],'')[[1]] != '-'] = t$V2[irow]:t$V3[irow]
    positions.b = rep(0, t$V7[irow])
    if(t$dir[irow] == 0){
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V4[irow]:t.base$V5[irow]  
    } else {
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = -(t.base$V5[irow]:t.base$V4[irow])
    }
    
    
    positions.q = positions.q[positions.b != 0]
    positions.b = positions.b[positions.b != 0]
    pos.corresp[abs(positions.b)] = positions.q * sign(positions.b)
    # 
    # positions.q.all = c(positions.q.all, positions.q)
    # positions.b.all = c(positions.b.all, positions.b)
  }
  return(pos.corresp)
}


getCorresponding2QueryPositions <- function(t, base.len, query.len){
  t.base = getBase(t, base.len)
  pos.corresp = rep(0, query.len)
  for(irow in 1:nrow(t)){
    # if(echo) print(irow)
    positions.q = rep(0, t$V7[irow])
    positions.q[strsplit(t$V8[irow],'')[[1]] != '-'] = t$V2[irow]:t$V3[irow]
    positions.b = rep(0, t$V7[irow])
    if(t$dir[irow] == 0){
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V4[irow]:t.base$V5[irow]  
    } else {
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V5[irow]:t.base$V4[irow]
    }
    
    
    
    positions.b = positions.b[positions.q != 0]
    positions.q = positions.q[positions.q != 0]
    pos.corresp[positions.q] = positions.b
    # 
    # positions.q.all = c(positions.q.all, positions.q)
    # positions.b.all = c(positions.b.all, positions.b)
  }
  return(pos.corresp)
}


getPositionQueryByPositionBaseCorresp <- function(t, pos1, pos2, base.len, query.len, echo=F){
  
  t.base = getBase(t, base.len)
  
  pos.corresp = rep(0, base.len)
  for(irow in 1:nrow(t)){
    if(echo) print(irow)
    positions.q = rep(0, t$V7[irow])
    positions.q[strsplit(t$V8[irow],'')[[1]] != '-'] = t$V2[irow]:t$V3[irow]
    positions.b = rep(0, t$V7[irow])
    if(t$dir[irow] == 0){
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V4[irow]:t.base$V5[irow]  
    } else {
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V5[irow]:t.base$V4[irow]
    }
    
    idx = positions.b != 0
    positions.b = positions.b[idx]
    positions.q = positions.q[idx]
    
    pos.corresp[positions.b] = positions.q
  }
  print('ok')
  # 
  idx = pos1 > pos2
  tmp = pos1[idx]
  pos1[idx] = pos2[idx]
  pos2[idx] = tmp
  # 
  # pos.tmp = -Inf
  # pos.prev = pos.corresp
  # for(i in 1:length(pos.corresp)){
  #   if(pos.corresp[i] == 0) {
  #     pos.prev[i] = pos.tmp
  #   } else {
  #     pos.tmp = max(pos.prev[i], pos.tmp)
  #   }
  # }
  # pos.tmp = +Inf
  # pos.next = pos.corresp
  # for(i in length(pos.corresp):1){
  #   if(pos.corresp[i] == 0) {
  #     pos.next[i] = pos.tmp
  #   } else {
  #     pos.tmp = min(pos.next[i], pos.tmp)
  #   }
  # }
  
  pos.base = cbind(pos.corresp[pos1], pos.corresp[pos2])
  # pos.flanking = cbind(pos.prev[pos1], pos.next[pos2])
  # 
  # return(cbind(pos.base, pos.flanking))
  
  return(pos.base)
}


removeOverlapdAfterBlastGap <- function(t){
  unique_names = unique(t$V1)
  t_new = c()
  rownames(t) <- NULL
  for(name in unique_names){
    #print(name)
    t_tmp = t[t$V1 == name,]
    t_tmp[,c('V2', 'V3')] = t_tmp[,c('V2', 'V3')] - min(t_tmp$V2) + 1
    n.max = max(t_tmp$V3)
    if(nrow(t_tmp) > 1){
      t_tmp = t_tmp[order(-t_tmp$V6),]
      cover = t_tmp$V3 - t_tmp$V2 + 1
      t_tmp = t_tmp[order(-cover),]
      cover = t_tmp$V3 - t_tmp$V2 + 1
      
      pos = rep(0, n.max)
      
      idx_remove = c()
      pos[t_tmp$V2[1]:t_tmp$V3[1]] = 1
      for(irow in 2:nrow(t_tmp)){
        if(sum(pos[t_tmp$V2[irow]:t_tmp$V3[irow]]) > 0.85 * cover[irow]){
          idx_remove = c(idx_remove, irow)
        } else {
          pos[t_tmp$V2[irow]:t_tmp$V3[irow]] = 1
        }
      }
      # if(length(idx_remove) > 0) t_tmp = t_tmp[-idx_remove,]
      idx_remain = setdiff(1:nrow(t_tmp), idx_remove)
      for(irow in idx_remove){
        idx = ((t_tmp$V2[irow] %in% t_tmp$V2[idx_remain]) & 
                 (t_tmp$V3[irow] %in% t_tmp$V3[idx_remain]) &
                 (t_tmp$V6[irow] %in% t_tmp$V6[idx_remain]))
        if(sum(idx) > 0) idx_remain = c(idx_remain, irow)
      }
      t_tmp = t_tmp[idx_remain,]
    }
    t_new = rbind(t_new, t[rownames(t_tmp),])
  }
  return(t_new)
}


removeLastN <- function(x, n){
  substr(x, 1, nchar(x)-n)
}

removeFirstN <- function(x, n){
  substr(x, n+1, nchar(x))
}

remainLastN <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

remainFirstN <- function(x, n){
  substr(x, 1, n)
}

sumGapLastN <- function(x, n){
  s <- remainLastN(x, n)
  return(sum(strsplit(s, '')[[1]] == '-'))
}

sumGapFirstN <- function(x, n){
  s <- substr(x, 1, n)
  return(sum(strsplit(s, '')[[1]] == '-'))
}


fixDirection <- function(t, base.fas.bw){
  t$dir = c()
  for(irow in 1:nrow(t)) {
    if(t[irow,'V5'] < t[irow,'V4']) {
      
      pos1 = t[irow,'V5']
      pos2 = t[irow,'V4']
      
      t[irow,'V5'] <- length(base.fas.bw) - pos1 + 1
      t[irow,'V4'] <- length(base.fas.bw) - pos2 + 1
      t[irow, 'dir'] = 1
    } else {
      t[irow, 'dir'] = 0
    }
  }
  return(t)
}

showt <- function(t, irow=NULL, chr=F, add = c()){
  idx = c(1:5,7,11, add)
  idx = intersect(idx, 1:ncol(t))
  irow = irow[irow <= nrow(t)]
  
  if(chr) idx = c(idx, 10)
  if(is.null(irow)){
    print(t[1:min(100, nrow(t)),idx])
  } else {
    print(t[irow,idx])
  }
}

glueByThreshold <- function(t, thresholds, base.fas.fw, base.fas.bw, query.fas.chr,
                            file.log = NULL, gap.open = 30, gap.ext = 0.5, maxval = 10^10,
                            base.overlap = T, query.overlap = T, log.append = F) {
  
  if(!is.null(file.log)) {
    write('', file=file.log, append=log.append)
  }
  t = orderT(t)
  for(threshold in thresholds) {
    
    irow = 1
    while(irow < nrow(t)) {
      # print(threshold)
      # print(irow)
      
      # if a base fragment of irow is already within some base fragment
      if(base.overlap){
        idx = (t[,'V4'] <= t[irow,'V4']) & (t[,'V5'] >= t[irow,'V5'])
        idx[t$dir != t$dir[irow]] = F
        idx[irow] = F
        if (sum(idx) > 0) {
          irow <- irow + 1
          next
        }        
      }

      
      # if a query fragment of irow is already within some query fragment
      if(query.overlap){
        idx = (t[,'V2'] <= t[irow,'V2']) & (t[,'V3'] >= t[irow,'V3'])
        idx[t$dir != t$dir[irow]] = F
        idx[irow] = F
        if (sum(idx) > 0) {
          irow <- irow + 1
          next
        }  
      }
      
      
      idx = abs(t[,'V2'] - (t[irow,'V3']))
      # idx = apply(cbind(abs(t[,'V2'] - (t[irow,'V3'])), abs(t[,'V4'] - (t[irow,'V5']))), 1, min)
      idx[irow] = maxval
      idx[t[,'V2'] <= t[irow,'V2']] = maxval
      idx[t[,'V4'] <= t[irow,'V4']] = maxval
      
      # Chromosomes should match
      idx[t[,'V10'] != t[irow,'V10']] = maxval
      
      # If new fragment is in irow fragment by query
      idx.inside = (t[,'V2'] >= t[irow,'V2']) & (t[,'V3'] <= t[irow,'V3'])
      idx[idx.inside] = maxval
      
      # If new fragment is in irow fragment by base
      idx.inside = (t[,'V4'] >= t[irow,'V4']) & (t[,'V5'] <= t[irow,'V5'])
      idx[idx.inside] = maxval
      
      idx[t$dir != t$dir[irow]] <- maxval
      
      idx[abs(t[, 'V2'] - t[irow, 'V3']) > threshold] = maxval
      idx[abs(t[, 'V4'] - t[irow, 'V5']) > threshold] = maxval
      
      if (min(idx) == maxval) {
        irow <- irow + 1
        next
      }
      
      idx = which(idx == min(idx))
      
      if (length(idx) == 0) {
        irow <- irow + 1
        next
      }
      idx = idx[1]
      
      if(threshold > 10^4){
        n.char = 50
        n.char.thresh = 100  
      } else {
        n.char = 9
        n.char.thresh = 50
      }
      
      n.char.thresh = min(min(n.char.thresh, (t[irow, 'V7'] - 2)), (t[idx, 'V7'] - 2) )
      
      d.tmp = 0
      while (d.tmp < n.char.thresh) {
        n.char = n.char + 1
        pos.query.start = max(1, t[irow, 'V3'] - n.char + sumGapLastN(t[irow, 'V8'], n.char))
        pos.query.end = t[idx, 'V2'] + n.char - sumGapFirstN(t[idx, 'V8'], n.char)
        
        
        pos.base.start = max(1, t[irow, 'V5'] - n.char + sumGapLastN(t[irow, 'V9'], n.char))
        pos.base.end = t[idx, 'V4'] + n.char - sumGapFirstN(t[idx, 'V9'], n.char)
        
        if(pos.base.end >= t[idx, 'V5']) break
        if(pos.query.end >= t[idx, 'V3']) break
        if(pos.base.start <= t[irow, 'V4']) break
        if(pos.query.end <= t[irow, 'V2']) break
        
        d.tmp = min(pos.query.end - pos.query.start, pos.base.end - pos.base.start)
      }
      
      if(t$dir[irow] == 0){
        base.fas = base.fas.fw
      } else {
        base.fas = base.fas.bw
      }
      
      
      seq.query <- toupper(query.fas.chr[(pos.query.start+1):(pos.query.end-1)])
      seq.base <- toupper(base.fas[(pos.base.start+1):(pos.base.end-1)])  
      
      globalAlign<- pairwiseAlignment(paste0(seq.query, collapse=''), 
                                      paste0(seq.base, collapse=''), gapOpening = gap.open, 
                                      gapExtension=gap.ext,
                                      type="global")
      globalAlign 
      
      if(!is.null(file.log)) {
        write(paste((pos.query.start+1), (pos.query.end-1), as.character(alignedPattern(globalAlign))),
              file=file.log, append=TRUE)
        write(paste((pos.base.start+1), (pos.base.end-1), as.character(alignedSubject(globalAlign))), 
              file=file.log, append=TRUE)
        write('\n', file=file.log, append=TRUE)
      }
      
      
      seq.query.new = paste(removeLastN(t[irow, 'V8'], n.char),  
                            alignedPattern(globalAlign), 
                            removeFirstN(t[idx, 'V8'], n.char), sep = '')
      
      seq.base.new = paste(removeLastN(t[irow, 'V9'], n.char),  
                           alignedSubject(globalAlign), 
                           removeFirstN(t[idx, 'V9'], n.char), sep = '')
      
      if(nchar(seq.query.new) != nchar(seq.base.new)) stop('aa')
      
      t1 <- t
      
      t[irow, 'V8'] <- seq.query.new
      t[irow, 'V9'] <- seq.base.new
      
      t[irow, 'V3'] <- t[idx, 'V3']
      t[irow, 'V5'] <- t[idx, 'V5']
      t[irow, 'V7'] <- nchar(seq.base.new)
      
      t <- t[-idx,]
      
      # s1 = seq.base.new
      # s2 = paste0(base.fas.fw[t[irow,'V4']:t[irow,'V5']], collapse = '')
      
      s.base = strsplit(seq.base.new,'')[[1]]
      idx.nogap.base = which(s.base != '-')
      if(length(idx.nogap.base) != abs(t[irow,'V5'] - t[irow,'V4']) + 1) stop('wrong in base')
      
      s.query = strsplit(seq.query.new,'')[[1]]
      idx.nogap.query = which(s.query != '-')
      if(length(idx.nogap.query) != abs(t[irow,'V3'] - t[irow,'V2']) + 1) stop('wrong in query')
      
      
      if(idx < irow) {
        irow = irow - 1
      }    
      
      rownames(t) <- NULL
    }
  }
  return(t)
}


glueZero <- function(t){
  t = t[order(t[,'V2']),]
  
  ipos = 1
  while(ipos < nrow(t)) {
    # print(ipos)
    idx = which((t[,'V4'] == (t[ipos,'V5'] + 1)) & (t[,'V2'] == (t[ipos,'V3'] + 1)))
    if (length(idx) == 0){
      ipos = ipos+1
      next
    }
    # if(idx != (ipos+1)) print(idx)  # glue with not the next record
    
    t[ipos, 'V3'] <- t[idx, 'V3']
    t[ipos, 'V5'] <- t[idx, 'V5']
    t[ipos, 'V7'] <- t[ipos, 'V7'] + t[idx, 'V7']
    t[ipos, 'V8'] <- paste(t[ipos, 'V8'], t[idx, 'V8'], sep = '')
    t[ipos, 'V9'] <- paste(t[ipos, 'V9'], t[idx, 'V9'], sep = '')
    t <- t[-idx,]
  }
  
  rownames(t) <- NULL
  return(t)
}

glueAlmostZero <- function(t, threshold, base.fas.fw, base.fas.bw, query.fas.chr){
  t = t[order(t[,'V2']),]
  rownames(t) <- NULL
  ipos = 1
  while(ipos < nrow(t)) {
    
    if(t$dir[ipos] != t$dir[ipos+1]){
      ipos = ipos + 1
      next
    }
    # print(ipos)
    d.q = t$V2[ipos + 1] - t$V3[ipos] - 1
    d.b = t$V4[ipos + 1] - t$V5[ipos] - 1
    
    if(!((min(d.q, d.b) == 0) && (max(d.q, d.b) <= threshold))){
      ipos = ipos + 1
      next
    }
    # print(ipos)
    #additional sequencies
    n.add = max(d.q, d.b)
    if(d.q != 0){
      s.b = rep('-', n.add)
      s.q = query.fas.chr[(t$V3[ipos]+1):(t$V2[ipos + 1]-1)]
    } else {
      s.q = rep('-', n.add)
      if(t$dir[ipos] == 0){
        s.b = base.fas.fw[(t$V5[ipos]+1):(t$V4[ipos + 1]-1)]  
      } else {
        s.b = base.fas.bw[(t$V5[ipos]+1):(t$V4[ipos + 1]-1)]  
      }
    }
    
    idx = ipos + 1
    t[ipos, 'V3'] <- t[idx, 'V3']
    t[ipos, 'V5'] <- t[idx, 'V5']
    t[ipos, 'V7'] <- t[ipos, 'V7'] + t[idx, 'V7'] + n.add
    t[ipos, 'V8'] <- paste(t[ipos, 'V8'], paste0(s.q, collapse = ''), t[idx, 'V8'], sep = '')
    t[ipos, 'V9'] <- paste(t[ipos, 'V9'], paste0(s.b, collapse = ''), t[idx, 'V9'], sep = '')
    t <- t[-idx,]
  }
  
  rownames(t) <- NULL
  return(t)
}


getCoverageBase <- function(t, len){
  idx = rep(F, len)
  for(irow in 1:nrow(t)){
    idx[t[irow,'V4']:t[irow,'V5']] = T
  }
  coverage = sum(idx) /  length(idx)
  print(coverage)
}

getCoverageQuery <- function(t, len){
  idx = rep(F, len)
  for(irow in 1:nrow(t)){
    idx[t[irow,'V2']:t[irow,'V3']] = T
  }
  coverage = sum(idx) /  length(idx)
  print(coverage)
}


removeShortOverlaps <- function(t, echo=F, len.min = 200, gap.thresh = 10000,
                                filter.q = T, filter.b = T){
  t <- t[order(t[,'V2']),]

  t.base = getBase(t, base.len)

  irow = 1
  while(irow < nrow(t)) {
    # if(irow == 31) break
    # if(t$dir[irow] != t$dir[irow+1]) {
    #   irow = irow + 1
    #   next
    # }
    # 
    # gap.q = t[irow + 1,'V2'] - t[irow,'V3'] - 1
    # gap.b = t[irow + 1,'V4'] - t[irow,'V5'] - 1
    # if(echo) print(c(irow, 
    #                  t.base[irow + 1,'V2'] - t.base[irow,'V3'] - 1, 
    #                  t.base[irow + 1,'V4'] - t.base[irow,'V5'] - 1, 
    #                  t$dir[irow], t$dir[irow+1]))
    
    q.beg = min(t[irow,'V2'], t[irow+1,'V2'])
    q.end = max(t[irow,'V3'], t[irow+1,'V3'])
    pos.q.tmp = rep(0, q.end - q.beg + 1)
    pos.q.tmp[c(t[irow,'V2'] : t[irow,'V3'])- q.beg + 1] = pos.q.tmp[c(t[irow,'V2'] : t[irow,'V3'])- q.beg + 1] + 1
    pos.q.tmp[c(t[irow+1,'V2'] : t[irow+1,'V3'])- q.beg + 1] = pos.q.tmp[c(t[irow+1,'V2'] : t[irow+1,'V3'])- q.beg + 1] + 1
    
    b.beg = min(t.base[irow,'V4'], t.base[irow+1,'V4'])
    b.end = max(t.base[irow,'V5'], t.base[irow+1,'V5'])
    pos.b.tmp = rep(0, b.end - b.beg + 1)
    pos.b.tmp[c(t.base[irow,'V4'] : t.base[irow,'V5'])- b.beg + 1] = pos.b.tmp[c(t.base[irow,'V4'] : t.base[irow,'V5'])- b.beg + 1] + 1
    pos.b.tmp[c(t.base[irow+1,'V4'] : t.base[irow+1,'V5'])- b.beg + 1] = pos.b.tmp[c(t.base[irow+1,'V4'] : t.base[irow+1,'V5'])- b.beg + 1] + 1
    
    
    gap.q = sum(pos.q.tmp == 0)
    gap.b = sum(pos.b.tmp == 0)
    if(gap.q == 0){
      gap.q = -sum(pos.q.tmp == 2)
    }
    if(gap.b == 0){
      gap.b = -sum(pos.b.tmp == 2)
    }
    if(echo) print(c(irow, 
                     gap.q,
                     gap.b,
                     t$dir[irow], t$dir[irow+1]))
    
    
    # if((gap.q < 0) && (gap.b < 0) && (gap.b < gap.q)) gap.q = 10
    

    if(filter.q & (gap.q < 0) & (gap.q > -gap.thresh) & 
       (t[irow, 'V2'] < t[irow+1, 'V2']) & (t[irow, 'V3'] < t[irow+1, 'V3'])){
      # overlap from one side
      s.q1 = strsplit(t[irow, 'V8'], '')[[1]]
      s.b1 = strsplit(t[irow, 'V9'], '')[[1]]
      s.i1 = rep(0, t[irow, 'V7'])
      s.i1[s.q1 != '-'] <- t[irow,'V2']:t[irow,'V3']
      s.q1 = s.q1[min(which(s.i1 >= (t[irow,'V3'] + gap.q + 1))) :t[irow, 'V7']]
      s.b1 = s.b1[min(which(s.i1 >= (t[irow,'V3'] + gap.q + 1))) :t[irow, 'V7']]
      
      
      s.q2 = strsplit(t[irow+1, 'V8'], '')[[1]]
      s.b2 = strsplit(t[irow+1, 'V9'], '')[[1]]
      s.i2 = rep(Inf, t[irow+1, 'V7'])
      s.i2[s.q2 != '-'] <- t[irow+1,'V2']:t[irow+1,'V3']
      
      s.q2 = s.q2[1:max(which(s.i2 <= (t[irow+1,'V2'] - gap.q - 1))) ]
      s.b2 = s.b2[1:max(which(s.i2 <= (t[irow+1,'V2'] - gap.q - 1)))]
      
      
      if(sum(c(s.q2, s.b2) == '-') > sum(c(s.q1, s.b1) == '-') ) {
        tmp = 1
      } else if(sum(c(s.q2, s.b2) == '-') < sum(c(s.q1, s.b1) == '-') ) {
        tmp = 2
      } else if (sum(s.q1 != s.b1) > sum(s.q2 != s.b2)){
        tmp = 2 
      } else {
        tmp = 1
      }
    } else if (filter.b & (gap.b < 0) & (gap.b > -gap.thresh) & 
               (t[irow, 'V4'] < t[irow+1, 'V4']) & (t[irow, 'V5'] < t[irow+1, 'V5'])) {
      
      # stop('a')
      # overlap from one side
      s.q1 = strsplit(t[irow, 'V8'], '')[[1]]
      s.b1 = strsplit(t[irow, 'V9'], '')[[1]]
      s.i1 = rep(0, t[irow, 'V7'])
      if(length(s.i1[s.b1 != '-']) != length(t[irow,'V4']:t[irow,'V5'])){
        message(irow)
        stop('aa')
      }
      s.i1[s.b1 != '-'] <- t[irow,'V4']:t[irow,'V5']
      s.q1 = s.q1[min(which(s.i1 >= (t[irow,'V5'] + gap.b + 1))) :t[irow, 'V7']]
      s.b1 = s.b1[min(which(s.i1 >= (t[irow,'V5'] + gap.b + 1))) :t[irow, 'V7']]
      
      
      s.q2 = strsplit(t[irow+1, 'V8'], '')[[1]]
      s.b2 = strsplit(t[irow+1, 'V9'], '')[[1]]
      s.i2 = rep(Inf, t[irow+1, 'V7'])
      s.i2[s.b2 != '-'] <- t[irow+1,'V4']:t[irow+1,'V5']
      s.q2 = s.q2[1:max(which(s.i2 <= (t[irow+1,'V4'] - gap.b - 1))) ]
      s.b2 = s.b2[1:max(which(s.i2 <= (t[irow+1,'V4'] - gap.b - 1)))]
      
      if(sum(c(s.q2, s.b2) == '-') > sum(c(s.q1, s.b1) == '-') ) {
        tmp = 1
      } else if(sum(c(s.q2, s.b2) == '-') < sum(c(s.q1, s.b1) == '-') ) {
        tmp = 2
      } else if (sum(s.q1 != s.b1) > sum(s.q2 != s.b2)){
        tmp = 2 
      } else {
        tmp = 1
      }
    } else {
      tmp = 0
    }
    if(tmp == 1){  # remain a part at the first piece
      # remove a part from the second piece
      t[irow+1, 'V2'] <- t[irow+1, 'V2'] + sum(s.q2 != '-')
      t[irow+1, 'V4'] <- t[irow+1, 'V4'] + sum(s.b2 != '-')
      
      t[irow+1, 'V8'] <- removeFirstN(t[irow+1, 'V8'], length(s.q2))
      t[irow+1, 'V9'] <- removeFirstN(t[irow+1, 'V9'], length(s.b2))
      
      t[irow+1, 'V7'] <- nchar(t[irow+1, 'V8'])
      
    } else if(tmp == 2){ # remain a part at the second piece
      # remove a part from the first piece
      t[irow, 'V3'] <- t[irow, 'V3'] - sum(s.q1 != '-')
      t[irow, 'V5'] <- t[irow, 'V5'] - sum(s.b1 != '-')
      
      t[irow, 'V8'] <- removeLastN(t[irow, 'V8'], length(s.q1))
      t[irow, 'V9'] <- removeLastN(t[irow, 'V9'], length(s.b1))
      t[irow, 'V7'] <- nchar(t[irow, 'V8'])
    }
    irow = irow + 1
  }
  t = t[t[,'V7'] > len.min,]
  return(t)
}

returnScoreFirstN <- function(t, jrow, n.overlap, by.base, echo=F){
  s.q = strsplit(t$V8[jrow], '')[[1]]
  s.b = strsplit(t$V9[jrow], '')[[1]]
  lens = c(sum(s.q != '-'), sum(s.b != '-'))
  pos = rep(0, t$V7[jrow])
  if(by.base){
    pos[s.b != '-'] = 1  
  } else {
    pos[s.q != '-'] = 1
  }
  pos = cumsum(pos)
  idx = max(which(pos == n.overlap))
  s.q = s.q[1:idx]
  s.b = s.b[1:idx]
  score = sum(s.q != s.b)
  if(echo){
    print(s.q)
    print(s.b)
  }
  return(c(score, lens))
}


returnScoreLastN <- function(t, jrow, n.overlap, by.base, echo=F){
  s.q = strsplit(t$V8[jrow], '')[[1]]
  s.b = strsplit(t$V9[jrow], '')[[1]]
  lens = c(sum(s.q != '-'), sum(s.b != '-'))
  pos = rep(0, t$V7[jrow])
  if(by.base){
    pos[s.b != '-'] = 1  
  } else {
    pos[s.q != '-'] = 1
  }
  pos = rev(cumsum(rev(pos)))
  idx = min(which(pos == n.overlap ))
  n.len = length(s.q)
  s.q = s.q[idx:n.len]
  s.b = s.b[idx:n.len]
  score = sum(s.q != s.b)
  if(echo){
    print(s.q)
    print(s.b)
  }
  return(c(score, lens))
}

removeAlnLastN <- function(t, jrow, n.overlap, by.base){
  s.q = strsplit(t$V8[jrow], '')[[1]]
  s.b = strsplit(t$V9[jrow], '')[[1]]
  lens = c(sum(s.q != '-'), sum(s.b != '-'))
  pos = rep(0, t$V7[jrow])
  if(by.base){
    pos[s.b != '-'] = 1  
  } else {
    pos[s.q != '-'] = 1
  }
  pos = rev(cumsum(rev(pos)))
  idx = min(which(pos == n.overlap ))
  n.len = length(s.q)
  s.q = s.q[-(idx:n.len)]
  s.b = s.b[-(idx:n.len)]
  t$V8[jrow] = paste0(s.q, sep = '', collapse = '')
  t$V9[jrow] = paste0(s.b, sep = '', collapse = '')
  t$V3[jrow] = t$V2[jrow] + sum(s.q != '-') - 1
  t$V5[jrow] = t$V4[jrow] + sum(s.b != '-') - 1
  t$V7[jrow] = length(s.q)
  return(t)
}

removeAlnFirstN <- function(t, jrow, n.overlap, by.base){
  s.q = strsplit(t$V8[jrow], '')[[1]]
  s.b = strsplit(t$V9[jrow], '')[[1]]
  lens = c(sum(s.q != '-'), sum(s.b != '-'))
  pos = rep(0, t$V7[jrow])
  if(by.base){
    pos[s.b != '-'] = 1  
  } else {
    pos[s.q != '-'] = 1
  }
  pos = cumsum(pos)
  idx = max(which(pos == n.overlap))
  s.q = s.q[-(1:idx)]
  s.b = s.b[-(1:idx)]
  t$V8[jrow] = paste0(s.q, sep = '', collapse = '')
  t$V9[jrow] = paste0(s.b, sep = '', collapse = '')
  t$V2[jrow] = t$V3[jrow] - sum(s.q != '-') + 1
  t$V4[jrow] = t$V5[jrow] - sum(s.b != '-') + 1
  t$V7[jrow] = length(s.q)
  return(t)
}

removeShortOverlaps2 <- function(t, base.len, echo=F, len.min = 200, gap.thresh = 10000,
                                filter.q = T, filter.b = T){
  
  # Overlap in query
  t <- t[order(t[,'V2']),]
  irow = 1
  while(irow < nrow(t)) {
    
    if(t$V2[irow+1] > t$V3[irow]){
      irow = irow + 1
      next
    }
    # stop(irow)
    n.overlap = t$V3[irow] - t$V2[irow+1] + 1
    if(n.overlap > (t$V3[irow] - t$V2[irow] + 1)){
      t = t[-irow,]
      next
    }
    
    if(n.overlap > (t$V3[irow + 1] - t$V2[irow + 1] + 1)){
      t = t[-irow,]
      next
    }
    
    # get score of alignment form the irow
    
    score0 = returnScoreLastN(t, irow, n.overlap, by.base = F, echo = echo)
    score1 = returnScoreFirstN(t, irow+1, n.overlap, by.base = F, echo = echo)
    
    
    # Cut the sequence
    idx.cut = -1
    if(score0[1] < score1[1]){
      idx.cut = 1
    } else if (score0[1] > score1[1]){
      idx.cut = 0
    } else {
      if(score0[2] > score1[2]){
        idx.cut = 1
      } else if (score0[2] < score1[2]){
        idx.cut = 0
      } else {
        if(score0[3] > score1[3]){
          idx.cut = 1
        } else if (score0[3] < score1[3]){
          idx.cut = 0
        } 
      }
    }
    # print(idx.cut)
    
    # break
    
    
    if(idx.cut == -1){ 
      # identical, remove both
      t = t[-c(irow, irow + 1),]
    } else if(idx.cut == 0){
      # remove from the irow
      t = removeAlnLastN(t, irow, n.overlap, by.base = F)
    } else { 
      # remove from the irow + 1
      t = removeAlnFirstN(t, irow+1, n.overlap, by.base = F)
    }
    
    # order again and start from the beginning
    t = t[t[,'V7'] > len.min,]
    # t = removeCompleteOverlaps2(t, base.len)
    t <- t[order(t$V2),]
    rownames(t) = NULL
    irow = 1 # why? just because
  }
  
  # ---------------------------------------------------------------
  
  
  # overlap in base
  t = t[t[,'V7'] > len.min,]
  t.base = getBase(t, base.len)
  t <- t[order(t.base$V4),]
  rownames(t) = NULL
  t.base = getBase(t, base.len)
  irow = 1
  while(irow < nrow(t)) {

    
    if(t.base$V4[irow+1] > t.base$V5[irow]){
      irow = irow + 1
      next
    }
    # print(irow)
    # if(231 == irow) stop()
    
    n.overlap = t.base$V5[irow] - t.base$V4[irow+1] + 1

    if(n.overlap > (t.base$V5[irow] - t.base$V4[irow] + 1)){
      t = t[-irow,]
      rownames(t) = NULL
      t.base = getBase(t, base.len)
      next
    }
    if(n.overlap > (t.base$V5[irow + 1] - t.base$V4[irow + 1] + 1)){
      t = t[-(irow+1),]
      rownames(t) = NULL
      t.base = getBase(t, base.len)
      next
    }
    
    
    # get score of alignment form the irow
    if (t.base$dir[irow] == 0) { # get alignment form the end part
      score0 = returnScoreLastN(t, irow, n.overlap, by.base = T, echo = echo)
    } else { # get alignment from the front part
      score0 = returnScoreFirstN(t, irow, n.overlap, by.base = T, echo = echo)
    } 
    
    # get score of alignment form the irow + 1
    if (t.base$dir[irow+1] == 0) { # get alignment form the end part
      score1 = returnScoreFirstN(t, irow+1, n.overlap, by.base = T, echo = echo)
    } else { # get alignment from the front part
      score1 = returnScoreLastN(t, irow+1, n.overlap, by.base = T, echo = echo)
    } 
    # print(score0)
    # print(score1)
    
    
    # Cut the sequence
    idx.cut = -1
    if(score0[1] < score1[1]){
      idx.cut = 1
    } else if (score0[1] > score1[1]){
      idx.cut = 0
    } else {
      if(score0[2] > score1[2]){
        idx.cut = 1
      } else if (score0[2] < score1[2]){
        idx.cut = 0
      } else {
        if(score0[3] > score1[3]){
          idx.cut = 1
        } else if (score0[3] < score1[3]){
          idx.cut = 0
        } 
      }
    }
    # print(idx.cut)

    # break
    
    
    if(idx.cut == -1){ 
      # identical, remove both
      t = t[-c(irow, irow + 1),]
    } else if(idx.cut == 0){
      # remove from the irow
      if(t$dir[irow] == 0){
        t = removeAlnLastN(t, irow, n.overlap, by.base = T)
      } else {
        t = removeAlnFirstN(t, irow, n.overlap, by.base = T)
      }
    } else { 
      # remove from the irow + 1
      if(t$dir[irow+1] == 0){
        t = removeAlnFirstN(t, irow+1, n.overlap, by.base = T)
      } else {
        t = removeAlnLastN(t, irow+1, n.overlap, by.base = T)
      }
    }
    
    # order again and start from the beginning
    t = t[t[,'V7'] > len.min,]
    t.base = getBase(t, base.len)
    t <- t[order(t.base$V4),]
    rownames(t) = NULL
    t.base = getBase(t, base.len)
    irow = 1 # why? just because
  }
  # ---------------------------------------------------------------
    

  t = t[t[,'V7'] > len.min,]
  return(t)
}

additionalLocalAlignments <- function(t, query.fas.chr, base.fas.fw, base.fas.bw, echo=F, 
                                      n.short=200, gap.max = 15000, gap.open = 30, gap.ext = 0.5, 
                                      file.log = NULL){
  if(!is.null(file.log)) {
    write('', file=file.log, append=F)
  }
  
  t <- t[order(t[,'V2']),]
  t.additional <- c()
  irow = 0
  while(irow < (nrow(t) - 2)) {
    irow = irow+1
    if (t[irow,'dir'] != t[irow+1,'dir']) next
    
    gap.q = t[irow + 1,'V2'] - t[irow,'V3'] - 1
    gap.b = t[irow + 1,'V4'] - t[irow,'V5'] - 1
    
    
    
    if ((gap.q <= 0) || (gap.b <= 0)) next
    if ((gap.q/gap.b > 10) ||(gap.b/gap.q > 10)) next
    if (max(gap.q, gap.b) > gap.max) next
    
    
    s.q = query.fas.chr[(t[irow,'V3']+1): (t[irow + 1,'V2']-1)]
    if(t$dir[irow] == 0) {
      s.b = base.fas.fw[(t[irow,'V5']+1): (t[irow + 1,'V4']-1)]
    } else {
      s.b = base.fas.bw[(t[irow,'V5']+1): (t[irow + 1,'V4']-1)]
    }
    
    s.q <- paste0(s.q, collapse='')
    s.b <- paste0(s.b, collapse='')
    
    localAlign <- pairwiseAlignment(s.q, s.b, 
                                    gapOpening = gap.open, gapExtension=gap.ext,
                                    type="local")
    
    if(!is.null(file.log)) {
      write(paste((t[irow,'V3']+1), (t[irow,'V5']+1), as.character(alignedPattern(localAlign))),
            file=file.log, append=TRUE)
      write(paste((t[irow,'V3']+1), (t[irow,'V5']+1), as.character(alignedSubject(localAlign))), 
            file=file.log, append=TRUE)
      write('\n', file=file.log, append=TRUE)
    }
    
    # localAlign <- pairwiseAlignment(s.q, s.b, type="global")
    
    s.q.range <- localAlign@pattern@range
    s.b.range <- localAlign@subject@range
    
    if(nchar(localAlign@pattern) < n.short) next
    
    if(echo) print(c(irow, gap.q, gap.b, nchar(localAlign@pattern)))
    
    t.tmp = data.frame(V1=paste('add', t[irow,'V1'], sep = '_'),
                       V2=t[irow,'V3'] + s.q.range@start,  # You do not need +1 or -1 here
                       V3=t[irow,'V3'] + s.q.range@start + s.q.range@width - 1,
                       V4=t[irow,'V5'] + s.b.range@start,
                       V5=t[irow,'V5'] + s.b.range@start + s.b.range@width - 1,
                       V6=100, V7=nchar(localAlign@pattern),
                       V8=as.character(alignedPattern(localAlign)), 
                       V9=as.character(alignedSubject(localAlign)), 
                       V10=t[irow,10], dir=t[irow,'dir'])
    t.additional <- rbind(t.additional, t.tmp)
  }
  
  t.additional[, 'V9'] <- as.character(t.additional[, 'V9'])
  
  t <- rbind(t, t.additional)
  
  checkCorrespToGenome(t, query.fas = query.fas.chr, 
                       base.fas.fw = base.fas.fw, 
                       base.fas.bw = base.fas.bw)
  t <- t[order(t[,'V2']),]
  
  return(t)
}

removeSomeOverlaps <- function(t, base.len, filter.query = T, filter.base = T, sort.len = T){
  rownames(t) <- NULL
  
  if(sort.len){
    t = t[order(-t$V7),]
  }
  
  # Get initial positions in base sequence
  t.base = getBase(t, base.len)
  
  
  idx.remove = c()
  # ------------------------------
  # Overlap in query
  if(filter.query){
    for(irow in 1:nrow(t)){
      if(irow %in% idx.remove) next
      idx = which((t.base[,'V2'] <= t.base[irow,'V3']) & 
                    (t.base[,'V3'] >= t.base[irow,'V2']))
      idx = setdiff(idx, irow)
      idx.remove = c(idx.remove, idx)
    }    
  }
  
  
  # Overlap in base
  if(filter.base) {
    for(irow in 1:nrow(t)){
      if(irow %in% idx.remove) next
      idx = which((t.base[,'V4'] <= t.base[irow,'V5']) & 
                    (t.base[,'V5'] >= t.base[irow,'V4']))
      idx = setdiff(idx, irow)
      idx.remove = c(idx.remove, idx)
    }    
  }
  
  if(length(idx.remove) == 0){
    return(t)
  }
  idx.remove = unique(idx.remove)
  t = t[-idx.remove,]
  t <- t[order(t[,'V2']),]
  rownames(t) <- NULL
  return(t)
}

removeCompleteOverlaps2 <- function(t, base.len, filter.query = T, filter.base = T, sort.len = F){
  rownames(t) <- NULL
  
  # Get initial positions in base sequence
  t.base = getBase(t, base.len)
  if(sort.len){
    t.base = t.base[order(t.base$V7),]
  }
  
  idx.remove = c()
  # ------------------------------
  # Overlap in query
  if(filter.query){
    for(irow in 1:nrow(t)){
      idx = which((t.base[,'V2'] >= t.base[irow,'V2']) & 
                    (t.base[,'V3'] <= t.base[irow,'V3']) )
      idx = setdiff(idx, irow)
      idx.remove = c(idx.remove, idx)
    }    
  }
  
  
  # Overlap in base
  if(filter.base) {
    for(irow in 1:nrow(t)){

      idx = which((t.base[,'V4'] >= t.base[irow,'V4']) & 
                    (t.base[,'V5'] <= t.base[irow,'V5']) )
      idx = setdiff(idx, irow)
      idx.remove = c(idx.remove, idx)
    }    
  }
  
  if(length(idx.remove) == 0){
    return(t)
  }
  
  t = t[-idx.remove,]
  t <- t[order(t[,'V2']),]
  rownames(t) <- NULL
  return(t)
}

removeCompleteOverlaps <- function(t, n.distant = 1000000, n.short = 200, 
                                   filter.query = T, filter.base = T,
                                   centromere.pos = NULL){
  rownames(t) <- NULL
  
  # Get initial positions in base sequence
  t.base = t
  tmp  = base.len - t.base[t.base$dir == 1,'V5'] + 1
  t.base[t.base$dir == 1,'V5'] = base.len - t.base[t.base$dir == 1,'V4'] + 1
  t.base[t.base$dir == 1,'V4'] = tmp
  
  t.base = t.base[rev(order(t.base['V7'])),]
  
  # ------------------------------

  
  # Overlap in query
  if(filter.query){
    irow = 1
    while(irow < nrow(t.base)) {
      # print(nrow(t.base))
      idx = which((t.base[,'V2'] >= t.base[irow,'V2']) & 
                    (t.base[,'V3'] <= t.base[irow,'V3']) )#& (t.base$dir == t.base$dir[irow]))
      idx = setdiff(idx, irow)
      idx = idx[t.base[irow,'V7'] > t.base[idx,'V7']]
      if(length(idx) >= 1) {
        t.base = t.base[-idx,] 
      } else {
        irow = irow + 1
      }
    }    
  }

  
  # Overlap in base
  if(filter.base) {
    irow = 1
    while(irow < nrow(t.base)) {
      # print(nrow(t.base))
      idx = which((t.base[,'V4'] >= t.base[irow,'V4']) & 
                    (t.base[,'V5'] <= t.base[irow,'V5']) )# & (t.base$dir == t.base$dir[irow]))
      idx = setdiff(idx, irow)
      idx = idx[t.base[irow,'V7'] > t.base[idx,'V7']]
      if(length(idx) >= 1) {
        t.base = t.base[-idx,] 
      } else {
        irow = irow + 1
      }
    }    
  }

  
  # idx = rep(F, base.len)
  # for(irow in 1:nrow(t.base)){
  #   idx[t.base[irow,'V4']:t.base[irow,'V5']] = T
  # }
  # 
  # coverage = sum(idx) / base.len
  # print(coverage)
  
  # idx = rep(F, length(query.fas.chr))
  # for(irow in 1:nrow(t.base)){
  #   idx[t.base[irow,'V2']:t.base[irow,'V3']] = T
  # }
  # 
  # coverage = sum(idx) / length(idx)
  # print(coverage)
  
  # if(!is.null(n.short)){
  #   t.base <- t[t.base[,'V7'] > n.short,]    
  # }
  
  # # remove distant outliers
  # if(!is.null(n.distant)){
  #   if(is.null(centromere.pos)){
  #     idx.distant <- c()
  #     for(irow in 1:nrow(t.base)){
  #       if(abs(t.base[irow,'V2'] - t.base[irow,'V4']) > n.distant){
  #         idx.distant <- c(idx.distant, irow)
  #       }
  #     }
  #     if(length(idx.distant) > 1)  t.base <- t.base[-idx.distant,]  
  #   } else {
  # 
  #     t.base = t.base[order(t.base[,'V2']),]
  #     idx = min(which(!((t.base[,'V2'] <= centromere.pos[1]) & (t.base[,'V3'] <= centromere.pos[2])) )) - 1
  #     while(t.base[idx,'V7'] > 10000) idx = idx + 1
  #     cent.start.idx = idx
  #     
  #     idx = max(which(!((t.base[,'V2'] >= centromere.pos[1]) & 
  #                         (t.base[,'V3'] >= centromere.pos[2])) )) + 1
  #     while(t.base[idx,'V7'] > 10000) idx = idx - 1
  #     cent.end.idx = idx
  #     
  #     t.arms = list( t.base[1:(cent.start.idx-1),], 
  #                    t.base[(cent.end.idx+1):nrow(t.base),])
  #     
  #     for(i.tmp in 1:2){
  #       #fit linear model
  #       # lm.res <- summary(lm('V2 ~ offset(1*V4)', t.arms[[i.tmp]][(t.arms[[i.tmp]][,'V7'] > 10^4)&
  #       #                                                     (t.arms[[i.tmp]]$dir == 0),]))
  #       # intercept <- lm.res$coefficients[1,1] + 1.5 * lm.res$coefficients[1,2]
  #       
  #       y = t.arms[[i.tmp]][(t.arms[[i.tmp]][,'V7'] > 10^4)&
  #                             +                         (t.arms[[i.tmp]]$dir == 0),]
  #       intersept <- max(abs(y[,'V2'] - y[,'V4']))
  #       
  #       idx.distant <- c()
  #       for(irow in 1:nrow(t.arms[[i.tmp]])){
  #         if(abs(t.arms[[i.tmp]][irow,'V2'] - t.arms[[i.tmp]][irow,'V4']) > intercept + n.distant){
  #           idx.distant <- c(idx.distant, irow)
  #         }
  #       }
  #       if(length(idx.distant) > 1)  t.arms[[i.tmp]] <- t.arms[[i.tmp]][-idx.distant,]  
  #     }
  #     t.base = t.base[sort(c(rownames(t.arms[[1]]),
  #                       rownames(t.arms[[2]]),
  #                       rownames(t.base[cent.start.idx:cent.end.idx,]))),]
  #     
  #   }
  #   
  # }
  
  t = t[rownames(t.base),]
  
  # plot(t.base[,'V4'], t.base[,'V2'])


  
  # plot(t[,'V4'], t[,'V2'])
  
  t <- t[order(t[,'V2']),]
  rownames(t) <- NULL
  
  return(t)
}


setDir <- function(t, base.len){
  # direction
  t$dir = c()
  idx.dir = t[,'V5'] < t[,'V4']
  pos1 = base.len - t[idx.dir,'V5'] + 1
  pos2 = base.len - t[idx.dir,'V4'] + 1
  
  t[idx.dir,'V5'] <- base.len - t[idx.dir,'V5'] + 1
  t[idx.dir,'V4'] <- base.len - t[idx.dir,'V4'] + 1
  t[, 'dir'] = 0
  t[idx.dir, 'dir'] = 1
  
  return(t)
}



getT <- function(t.file, query.fas.chr, base.fas.fw, base.fas.bw,
                 thresholds = c(100, 500, 1000, 1500, 2000), echo=T){
  # ------- Read blast results -------
  
  if(echo) message(paste0(c('Reading', t.file, '...'), collapse = ' '))
  base.len = length(base.fas.fw)
  t = read.table(t.file, stringsAsFactors = F, header = F)
  
  t <- setDir(t, base.len)
  
  # Get right positions
  start.pos = as.numeric(sapply(strsplit(t[,1], "\\|"), "[", 4)) - 1
  t[,2:3] = t[,2:3] + start.pos
  
  rownames(t) <- NULL
  
  # check
  for(irow in 1:nrow(t)) {
    if(nchar(t[irow,'V8']) != nchar(t[irow,'V9'])) stop('aaa')
  }
  
  source("/Users/anna/OneDrive/pushkin/cryptic/nanopore/scripts/synteny_infer.R")
  checkCorrespToGenome(t, query.fas = query.fas.chr, 
                       base.fas.fw = base.fas.fw, 
                       base.fas.bw = base.fas.bw)
  
  
  # ------- Grow up the alignment -------
  
  ## ---- First glue ----
  
  if(echo) message('First glue (exact)...')
  t <- glueZero(t)
  
  ## ---- Pairwise global alignments of gaps ----
  
  if(echo) message('Glue with thresholds')
  t = glueByThreshold(t, thresholds, query.fas = query.fas.chr, 
                      base.fas.fw = base.fas.fw, 
                      base.fas.bw = base.fas.bw, file.log=file.log)
  
  ## ---- Remove complete overlap from both sides (base an query) ----
  
  if(echo) message('Remove complete overlaps')
  t <- removeCompleteOverlaps(t)
  
  ## ---- Remove short overlaps and Local alignments ----
  
  # Remove short overlaps
  if(echo) message('Remove short overlaps')
  t <- removeShortOverlaps(t, echo = F)
  
  # Additional alignments
  if(echo) message('Additional alignments')
  nrow.tmp = nrow(t) - 1
  while(nrow(t) != nrow.tmp) {
    nrow.tmp = nrow(t)
    # print(nrow.tmp)
    t <- additionalLocalAlignments(t, query.fas.chr, base.fas.fw, base.fas.bw, echo = F, n.short=100,
                                   file.log = file.log)
    t <- glueZero(t)
    
    checkCorrespToGenome(t, query.fas = query.fas.chr, 
                         base.fas.fw = base.fas.fw, 
                         base.fas.bw = base.fas.bw)
    
    t = glueByThreshold(t, thresholds, query.fas = query.fas.chr, 
                        base.fas.fw = base.fas.fw, 
                        base.fas.bw = base.fas.bw, file.log=file.log)
  }
  
  return(t)
}


plotSyntenyBlocks <-function(t, base.len, idx=NULL, hlines=NULL, vlines=NULL){
  df = c()
  for(i in 1:nrow(t)) {
    if(t$dir[i] == 0) {
      df = rbind(df, c(t[i, 'V2'], t[i, 'V4'], i, 0))
      df = rbind(df, c(t[i, 'V3'], t[i, 'V5'], i, 0))
    } else {
      df = rbind(df, c(t[i, 'V2'], base.len - t[i, 'V4'] + 1, i, 1))
      df = rbind(df, c(t[i, 'V3'], base.len - t[i, 'V5'] + 1, i, 1))
    }
  }
  
  df = as.data.frame(df)
  df$clr <- df$V3 %% 2;
  
  if(!is.null(idx)) df$V4[c(idx*2, (idx*2-1))] = 2
  
  p <- ggplot(df, aes(x = V1, y=V2, color = as.factor(V4), group=as.factor(V3)  )) + 
    geom_line(show.legend = FALSE) + theme_bw() + xlab('query') + ylab('base')
  
  if(!is.null(hlines)){
    p <- p + geom_hline(yintercept=hlines, linetype="dashed", color='grey50')
  }
  
  if(!is.null(vlines)){
    p <- p + geom_vline(xintercept=vlines, linetype="dashed", color='grey50')
  }
  
  return(p)
}


orderT <- function(t, order = 'query'){
  if(order == 'query'){
    t = t[order(t[,'V2']),]
  } else if(order == 'base'){
    t = t[order(t[,'V4']),]
  } else {
    message('Not ordered!')
  }
  rownames(t) <- NULL
  return(t)
}

getBlastStat <- function(t.reb, cover.cutoff = 0.9, 
                         flag.max.cover = T,
                         flag.tot.cover = T,
                         flag.pers.cover = T,
                         flag.max.cover.one=F,
                         suff.res = '') {
  
  reb.unique = unique(t.reb[, 'V1'])
  t.reb = t.reb[!duplicated(t.reb[,1:10]),]

  reb.res = c()
  names.res = c()
  
  for(i.reb in 1:length(reb.unique)){
    reb.name = reb.unique[i.reb]
    
    tmp.res = c()
    names.res = c()
    
    len.tmp = as.numeric(strsplit(reb.name, '\\|')[[1]][4])
    n.names.init = length(names.res)
    
    # if record is here
    
    t.tmp = t.reb[t.reb[,'V1'] == reb.name,, drop=F]
    
    
    if(flag.max.cover){
      idx = which(t.tmp[,'V7'] >= cover.cutoff * len.tmp)
      tmp = length(idx)
      if(tmp < 1) {
        tmp = max(t.tmp[,'V7']) / len.tmp
      }
      tmp.res = c(tmp.res, tmp)
      names.res = c(names.res, 'max_coverage')
    }
    
    if(flag.max.cover.one){
      

      tmp = max(t.tmp[,'V3'] - t.tmp[,'V2'] + 1) / len.tmp
      tmp.res = c(tmp.res, tmp)
      names.res = c(names.res, 'max_coverage_one')
    }
    
    if(flag.tot.cover){
      t.coverage = rep(0, len.tmp)
      for(i.tmp in 1:nrow(t.tmp)) {
        t.coverage[t.tmp[i.tmp, 'V2']: t.tmp[i.tmp, 'V3']] = 1
      }
      t.coverage = sum(t.coverage) / len.tmp 
      tmp.res = c(tmp.res, t.coverage)
      names.res = c(names.res, 'tot_coverage')
    }
    
    if(flag.pers.cover){
      t7 = t.tmp[,'V7'] / len.tmp
      t7[t7 > 1] = 1
      tmp = table(c(floor(t7*10), 0:10)) - 1  
      tmp.res = c(tmp.res, rev(tmp))
      names.res = c(names.res, paste('pers', rev(0:10)/10, sep = '_'))
    }
    
    
    if(nchar(suff.res) > 0){
      names.res[(n.names.init+1):length(names.res)] <- 
        paste(names.res[(n.names.init+1):length(names.res)], suff.res, sep = '_')
    }
    
    names(tmp.res) <- names.res
    reb.res = rbind(reb.res, tmp.res)
  }
  rownames(reb.res) <- reb.unique
  
  return(reb.res)
}


# ignore sequences of this size
# thresh - threshold fo remove reduldant sequences
getMobilome <- function(t, base.len, allowed_gap = 5, m.min.length = 15, 
                        thresh = 0.8, m.max.length = Inf, echo=F){
  

  m1 = c()
  m2 = c()
  pos.m1 = c()
  pos.m2 = c()
  for(irow in 1:nrow(t)){
    if(echo) print(irow)
    s1 = strsplit(t[irow,'V8'], '')[[1]]
    s2 = strsplit(t[irow,'V9'], '')[[1]]
    
    pos = rep(0, length(s1))
    pos1 = pos; pos1[s1 != '-'] = t[irow,'V2']: t[irow,'V3']
    pos2 = pos; pos2[s2 != '-'] = t[irow,'V4']: t[irow,'V5']
    
    icol = 0
    
    # icol_gap = which((s1 == '-') |(s2 =='_'))
    # icol_gap
    
    while(icol < length(s1)){
      icol = icol + 1
      
      if(s1[icol] == '-') {  # start gaps on query
        i.start = icol
        s = c()
        while((sum(s1[icol:min((icol+allowed_gap), length(s1))] == '-') > 0) && (icol <= length(s1))){
          s = c(s, s2[icol])
          icol = icol + 1
        }
        i.end = icol - 1
        
        if((i.end - i.start + 1) < m.min.length) next
        if((i.end - i.start + 1) > m.max.length) next
        
        # Nucleotide frequences
        f = table(s)
        if(max(f) > thresh * length(s)) next
        
        # Di-Nucleotide frequences
        cf = apply(combn(f, 2), 2, sum)
        if(max(cf) > thresh * length(s)) next
        
        if(t[irow, 'dir'] == 1){
          pos2.start <- base.len - pos2[i.end] + 1
          pos2.end <- base.len - pos2[i.start] + 1
        } else {
          pos2.start <- pos2[i.start]
          pos2.end <- pos2[i.end]
        }
        
        if(i.start == 1){
          # print(dim(m2))
          if(length(c(pos2.start, pos2.end, 0, length(s), paste0(s, collapse = ''))) != 5) stop('ddd')
          m2 = rbind(m2, c(pos2.start, pos2.end, 0, length(s), paste0(s, collapse = ''))) 
        } else {
          # print(dim(m2))
          if(length(c(pos2.start, pos2.end, pos1[i.start-1], length(s), paste0(s, collapse = ''))) != 5) stop('eee')
          m2 = rbind(m2, c(pos2.start, pos2.end, pos1[i.start-1], length(s), paste0(s, collapse = ''))) 
        }
        
      } else if(s2[icol] == '-') {  # start gaps on query
        i.start = icol
        s = c()
        while((sum(s2[icol:min((icol+allowed_gap), length(s2))] == '-') > 0) && (icol <= length(s2))){
          s = c(s, s1[icol])
          icol = icol + 1
        }
        i.end = icol - 1
        if((i.end - i.start + 1) < m.min.length) next
        if((i.end - i.start + 1) > m.max.length) next
        
        # Nucleotide frequences
        f = table(s)
        if(max(f) > thresh * length(s)) next
        
        # Di-Nucleotide frequences
        cf = apply(combn(f, 2), 2, sum)
        if(max(cf) > thresh * length(s)) next
        
        if(i.start == 1){
          # print(dim(m1))
          if(length(c(pos1[i.start], pos1[i.end], 0, length(s), paste0(s, collapse = ''))) != 5) stop('bb')
          m1 = rbind(m1, c(pos1[i.start], pos1[i.end], 0, length(s), paste0(s, collapse = '')))  
        } else {
          
          if(t[irow, 'dir'] == 1){
            pos2.start <- base.len - pos2[i.start-1] + 1
          } else {
            pos2.start = pos2[i.start-1]
          }
          # print(dim(m1))
          if(length(c(pos1[i.start], pos1[i.end], pos2.start, length(s), paste0(s, collapse = ''))) != 5) stop('aa')
          m1 = rbind(m1, c(pos1[i.start], pos1[i.end], pos2.start, length(s), paste0(s, collapse = '')))
        }
        
      }
    
      
      if(sum(is.na(m1)) >0) stop('NA1')
      if(sum(is.na(m2)) >0) stop('NA2')
    }
  }
  
  colnames(m1) <- c('pos.q.start', 'pos.q.end', 'pos.b.before', 'length', 'seq')
  colnames(m2) <- c('pos.b.start', 'pos.b.end', 'pos.q.before', 'length', 'seq')
  
  m1 = as.data.frame(m1)
  m2 = as.data.frame(m2)
  
  for(i in 1:4){
    m1[,i] =  as.numeric(as.character(m1[,i]))
    m2[,i] =  as.numeric(as.character(m2[,i]))
  }
  
  return(list(query = m1, base = m2))
}

addRowsByNames <- function(reb, s.fasta){
  
  init.names = rownames(reb)
  lost.names = setdiff(names(s.fasta), init.names)
  for(name in lost.names){
    reb = rbind(reb, rep(0, ncol(reb)))
  }
  rownames(reb) <- c(init.names, lost.names)
  reb <- reb[names(s.fasta),, drop=F]
  return(reb)
  
}

showd <- function(t){
  t.base
  links.q = t[2:(nrow(t)),'V2'] - t[1:(nrow(t)-1),'V3'] - 1
  links.b = t[2:(nrow(t)),'V4'] - t[1:(nrow(t)-1),'V5'] - 1
  print(cbind(links.q, links.b, t$dir[2:(nrow(t))], t$dir[1:(nrow(t)-1)]))
}


njColumns <- function(mx){
  n.tmp = ncol(mx)
  d = matrix(0, nrow = n.tmp, ncol = n.tmp, 
             dimnames = list(colnames(mx), colnames(mx)))
  for(i in 1:n.tmp){
    for(j in 1:n.tmp){
      d[i,j] = sum(mx[,i] != mx[,j])
    }
  }
  d1 <- as.dist(d)
  tree <- nj(d1)
  return(tree)
}

hclustColumns <- function(mx){
  n.tmp = ncol(mx)
  d = matrix(0, nrow = n.tmp, ncol = n.tmp, 
             dimnames = list(colnames(mx), colnames(mx)))
  for(i in 1:n.tmp){
    for(j in 1:n.tmp){
      d[i,j] = sum(mx[,i] != mx[,j])
    }
  }
  d1 <- as.dist(d)
  hc <- hclust(d1)
  return(hc)
}

s2blocks <- function(s, block.len, echo=F){
  blocks = c()
  positions = c()
  block.n = floor(length(s) / block.len)
  if(block.n == 0){
    block.n = 1
  }
  for(i in 1:block.n){
    pos.s = 1 + (i-1)*block.len
    if(i == block.n){
      pos.e = length(s)
    } else {
      pos.e = i*block.len
    }
    if(echo) print(c(pos.s, pos.e, pos.e - pos.s + 1))
    positions = rbind(positions, c(pos.s, pos.e))
    blocks = c(blocks, paste0(s[pos.s:pos.e], collapse = ''))
  }
  return(list(b=blocks, p=positions))
}


alignBlocks <- function(t.arm, irow.q, irow.b, 
                        query.fas.chr, base.fas.fw,
                        base.fas.bw, block.len, gap.open=10, gap.ext=0.5,
                        echo = F, g = NULL, rev.compl=F, diag.cutoff = NULL){
  
  t.arm.base = t.arm
  tmp  = base.len - t.arm.base[t.arm.base$dir == 1,'V5'] + 1
  t.arm.base[t.arm.base$dir == 1,'V5'] = base.len - t.arm.base[t.arm.base$dir == 1,'V4'] + 1
  t.arm.base[t.arm.base$dir == 1,'V4'] = tmp
  
  
  if(echo) print(c(irow.q, irow.b))
  s.q = query.fas.chr[(t.arm[irow.q,'V3']+1) : (t.arm[irow.q+1,'V2']-1)]
  if(t.arm.base[irow.b,'V4'] > t.arm.base[irow.b+1,'V5']){
    s.b = base.fas.fw[(t.arm.base[irow.b+1,'V5']+1) : (t.arm.base[irow.b,'V4']-1)]
    pos.start = t.arm.base[irow.b+1,'V5']
  } else {
    s.b = base.fas.fw[(t.arm.base[irow.b,'V5']+1) : (t.arm.base[irow.b+1,'V4']-1)]
    pos.start = t.arm.base[irow.b,'V5']
  }
  # if(t.arm$dir[irow.b] == 0){
    # s.b = base.fas.fw[(t.arm.base[irow.b,'V5']+1) : (t.arm.base[irow.b+1,'V4']-1)]
  # } else {
    # s.b = base.fas.bw[(t.arm[irow.b,'V5']+1) : (t.arm[irow.b+1,'V4']-1)]
  # }
  
  if(rev.compl){
    s.b = reverseComplement(s.b)
  }
  
  # Split into blocks
  blocks.q <- s2blocks(s.q, block.len)
  blocks.b <- s2blocks(s.b, block.len)
  
  # Alignments
  t.tmp = c()
  mx = matrix(0, nrow = length(blocks.q$b), ncol = length(blocks.b$b))
  for(i.q in 1:length(blocks.q$b)){
    for(i.b in 1:length(blocks.b$b)){
      
      
      if(str_count(blocks.q$b[i.q], "n") / nchar(str_count(blocks.q$b[i.q], "n")) > 0.1) next
      if(str_count(blocks.b$b[i.b], "n") / nchar(str_count(blocks.b$b[i.b], "n")) > 0.1) next
      
      if(length(blocks.q$b) >= length(blocks.b$b)){
        if(!is.null(diag.cutoff)){
          if((i.q - i.b) > diag.cutoff) next
          if(-(i.q - i.b) > length(blocks.q$b) + diag.cutoff) next
        }
      }
      # print(c(i.q, i.b))
      
      if(length(blocks.q$b) <= length(blocks.b$b)){
        if(!is.null(diag.cutoff)){
          if((i.b - i.q) > diag.cutoff) next
          if((i.b - i.q) > length(blocks.b$b) + diag.cutoff) next
        }
        
      }
      
      
      # print(c(i.q, i.b))
      
      localAlign <- pairwiseAlignment(blocks.q$b[i.q], blocks.b$b[i.b], 
                                      gapOpening = gap.open, gapExtension=gap.ext,
                                      type="local")
      pos.s.q = localAlign@pattern@range@start + blocks.q$p[i.q] - 1
      pos.e.q = localAlign@pattern@range@start + blocks.q$p[i.q] - 1 +
        localAlign@pattern@range@width - 1
      
      pos.s.b = localAlign@subject@range@start + blocks.b$p[i.b] - 1
      pos.e.b = localAlign@subject@range@start + blocks.b$p[i.b] - 1 +
        localAlign@subject@range@width - 1
      
      s.q.aln <- as.character(alignedPattern(localAlign))
      s.b.aln <- as.character(alignedSubject(localAlign))
      
      if(!rev.compl){
        dir.tmp = 0
      } else {
        dir.tmp = 1
      }
      
      t.tmp = rbind(t.tmp, 
                    data.frame(V1=paste0(c('betw', i.q, i.b), collapse = '|'),
                               V2=pos.s.q, V3=pos.e.q, V4=pos.s.b, V5=pos.e.b, V6=100, V7=nchar(s.q.aln), 
                               V8=s.q.aln, V9=s.b.aln, V10=t.arm[irow.q,'V10'], 
                               dir=dir.tmp, stringsAsFactors = FALSE))
      
      mx[i.q, i.b] = localAlign@pattern@range@width / nchar(blocks.q$b[i.q])
    }
    
  }
  t.tmp1 = t.tmp
  
  if(is.null(t.tmp)) return(list(t = t.tmp, mx = NULL))
  
  # Visializations
  
  # df <- melt(mx)
  # ggplot(data = df, aes(x=Var1, y=Var2, fill=value)) + 
  #   geom_tile() + theme_minimal()
  
  
  # Glue up things
  t.tmp = t.tmp1
  
  t.tmp[,c('V2','V3')] <- t.tmp[,c('V2','V3')] + t.arm[irow.q,'V3']
  if(rev.compl == FALSE){
    t.tmp[,c('V4','V5')] <- t.tmp[,c('V4','V5')] + pos.start
  } else {
    t.tmp[,c('V4','V5')] <- t.tmp[,c('V4','V5')]  + base.len - pos.start - length(s.b)
  }
  
  
  if(nrow(t.tmp) > 1) t.tmp = t.tmp[order(t.tmp[,'V2']),]
  t.tmp <- t.tmp[t.tmp[,'V7'] > 30,]
  # print(nrow(t.tmp))
  # showt(t.tmp)
  if(nrow(t.tmp) == 0)  return(list(t = NULL, mx = mx))
  if(nrow(t.tmp) > 1) {
    t.tmp <- glueZero(t.tmp)
    t.tmp = t.tmp[order(t.tmp[,'V2']),]
    t.tmp = glueByThreshold(t.tmp, block.len/2, query.fas = query.fas.chr, 
                            base.fas.fw = base.fas.fw, 
                            base.fas.bw = base.fas.bw, gap.open = gap.open)
    t.tmp <- removeCompleteOverlaps(t.tmp)
  }
  # showt(t.tmp)
  
  
  # plotSyntenyBlocks(t.tmp)
  

  
  checkCorrespToGenome(t.tmp, query.fas = query.fas.chr, 
                       base.fas.fw = base.fas.fw, 
                       base.fas.bw = base.fas.bw)
  return(list(t = t.tmp, mx = mx))
}


isIntersect <- function(x1, x2, y1, y2){
  return(!(( x2 < y1 ) || (x1 > y2)))
}


nIntersect <- function(x1, x2, y1, y2){
  if(!isIntersect(x1, x2, y1, y2)) return(0)
  return(min(x2 - y1 + 1, y2 - x1 + 1))
}

tIntersect <- function(t, i, j){
  return(nIntersect(t[i,'V2'], t[i,'V3'], t[j,'V2'], t[j,'V3']))
}


refineArms <- function(t.arm, block.len = 1000, diag.cutoff = NULL){
  
  # gap.open = 10
  # gap.ext = 0.5
  # irow = 6
  
  
  irow = 0
  while(irow < (nrow(t.arm)-1)){
    irow = irow + 1
    
    if(t.arm$dir[irow] == 1){
      rev.set = c(TRUE, FALSE)
    } else {
      rev.set = c(FALSE, TRUE)
    }
    for(rev.compl in rev.set){
      
      t.arm.base = t.arm
      tmp  = base.len - t.arm.base[t.arm.base$dir == 1,'V5'] + 1
      t.arm.base[t.arm.base$dir == 1,'V5'] = base.len - t.arm.base[t.arm.base$dir == 1,'V4'] + 1
      t.arm.base[t.arm.base$dir == 1,'V4'] = tmp
      
      
      # print(c(irow, rev.compl))
      
      irow.q = irow
      irow.b = irow
      
      if(t.arm[irow.q+1,'V2'] - t.arm[irow.q,'V3'] - 1 < block.len) next
      if(t.arm[irow.q+1,'V2'] - t.arm[irow.q,'V3'] - 1 > 100 * block.len) next
      
      if(t.arm.base[irow.b,'V4'] > t.arm.base[irow.b+1,'V5'] ){
        if(t.arm.base[irow.b,'V4'] - t.arm.base[irow.b+1,'V5'] - 1 < block.len) next
        if(t.arm.base[irow.b,'V4'] - t.arm.base[irow.b+1,'V5'] - 1 > 100 * block.len) next
      } else {
        if(t.arm.base[irow.b+1,'V4'] - t.arm.base[irow.b,'V5'] - 1 < block.len) next
        if(t.arm.base[irow.b+1,'V4'] - t.arm.base[irow.b,'V5'] - 1 > 100 * block.len) next
      }
      
      
      # Perform Alignment
      aln.blocks <- alignBlocks(t.arm, irow.q, irow.b, 
                                query.fas.chr, base.fas.fw,
                                base.fas.bw, block.len, rev.compl = rev.compl,
                                diag.cutoff = diag.cutoff)
      t.tmp = aln.blocks$t
      
      if(is.null(t.tmp)) next
      
      # # Visualisation
      # df <- melt(aln.blocks$mx)
      # ggplot(data = df, aes(x=Var1, y=Var2, fill=value)) +
      #   geom_tile() + theme_minimal()
      
      
      # Split with previous
      t.tmp = rbind(t.arm[c(irow, irow+1),], t.tmp)
      
      t.tmp = glueByThreshold(t.tmp, thresholds, query.fas = query.fas.chr, 
                              base.fas.fw = base.fas.fw, 
                              base.fas.bw = base.fas.bw, file.log=file.log, log.append = T)
      
      t.tmp = t.tmp[t.tmp[,'V7'] > 500,]
      
      if(sum(duplicated(rbind(t.tmp, t.arm[c(irow, irow+1),]))) == 2) next
      
      t.arm = rbind(t.arm[-c(irow, irow+1),], t.tmp)
      t.arm <- t.arm[order(t.arm[,'V2']),]
      checkCorrespToGenome(t.arm, query.fas = query.fas.chr, 
                           base.fas.fw = base.fas.fw, 
                           base.fas.bw = base.fas.bw)
      
      
      #showt(t.arm, irow:(irow+5))
      
      # irow = irow - 1 + nrow(t.tmp) - 1
      irow  = irow - 1
      break 
      

      # t.arm = afterRefinement(t.arm)
    }
  }
  return(t.arm)
}


afterRefinement <- function(t, min.len = 500){
  t = t[t[,'V7'] > min.len,]
  t <- removeCompleteOverlaps(t, n.distant = NULL, filter.base = F)
  t = t[order(t[, 'V2']),]
  t <- removeShortOverlaps(t, echo = F, gap.thresh = Inf)
  t = glueZero(t)
  t = t[order(t[,'V2']),]
  
  return(t)
}

findPositionInBase <- function(t, pos.search, base.len){
  
  t.base = t
  tmp  = base.len - t.base[t.base$dir == 1,'V5'] + 1
  t.base[t.base$dir == 1,'V5'] = base.len - t.base[t.base$dir == 1,'V4'] + 1
  t.base[t.base$dir == 1,'V4'] = tmp
  
  irows = which((t.base[,'V4'] <= pos.search) & (t.base[,'V5'] >= pos.search))
  if(length(irows) == 0) return(NA)
    
  pos.q = c()
  for(irow in irows){
    s1 = strsplit(t[irow,'V8'], '')[[1]]
    s2 = strsplit(t[irow,'V9'], '')[[1]]
    pos = rep(0, length(s1))
    pos1 = pos; pos1[s1 != '-'] = t[irow,'V2']: t[irow,'V3']
    if(t[irow, 'dir'] == 0){
      pos2 = pos; pos2[s2 != '-'] = t[irow,'V4']: t[irow,'V5']
    } else {
      pos2 = pos; pos2[s2 != '-'] = rev(t.base[irow,'V4']: t.base[irow,'V5'])
    }
    pos.q = c(pos.q, pos1[pos2 == pos.search] * (-1)^t[irow,'dir'])
    
  }
  if(is.null(pos.q)) return(NA)
  return(pos.q)  
  
}


writeSynteny <- function(t, file.log){
  write('', file=file.log)
  for(irow in 1:nrow(t)){
    write(paste0('>', t$V1[irow], '_query', collapse=''), file=file.log, append=T)
    write(t$V8[irow], file=file.log, append=T)
    write(paste0('>', t$V1[irow], '_base', collapse=''), file=file.log, append=T)
    write(t$V9[irow], file=file.log, append=T)
    write('', file=file.log, append=T)
  }
}


writeFasta <- function(s, file.log, append=FALSE, gap.line=T){
  if(!append) write('', file=file.log)
  for(i in 1:length(s)){
    write(paste0('>', names(s)[i], collapse=''), file=file.log, append=T)
    write(s[i], file=file.log, append=T)
    if(gap.line) write('', file=file.log, append=T)
  }
  if(!gap.line) write('', file=file.log, append=append)
}

getBase <- function(t, base.len){
  t.base = t
  tmp  = base.len - t.base[t.base$dir == 1,'V5'] + 1
  t.base[t.base$dir == 1,'V5'] = base.len - t.base[t.base$dir == 1,'V4'] + 1
  t.base[t.base$dir == 1,'V4'] = tmp
  
  return(t.base)
}

bestLocalAln <-function(s1, s2, gapExtension=0.05, gapOpening=5){
  localAlign1 <- pairwiseAlignment(s1, 
                                   s2, 
                                   gapExtension=gapExtension, gapOpening=gapOpening,
                                   type="local")
  
  localAlign2 <- pairwiseAlignment(s1, 
                                   reverseComplement(s2), 
                                   gapExtension=gapExtension, gapOpening=gapOpening,
                                   type="local")
  if(localAlign2@score > localAlign1@score){
    localAlign = localAlign2
  } else {
    localAlign = localAlign1
  }
  return(localAlign)
}

getGroups <- function(t){
  tbl = table(c(t[,'V1'], t[,'V8']))
  tbl = tbl[order(-tbl)]
  t.unique = names(tbl)
  t.group = rep(0, length(t.unique))
  names(t.group) = t.unique
  i.gr = 1
  
  while(length(t.unique) > 0){
    # print(length(t.unique))
    relatives = t.unique[1]
    n.relatives = 0
    while(n.relatives != length(relatives)){
      n.relatives = length(relatives)
      # print(n.relatives)
      relatives = unique(c(relatives, 
                           t$V8[t$V1 %in% relatives], 
                           t$V1[t$V8 %in% relatives]))
    }
    t.group[relatives] = i.gr
    i.gr = i.gr + 1
    t.unique = t.unique[!(t.unique %in% relatives)]
    t = t[(t$V1 %in% t.unique),]
  }
  return(t.group)
}

getSeqsFasta <- function(file){
  m.seq.types = read.fasta(file)
  m.seqs = sapply(names(m.seq.types), function(s) paste0(m.seq.types[[s]], collapse = ''), USE.NAMES = T)
  m.seqs = toupper(m.seqs)
  return(m.seqs)
}


getMobilome2 <- function(t, base.len, allowed_gap = 5, m.min.length = 15, 
                        thresh = 0.8, m.max.length = Inf, echo=F){
  
  m1 = c()
  m2 = c()
  pos.m1 = c()
  pos.m2 = c()
  for(irow in 1:nrow(t)){
    # if(echo) print(irow)
    s1 = strsplit(t[irow,'V8'], '')[[1]]
    s2 = strsplit(t[irow,'V9'], '')[[1]]
    
    pos = rep(0, length(s1))
    pos1 = pos; pos1[s1 != '-'] = t[irow,'V2']: t[irow,'V3']
    pos2 = pos; pos2[s2 != '-'] = t[irow,'V4']: t[irow,'V5']
    
    icol = 0
    n = length(s1)
    
    
    # Gaps in s1    
    all_ends = which((s1[1:(n-1)] == '-') & (s1[2:n] != '-'))
    all_starts = which((s1[1:(n-1)] != '-') & (s1[2:n] == '-')) + 1
    if(s1[1] == '-') all_starts = c(1, all_starts)
    if(s1[n] == '-') all_ends = c(all_ends, n)
    if(length(all_starts) != length(all_ends)) stop('problem')
    if(length(all_ends > 0)){
      gaps1 = cbind(all_starts, all_ends, 1)
      # print(nrow(gaps1))
      
      jrow = 2
      while (jrow <= nrow(gaps1)){
        if(gaps1[jrow,1] - gaps1[jrow-1,2] - 1 <= allowed_gap){
          gaps1[jrow-1, 2] = gaps1[jrow, 2]
          gaps1 = gaps1[-jrow,]
          if(is.null(nrow(gaps1))){
            gaps1 = c()
            break
          }
        } else {
          jrow = jrow + 1
        }
      }
      # print(nrow(gaps1))
    } else {
      gaps1 = c()
    }
    
    # Gaps in s2
    all_ends = which((s2[1:(n-1)] == '-') & (s2[2:n] != '-'))
    all_starts = which((s2[1:(n-1)] != '-') & (s2[2:n] == '-')) + 1
    if(s2[1] == '-') all_starts = c(1, all_starts)
    if(s2[n] == '-') all_ends = c(all_ends, n)
    if(length(all_starts) != length(all_ends)) stop('problem')
    if(length(all_ends > 0)){
      gaps2 = cbind(all_starts, all_ends, 2)
      # print(nrow(gaps2))
      jrow = 2
      while (jrow <= nrow(gaps2)){
        if(gaps2[jrow,1] - gaps2[jrow-1,2] - 1 <= allowed_gap){
          gaps2[jrow-1, 2] = gaps2[jrow, 2]
          gaps2 = gaps2[-jrow,]
          if(is.null(nrow(gaps2))){
            gaps2 = c()
            break
          }
        } else {
          jrow = jrow + 1
        }
      }
    } else {
      gaps2 = c()
    }
    # print(nrow(gaps2))
    
    
    gaps = rbind(gaps1, gaps2)
    if(is.null(gaps)) next
    gaps = gaps[order(gaps[,1]),, drop=F]
    gaps = gaps[gaps[,2] - gaps[,1] + 1 >= 15,, drop=F]
    
    if(is.null(gaps)) next
    if(nrow(gaps) == 0) next
    
    gaps[,2] - gaps[,1] + 1
    
    # print(c(irow, nrow(gaps)))
    
    # print(c(sum(gaps[,3] == 1), sum(gaps[,3] == 2)))
    # print(sum(gaps[,3] == 2))
    
    for(jrow in 1:nrow(gaps)){
      if(gaps[jrow, 3] == 1){ # start gaps on query
        i.start = gaps[jrow, 1]
        i.end = gaps[jrow, 2]
        s = toupper(s2[i.start:i.end])
        # print(length(s))
        # Nucleotide frequences
        f = table(s)
        if(max(f) > thresh * length(s))next
        
        # Di-Nucleotide frequences
        cf = apply(combn(f, 2), 2, sum)
        if(max(cf) > thresh * length(s)) next
        
        if(t[irow, 'dir'] == 1){
          pos2.start <- base.len - pos2[i.end] + 1
          pos2.end <- base.len - pos2[i.start] + 1
        } else {
          pos2.start <- pos2[i.start]
          pos2.end <- pos2[i.end]
        }
        
        if(i.start == 1){
          # print(dim(m2))
          if(length(c(pos2.start, pos2.end, 0, length(s), paste0(s, collapse = ''))) != 5) stop('ddd')
          m2 = rbind(m2, c(pos2.start, pos2.end, 0, length(s), paste0(s, collapse = ''))) 
        } else {
          # print(dim(m2))
          if(length(c(pos2.start, pos2.end, pos1[i.start-1], length(s), paste0(s, collapse = ''))) != 5) stop('eee')
          m2 = rbind(m2, c(pos2.start, pos2.end, pos1[i.start-1], length(s), paste0(s, collapse = ''))) 
        }
      } else {
        i.start = gaps[jrow, 1]
        i.end = gaps[jrow, 2]
        s = toupper(s1[i.start:i.end])
        # print(length(s))
        
        # Nucleotide frequences
        f = table(s)
        if(max(f) > thresh * length(s)) {
          # print(c('s:', s))
          next
        }
        
        # Di-Nucleotide frequences
        cf = apply(combn(f, 2), 2, sum)
        if(max(cf) > thresh * length(s)) {
          # print(c('s:', s))
          next
        }
        
        if(i.start == 1){
          # print(dim(m1))
          if(length(c(pos1[i.start], pos1[i.end], 0, length(s), paste0(s, collapse = ''))) != 5) stop('bb')
          m1 = rbind(m1, c(pos1[i.start], pos1[i.end], 0, length(s), paste0(s, collapse = '')))  
        } else {
          if(t[irow, 'dir'] == 1){
            pos2.start <- base.len - pos2[i.start-1] + 1
          } else {
            pos2.start = pos2[i.start-1]
          }
          # print(dim(m1))
          if(length(c(pos1[i.start], pos1[i.end], pos2.start, length(s), paste0(s, collapse = ''))) != 5) stop('aa')
          m1 = rbind(m1, c(pos1[i.start], pos1[i.end], pos2.start, length(s), paste0(s, collapse = '')))
        }        
      }
    } # for(jrow in nrow(gaps))
    # print('----')
    # print(c(dim(m1), dim(m2)))
    # print('----')
    
    # if('5892320' %in% m2[,1]) stop('')
  }

if(is.null(dim(m1)) || is.null(dim(m1))) return(NULL)
#print(dim(m1))
#print(dim(m2))
  
  colnames(m1) <- c('pos.q.start', 'pos.q.end', 'pos.b.before', 'length', 'seq')
  colnames(m2) <- c('pos.b.start', 'pos.b.end', 'pos.q.before', 'length', 'seq')
  
  m1 = as.data.frame(m1)
  m2 = as.data.frame(m2)
  
  for(i in 1:4){
    m1[,i] =  as.numeric(as.character(m1[,i]))
    m2[,i] =  as.numeric(as.character(m2[,i]))
  }
  
  colnames(m1) <- c('pos.q.start', 'pos.q.end', 'pos.b.before', 'length', 'seq')
  colnames(m2) <- c('pos.b.start', 'pos.b.end', 'pos.q.before', 'length', 'seq')
  
  m1 = as.data.frame(m1)
  m2 = as.data.frame(m2)
  
  for(i in 1:4){
    m1[,i] =  as.numeric(as.character(m1[,i]))
    m2[,i] =  as.numeric(as.character(m2[,i]))
  }
  
  return(list(query = m1, base = m2))
}


getAlnByPositionQuery <- function(t, pos1, pos2){
  idx1 = which(pos1 >= t$V2)
  idx2 = which(pos2 <= t$V3)
  idx = intersect(idx1, idx2)
  if(length(idx) == 0){
    message('No block found')
    return(c(NA, NA))
  }
  if(length(idx) > 1){
    message('Several blocks found')
    return(c(NA, NA))
  }
  
  positions = rep(0, t$V7[idx])
  positions[strsplit(t$V8[idx],'')[[1]] != '-'] = t$V2[idx]:t$V3[idx]
  pos.aln = c(which(positions == pos1), which(positions == pos2))
  seq1 = substr(t$V8[idx],pos.aln[1],pos.aln[2])
  seq2 = substr(t$V9[idx],pos.aln[1],pos.aln[2])
  return(c(seq1, seq2))
}

getAlnByPositionBase <- function(t, pos1, pos2){
  idx1 = which(pos1 >= t$V4)
  idx2 = which(pos2 <= t$V5)
  idx = intersect(idx1, idx2)
  if(length(idx) == 0){
    message('No block found')
    return(c(NA, NA))
  }
  if(length(idx) > 1){
    message('Several blocks found')
    return(c(NA, NA))
  }
  
  positions = rep(0, t$V7[idx])
  positions[strsplit(t$V9[idx],'')[[1]] != '-'] = t$V4[idx]:t$V5[idx]
  pos.aln = c(which(positions == pos1), which(positions == pos2))
  seq1 = substr(t$V8[idx],pos.aln[1],pos.aln[2])
  seq2 = substr(t$V9[idx],pos.aln[1],pos.aln[2])
  return(c(seq1, seq2))
}

getPositionBaseByPositionQuery <- function(t, pos1, pos2){
  idx1 = which(pos1 >= t$V2)
  idx2 = which(pos2 <= t$V3)
  idx = intersect(idx1, idx2)
  if(length(idx) == 0){
    message('No block found')
    return(c(NA, NA))
  }
  if(length(idx) > 1){
    message('Several blocks found')
    return(c(NA, NA))
  }
  
  positions.q = rep(0, t$V7[idx])
  positions.q[strsplit(t$V8[idx],'')[[1]] != '-'] = t$V2[idx]:t$V3[idx]
  
  positions.b = rep(0, t$V7[idx])
  positions.b[strsplit(t$V9[idx],'')[[1]] != '-'] = t$V4[idx]:t$V5[idx]
  
  
  pos.aln = c(which(positions.q == pos1), which(positions.q == pos2))
  # print(pos.aln)
  
  pos.b1 = max(positions.b[1:pos.aln[1]])
  positions.b[positions.b == 0] = +Inf
  pos.b2 = min(positions.b[pos.aln[2]:t$V7[idx]])
  
  return(c(pos.b1, pos.b2))
}

getPositionBaseByPositionQuery1 <- function(t, pos, base.len, query.len, echo=F){
  
  t.base = getBase(t, base.len)
  
  pos.corresp = rep(0, query.len)
  for(irow in 1:nrow(t)){
    if(echo) print(irow)
    positions.q = rep(0, t$V7[irow])
    positions.q[strsplit(t$V8[irow],'')[[1]] != '-'] = t$V2[irow]:t$V3[irow]
    positions.b = rep(0, t$V7[irow])
    if(t$dir[irow] == 0){
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V4[irow]:t.base$V5[irow]  
    } else {
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V5[irow]:t.base$V4[irow]
    }
    
    
    positions.b = positions.b[positions.q != 0]
    pos.corresp[t$V2[irow]:t$V3[irow]] = positions.b
    # 
    # positions.q.all = c(positions.q.all, positions.q)
    # positions.b.all = c(positions.b.all, positions.b)
  }
  return(pos.corresp[pos])
}


getPositionBaseByPositionQuery2 <- function(t, pos1, pos2, base.len, query.len, echo=F){
  
  t.base = getBase(t, base.len)
  
  pos.corresp = rep(0, query.len)
  for(irow in 1:nrow(t)){
    if(echo) print(irow)
    positions.q = rep(0, t$V7[irow])
    positions.q[strsplit(t$V8[irow],'')[[1]] != '-'] = t$V2[irow]:t$V3[irow]
    positions.b = rep(0, t$V7[irow])
    if(t$dir[irow] == 0){
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V4[irow]:t.base$V5[irow]  
    } else {
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V5[irow]:t.base$V4[irow]
    }
    
    
    positions.b = positions.b[positions.q != 0]
    pos.corresp[t$V2[irow]:t$V3[irow]] = positions.b
    # 
    # positions.q.all = c(positions.q.all, positions.q)
    # positions.b.all = c(positions.b.all, positions.b)
  }
  
  idx = pos1 > pos2
  tmp = pos1[idx]
  pos1[idx] = pos2[idx]
  pos2[idx] = tmp
  
  pos.base = c()
  for(i in 1:length(pos1)){
    pos.tmp = pos.corresp[pos1[i]:pos2[i]]
    pos.tmp = pos.tmp[pos.tmp!=0]
    if(length(pos.tmp) == 0) {
      x1 = pos.corresp[max(which(pos.corresp[1:pos1[i]] != 0))]
      x2 = pos.corresp[min(which(pos.corresp[pos2[i]:length(pos.corresp)] != 0)) + pos2[i] - 1]
      if(x1 <= x2){
        pos.base = rbind(pos.base, c(x1, x2, 0))
      } else {
        pos.base = rbind(pos.base, c(x2, x1, 1))
      }
      # pos.base = rbind(pos.base, c(0,0, 0))
    } else {
      pos.base = rbind(pos.base, c(min(pos.tmp),max(pos.tmp), pos.tmp[length(pos.tmp)] < pos.tmp[1]))
    }
  }
  return(pos.base)
}


getPositionBaseByPositionQuery3 <- function(t, pos1, pos2, base.len, query.len, echo=F){
  
  t.base = getBase(t, base.len)
  
  pos.corresp = rep(0, query.len)
  for(irow in 1:nrow(t)){
    if(echo) print(irow)
    positions.q = rep(0, t$V7[irow])
    positions.q[strsplit(t$V8[irow],'')[[1]] != '-'] = t$V2[irow]:t$V3[irow]
    positions.b = rep(0, t$V7[irow])
    if(t$dir[irow] == 0){
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V4[irow]:t.base$V5[irow]  
    } else {
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V5[irow]:t.base$V4[irow]
    }
    
    idx = positions.q != 0
    positions.b = positions.b[idx]
    positions.q = positions.q[idx]
    
    pos.corresp[positions.q] = positions.b
    if(sum(is.na(pos.corresp))) stop()
  }
  
  idx = pos1 > pos2
  if(sum(idx) > 0){
    tmp = pos1[idx]
    pos1[idx] = pos2[idx]
    pos2[idx] = tmp
  }
  

  pos.tmp = -Inf
  pos.prev = pos.corresp
  for(i in 1:length(pos.corresp)){
    if(pos.corresp[i] == 0) {
      pos.prev[i] = pos.tmp
    } else {
      pos.tmp = max(pos.prev[i], pos.tmp)
    }
  }
  pos.tmp = +Inf
  pos.next = pos.corresp
  for(i in length(pos.corresp):1){
    if(pos.corresp[i] == 0) {
      pos.next[i] = pos.tmp
    } else {
      pos.tmp = min(pos.next[i], pos.tmp)
    }
  }
  # pos.base = cbind(pos.prev[pos1], pos.next[pos2])
  
  pos.base = cbind(pos.corresp[pos1], pos.corresp[pos2])
  pos.base = cbind(pos.base, 0)
  idx = which((pos.base[,1] == 0) | (pos.base[,2] == 0))
  idx.good = setdiff(1:nrow(pos.base), idx)

  idx.reverse = idx.good[pos.base[idx.good,1] > pos.base[idx.good,2]]
  pos.base[idx.reverse,3] = 1
  for(i in idx){
    if(echo) print(i)

    if((pos.next[pos2[i]] == pos.next[pos1[i]]) & (pos.prev[pos2[i]] == pos.prev[pos1[i]])) {
      x1 = pos.prev[pos1[i]]
      x2 = pos.next[pos2[i]]
      if(x1 <= x2){
        pos.base[i,] = c(x1, x2, 0)
      } else {
        pos.base[i,] = c(x2, x1, 1)
      }
      # pos.base = rbind(pos.base, c(0,0, 0))
    } else {
      pos.tmp = c(pos.prev[pos1[i]], pos.next[pos2[i]])
      pos.base[i,] = c(min(pos.tmp),max(pos.tmp), pos.tmp[length(pos.tmp)] < pos.tmp[1])
    }
  }
  
  
  # pos.base = cbind(pos.corresp[pos1], pos.corresp[pos2])
  # pos.base = cbind(pos.base, 0)
  # idx = which((pos.base[,1] == 0) | (pos.base[,2] == 0))
  # idx.good = setdiff(1:nrow(pos.base), idx)
  # 
  # idx.reverse = idx.good[pos.base[idx.good,1] > pos.base[idx.good,2]]
  # pos.base[idx.reverse,3] = 1
  # for(i in idx){
  #   if(echo) print(i)
  #   pos.tmp = pos.corresp[pos1[i]:pos2[i]]
  #   pos.tmp = pos.tmp[pos.tmp!=0]
  #   if(length(pos.tmp) == 0) {
  #     x1 = pos.corresp[max(which(pos.corresp[1:pos1[i]] != 0))]
  #     x2 = pos.corresp[min(which(pos.corresp[pos2[i]:length(pos.corresp)] != 0)) + pos2[i] - 1]
  #     if(is.na(x1)) x1 = -Inf
  #     if(is.na(x2)) x2 = Inf
  #     if(x1 <= x2){
  #       pos.base[i,] = c(x1, x2, 0)
  #     } else {
  #       pos.base[i,] = c(x2, x1, 1)
  #     }
  #     # pos.base = rbind(pos.base, c(0,0, 0))
  #   } else {
  #     pos.base[i,] = c(min(pos.tmp),max(pos.tmp), pos.tmp[length(pos.tmp)] < pos.tmp[1])
  #   }
  # }
  return(pos.base)
}

getPositionQueryByPositionBase2 <- function(t, pos1, pos2, base.len, query.len, echo=F){
  
  t.base = getBase(t, base.len)
  
  pos.corresp = rep(0, base.len)
  for(irow in 1:nrow(t)){
    if(echo) print(irow)
    positions.q = rep(0, t$V7[irow])
    positions.q[strsplit(t$V8[irow],'')[[1]] != '-'] = t$V2[irow]:t$V3[irow]
    positions.b = rep(0, t$V7[irow])
    if(t$dir[irow] == 0){
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V4[irow]:t.base$V5[irow]  
    } else {
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V5[irow]:t.base$V4[irow]
    }
    
    idx = positions.b != 0
    positions.b = positions.b[idx]
    positions.q = positions.q[idx]
    
    pos.corresp[positions.b] = positions.q
  }
  
  idx = pos1 > pos2
  tmp = pos1[idx]
  pos1[idx] = pos2[idx]
  pos2[idx] = tmp
  
  pos.base = cbind(pos.corresp[pos1], pos.corresp[pos2])
  pos.base = cbind(pos.base, 0)
  idx = which((pos.base[,1] == 0) | (pos.base[,2] == 0))
  idx.good = setdiff(1:nrow(pos.base), idx)
  
  idx.reverse = idx.good[pos.base[idx.good,1] > pos.base[idx.good,2]]
  pos.base[idx.reverse,3] = 1
  
  for(i in idx){
    if(echo) print(i)
    pos.tmp = pos.corresp[pos1[i]:pos2[i]]
    pos.tmp = pos.tmp[pos.tmp!=0]
    if(length(pos.tmp) == 0) {
      x1 = pos.corresp[max(which(pos.corresp[1:pos1[i]] != 0))]
      x2 = pos.corresp[min(which(pos.corresp[pos2[i]:length(pos.corresp)] != 0)) + pos2[i] - 1]
      if(is.na(x1)) x1 = -Inf
      if(is.na(x2)) x2 = Inf
      if(x1 <= x2){
        pos.base[i,] = c(x1, x2, 0)
      } else {
        pos.base[i,] = c(x2, x1, 1)
      }
      # pos.base = rbind(pos.base, c(0,0, 0))
    } else {
      pos.base[i,] = c(min(pos.tmp),max(pos.tmp), pos.tmp[length(pos.tmp)] < pos.tmp[1])
    }
  }
  return(pos.base)
}

getPositionQueryByPositionBase1 <- function(t, pos, base.len, query.len, echo=F){
  
  t.base = getBase(t, base.len)
  
  pos.corresp = rep(0, query.len)
  for(irow in 1:nrow(t)){
    if(echo) print(irow)
    positions.q = rep(0, t$V7[irow])
    positions.q[strsplit(t$V8[irow],'')[[1]] != '-'] = t$V2[irow]:t$V3[irow]
    positions.b = rep(0, t$V7[irow])
    if(t$dir[irow] == 0){
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V4[irow]:t.base$V5[irow]  
    } else {
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V5[irow]:t.base$V4[irow]
    }
    
    positions.q = positions.q[positions.b != 0]
    pos.corresp[t.base$V4[irow]:t.base$V5[irow]] = positions.q
    # 
    # positions.q.all = c(positions.q.all, positions.q)
    # positions.b.all = c(positions.b.all, positions.b)
  }
  
  
  return(pos.corresp[pos])
}


getAlignmentByPositionBase1 <- function(t, pos1, pos2, base.len, echo=F){
  
  t.base = getBase(t, base.len)
  
  pos.corresp = rep(0, query.len)
  ipos1 = which((t.base$V4 <= pos1) & (t.base$V5 >= pos1) )
  ipos2 = which((t.base$V4 <= pos2) & (t.base$V5 >= pos2) )
  
  if(ipos1 == ipos2){
    irow = ipos1
    positions.q = rep(0, t$V7[irow])
    positions.q[strsplit(t$V8[irow],'')[[1]] != '-'] = t$V2[irow]:t$V3[irow]
    positions.b = rep(0, t$V7[irow])
    if(t$dir[irow] == 0){
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V4[irow]:t.base$V5[irow]  
    } else {
      positions.b[strsplit(t$V9[irow],'')[[1]] != '-'] = t.base$V5[irow]:t.base$V4[irow]
    }
    
    idx1 = which(positions.b == pos1)
    idx2 = which(positions.b == pos2)
    
    ipos1 = positions.q[idx1]
    ipos2 = positions.q[idx2]
    
    message(paste(ipos1, ipos2, sep = ':'))
    
    if(idx2 < idx1){
      tmp = idx1
      idx1 = idx2
      idx2 = tmp
    }
    
    s1 = paste0(strsplit(t$V8[irow], '')[[1]][idx1:idx2], sep = '', collapse = '')
    s2 = paste0(strsplit(t$V9[irow], '')[[1]][idx1:idx2], sep = '', collapse = '')
    
    return(list(s1 = s1, s2 = s2))
  } else {
    message('Additional code is required')
  }
  
  return(NULL)
}

getBeginEndOfGaps <- function(pos.query.free){
  query.len = length(pos.query.free)
  pos.query.beg <- which((pos.query.free[-1] == 0) & (pos.query.free[-query.len] == 1)) + 1
  pos.query.end <- which((pos.query.free[-1] == 1) & (pos.query.free[-query.len] == 0)) # - 1
  if(pos.query.free[1] != 1) pos.query.beg = c(1, pos.query.beg)
  if(pos.query.free[query.len] != 1) pos.query.end = c(pos.query.end, query.len)

  if(length(pos.query.end) != length(pos.query.beg)) stop('Problem with begind and ends: 1')
  # pos.query.end[pos.query.end - pos.query.beg < 0]
  # pos.query.beg[pos.query.end - pos.query.beg < 0]
  
  if(min(pos.query.end - pos.query.beg) < -1) stop('Problem with begind and ends: 2')
  idx = pos.query.end - pos.query.beg > 15
  
  pos.query.end = pos.query.end[idx]
  pos.query.beg = pos.query.beg[idx]
  return(list(beg = pos.query.beg, end = pos.query.end))
}

getFeaturesFromRow <- function(t, irow, acc.name, base.name, gap.window = 20, snp.window = 5,
                               indel.len = 15){

  s1 = strsplit(t$V8[irow], '')[[1]]
  s2 = strsplit(t$V9[irow], '')[[1]]
  s.dif = s1 == s2
  n = length(s.dif)
  
  pos1 = rep(0,n)
  pos1[s1 != '-'] = t[irow, 'V2']:t[irow, 'V3']
  pos2 = rep(0,n)
  pos2[s2 != '-'] = t[irow, 'V4']:t[irow, 'V5']
  
  all_ends = which((s.dif[1:(n-1)] == F) & (s.dif[2:n] != F))
  all_starts = which((s.dif[1:(n-1)] != F) & (s.dif[2:n] == F)) + 1
  if(s.dif[1] == F) all_starts = c(1, all_starts)
  if(s.dif[n] == F) all_ends = c(all_ends, n)
  if(length(all_starts) != length(all_ends)) stop('problem')
  
  
  gaps <- cbind(all_starts, all_ends, all_ends - all_starts + 1)
  gaps <- cbind(gaps, gaps[,3] >= indel.len)
  idx.gaps <- which(gaps[,4] == 1)
  idx.gaps <- idx.gaps[order(gaps[idx.gaps,3])]
  
  if(nrow(gaps) == 0){
    return(list(snps=NULL, indels=NULL))
  }
  
  gaps.add <- c()
  for(igap in idx.gaps){
    if(is.na(gaps[igap, 1])) next
    gap.dist = min(round(gaps[igap, 3] / 2), gap.window)
    gap.range = gaps[igap, 1:2] + c(-gap.dist, gap.dist)
    idx.in.gap = which(((gaps[,1] >= gap.range[1]) & (gaps[,1] <= gap.range[2])) | 
                         ((gaps[,2] >= gap.range[1]) & (gaps[,2] <= gap.range[2])))
    # idx.in.gap = setdiff(idx.in.gap, igap)
    idx.in.gap = intersect(idx.in.gap, idx.gaps)
    if(length(idx.in.gap) == 1) next
    # break
    tmp <- c(min(gaps[idx.in.gap,1]), max(gaps[idx.in.gap,2]))
    gaps.add <- rbind(gaps.add, c(tmp, tmp[2] - tmp[1], 1,igap))
    gaps[idx.in.gap,] <- NA
    gaps <- rbind(gaps, c(tmp, tmp[2] - tmp[1] + 1, 1))
  }
  
  gaps = gaps[!is.na(gaps[,1]),, drop=F]
  gaps = gaps[order(gaps[,1]),, drop=F]
  
  
  # remove overlaps
  idx.overlap <- c()
  for(igap in 1:nrow(gaps)){
    idx.tmp = (gaps[igap,1] >= gaps[,1]) & (gaps[igap,2] <= gaps[,2])
    if(sum(idx.tmp) == 1) next
    # idx.overlap <- c(idx.overlap, setdiff(which(idx.tmp), igap))
    idx.overlap <- c(idx.overlap, igap)
  }
  if(length(idx.overlap) > 0) gaps <- gaps[-idx.overlap,, drop=F]
  
  
  # add SNPs in distance
  idx.add.snp <- c()
  for(igap in 1:nrow(gaps)){
    if(gaps[igap,4] == 1) next
    idx.tmp = which((gaps[igap,1] >= (gaps[,1] - snp.window)) & (gaps[igap,2] <= (gaps[,2] + snp.window)) & 
                      (gaps[,4] == 1))
    if(length(idx.tmp) == 0) next
    idx.tmp = idx.tmp[gaps[idx.tmp,3] == max(gaps[idx.tmp,3])]
    idx.tmp = idx.tmp[1]
    idx.add.snp <- rbind(idx.add.snp, c(igap, idx.tmp))
  }
  
  unique.add <- unique(idx.add.snp[,2])
  for(igap in unique.add){
    idx.group <- unique(c(idx.add.snp[idx.add.snp[,2] == igap,1], igap))
    tmp <- c(min(gaps[idx.group,1], na.rm=TRUE), max(gaps[idx.group,2], na.rm=TRUE))
    gaps <- rbind(gaps, c(tmp, tmp[2] - tmp[1] + 1, 1))
    gaps[idx.group,] <- NA
  }
  
  gaps = gaps[!is.na(gaps[,1]),, drop=F]
  gaps = gaps[order(gaps[,1]),, drop=F]
  
  
  gaps.indel <- gaps[gaps[,4] == 1, ,drop=F]
  ng <- nrow(gaps)
  if(ng == 1){
    gaps <- cbind(gaps, cbind(Inf, Inf))
    
  } else {
    gaps <- cbind(gaps, c(gaps[2:ng,2] - gaps[1:(ng-1),1], Inf))
    gaps <- cbind(gaps, c(Inf, gaps[2:ng,2] - gaps[1:(ng-1),1]))
  }
  
  
  
  gaps.snps = gaps[(gaps[,3] == 1) & (gaps[,5] > snp.window) & (gaps[,6] > snp.window), , drop=F]
  
  seqs.indels <- c()
  seqs.indels.names <- c()
  for(igap in 1:nrow(gaps.indel)){
    if(igap > nrow(gaps.indel)) break
    s1.tmp = s1[gaps.indel[igap,1]:gaps.indel[igap,2]]
    s2.tmp = s2[gaps.indel[igap,1]:gaps.indel[igap,2]]
    
    pos1.tmp = pos1[c((gaps.indel[igap,1]-1),(gaps.indel[igap,2]+1))]
    pos2.tmp = pos2[c((gaps.indel[igap,1]-1),(gaps.indel[igap,2]+1))]
    if(sum(c(pos1.tmp, pos2.tmp) == 0) != 0) stop('AAAA')
    
    s1.tmp = paste0(s1.tmp[s1.tmp!='-'], collapse = '')
    s2.tmp = paste0(s2.tmp[s2.tmp!='-'], collapse = '')
    
    if((nchar(s1.tmp) < indel.len) | (nchar(s2.tmp) < indel.len)) {
      indel.type='single'
    } else {
      indel.type='double'
    }
    if(nchar(s1.tmp) >= indel.len){
      seqs.indels <- c(seqs.indels, s1.tmp)
      name.tmp <- paste0(acc.name, '_', base.name, '|' ,
                         paste0(pos1.tmp, sep='|', collapse = ''),
                         paste0(pos2.tmp, sep='|', collapse = ''), 'q', '|', indel.type, collapse = '')
      seqs.indels.names <- c(seqs.indels.names, name.tmp)
    }
    if(nchar(s2.tmp) >= indel.len){
      seqs.indels <- c(seqs.indels, s2.tmp)
      name.tmp <- paste0(acc.name, '_', base.name, '|' ,
                         paste0(pos1.tmp, sep='|', collapse = ''),
                         paste0(pos2.tmp, sep='|', collapse = ''), 'b', '|', indel.type, collapse = '')
      seqs.indels.names <- c(seqs.indels.names, name.tmp)
    }
  }
  
  names(seqs.indels) <- seqs.indels.names
  
  gaps.snps <- cbind(gaps.snps, cbind(pos1[gaps.snps[,1]] , pos2[gaps.snps[,2]]))
  gaps.snps <- cbind(gaps.snps, cbind(s1[gaps.snps[,1]] , s2[gaps.snps[,2]]))
  gaps.snps <- gaps.snps[rowSums(gaps.snps == '-') == 0,, drop=F]
  if(sum(gaps.snps[,ncol(gaps.snps)] == gaps.snps[,ncol(gaps.snps) - 1]) != 0) stop('BBBBBBB')
  gaps.snps <- as.data.frame(gaps.snps[,(ncol(gaps.snps)-3):ncol(gaps.snps),drop=F], stringsAsFactors=F)
  gaps.snps[,1] <- as.numeric(gaps.snps[,1])
  gaps.snps[,2] <- as.numeric(gaps.snps[,2])
  return(list(snps=gaps.snps, indels=seqs.indels))
}

cleanTails <- function(t, mismatch.window = 4, echo=F){
  for(irow in 1:nrow(t)){
    # if(echo) print(irow)
    s1 = strsplit(t$V8[irow], '')[[1]]
    s2 = strsplit(t$V9[irow], '')[[1]]
    s.dif = (s1 != s2)*1
    
    s.dif.cumm <- s.dif
    for(istep in 1:(mismatch.window-1) ){
      s.dif.cumm <- s.dif.cumm + c(s.dif[-(1:istep)],rep(0, istep))
    }
    
    idx.dif <- which(s.dif.cumm == 0)
    idx.start <- min(idx.dif)
    
    if(idx.start != 1){
      if(echo){
        print(c('beg',irow))
        # break
        print(remainFirstN(t[irow, 'V8'], (idx.start-1)))
        print(remainFirstN(t[irow, 'V9'], (idx.start-1)))
      }
      
      n.nt.start.q <- sum(s1[1:(idx.start-1)] != '-')
      n.nt.start.b <- sum(s2[1:(idx.start-1)] != '-')
      t[irow, 'V2'] <- t[irow, 'V2'] + n.nt.start.q
      t[irow, 'V4'] <- t[irow, 'V4'] + n.nt.start.b
      t[irow, 'V8'] <- removeFirstN(t[irow, 'V8'], (idx.start-1))
      t[irow, 'V9'] <- removeFirstN(t[irow, 'V9'], (idx.start-1))
    } 
    
    s.dif.cumm <- s.dif
    for(istep in 1:(mismatch.window-1) ){
      s.dif.cumm <- s.dif.cumm + c(rep(0, istep),s.dif[1:(length(s.dif)-istep)])
    }
    idx.dif <- which(s.dif.cumm == 0)
    idx.end <- max(idx.dif)
    
    n.seq = length(s1)
    if(idx.end != n.seq){
      if(echo){
        print(c('end',irow))
        print(remainLastN(t[irow, 'V8'], (n.seq - idx.end)))
        print(remainLastN(t[irow, 'V9'], (n.seq - idx.end)))
      }
      # break 
      n.nt.end.q <- sum(s1[(idx.end+1):n.seq] != '-')
      n.nt.end.b <- sum(s2[(idx.end+1):n.seq] != '-')
      t[irow, 'V3'] <- t[irow, 'V3'] - n.nt.end.q
      t[irow, 'V5'] <- t[irow, 'V5'] - n.nt.end.b
      t[irow, 'V8'] <- removeLastN(t[irow, 'V8'], (n.seq - idx.end))
      t[irow, 'V9'] <- removeLastN(t[irow, 'V9'], (n.seq - idx.end))
      
      
    } 
    t[irow, 'V7'] <- nchar(t[irow, 'V8'])
  }
  return(t)
}

getIndelsBetween <- function(t, base.fas.fw, query.fas, acc.name, base.name,
                             indel.len = 15, indel.len.max = 30000){
  
  
  t <- t[order(t$V2),]
  t.base <- getBase(t, base.len)
  seqs.indels <- c()
  seqs.indels.names <- c()
  for(irow in 1:(nrow(t)-1)){
    if(t.base$dir[irow] != t.base$dir[irow+1]) next
    
    if(t.base$dir[irow] == 0){
      jrow = min(which(t.base$V4 > t.base$V5[irow]))
      # print(c(irow, jrow))
      if(jrow != (irow+1) ) next
      pos1 = c(t.base[irow, 'V3'], t.base[irow+1, 'V2'])
      pos2 = c(t.base[irow, 'V5'], t.base[irow+1, 'V4'])
      
      
    } else {
      tmp = abs(t.base$V4[irow] - t.base$V5)
      jrow = which(tmp == min(tmp))
      # print(c(irow, jrow))
      if(jrow != (irow+1) ) next
      pos1 = c(t.base[irow, 'V3'], t.base[irow+1, 'V2'])
      pos2 = c(t.base[irow+1, 'V4'], t.base[irow, 'V4'])
    }
    
    # print(irow)
    d1 = pos1[2] - pos1[1]-1
    d2 = pos2[2] - pos2[1]-1
    
    if((d1 < indel.len) | (d2 < indel.len)) {
      indel.type='single'
    } else {
      indel.type='double'
    }  
    
    if((d1 >= indel.len) & (d1 <= indel.len.max)){
      s1.tmp = toupper(paste(query.fas[(pos1[1]+1):(pos1[2]-1)], collapse = ''))
      seqs.indels <- c(seqs.indels, s1.tmp)
      name.tmp <- paste0(acc.name, '_', base.name, '|' ,
                         paste0(pos1, sep='|', collapse = ''),
                         paste0(pos2, sep='|', collapse = ''), 'q_bw', '|', indel.type, collapse = '')
      seqs.indels.names <- c(seqs.indels.names, name.tmp)
    }
    if((d2 >= indel.len) & (d2 <= indel.len.max)){
      s2.tmp = toupper(paste(base.fas.fw[(pos2[1]+1):(pos2[2]-1)],collapse = ''))
      seqs.indels <- c(seqs.indels, s2.tmp)
      name.tmp <- paste0(acc.name, '_', base.name, '|' ,
                         paste0(pos1, sep='|', collapse = ''),
                         paste0(pos2, sep='|', collapse = ''), 'b_bw', '|', indel.type, collapse = '')
      seqs.indels.names <- c(seqs.indels.names, name.tmp)
    }  
  }
  names(seqs.indels) <- seqs.indels.names
  return(seqs.indels)
}

removeOverlapsExtra <- function(t, query.len, base.len, to.stop=F, n.over = 20){
  t.base = getBase(t, base.len)
  pos.query = rep(0, query.len)
  pos.base = rep(0, base.len)
  for(irow in 1:nrow(t.base)){
    pos.query[t.base$V2[irow]:t.base$V3[irow]] = pos.query[t.base$V2[irow]:t.base$V3[irow]] + 1
    pos.base[t.base$V4[irow]:t.base$V5[irow]] = pos.base[t.base$V4[irow]:t.base$V5[irow]] + 1
  }
  
  idx.over.q = which(pos.query > 1)
  idx.over.b = which(pos.base > 1)
  
  irow.over.q = c()
  irow.over.b = c()
  for(irow in 1:nrow(t.base)){
    if(length(intersect(t.base$V2[irow]:t.base$V3[irow], idx.over.q)) > 0) irow.over.q <- c(irow.over.q, irow)
    if(length(intersect(t.base$V4[irow]:t.base$V5[irow], idx.over.b)) > 0) irow.over.b <- c(irow.over.b, irow)
  }

  if(to.stop){  
    print(irow.over.q)
    print(irow.over.b)
    if(!is.null(irow.over.q)) stop('aaa')
    if(!is.null(irow.over.b)) stop('bbb')
  }
  
  idx.remove = c()
  
  # print(irow.over.b)
  if(!is.null(irow.over.b)){
    irow.over = irow.over.b
    V2 = 'V4'
    V3 = 'V5'
    for(i.q1 in irow.over){
      for(i.q2 in irow.over){
        if(i.q1 == i.q2) next
        
        p1 = max(t.base[i.q1, V2], t.base[i.q2, V2]) 
        p2 = min(t.base[i.q1, V3], t.base[i.q2, V3])
          
        # if intersection
        if(p1 < p2){
          # print(p2 - p1)
          if((p2 - p1) < n.over) next
          
          # print(c(i.q1, i.q2))
          # print('yes')
          l1 = t.base[i.q1, V3] - t.base[i.q1, V2]
          l2 = t.base[i.q2, V3] - t.base[i.q2, V2]
          if(l1 > l2){
            idx.remove = c(idx.remove, i.q2)
          } else if (l1 < l2){
            idx.remove = c(idx.remove, i.q1)
          } else {
            idx.remove = c(idx.remove, c(i.q2, i.q1))
          }
        }
      }
    }
  }
  
  # print(irow.over.q)
  if(!is.null(irow.over.q)){
    irow.over = irow.over.q
    V2 = 'V2'
    V3 = 'V3'
    idx.remove = c()
    for(i.q1 in irow.over){
      for(i.q2 in irow.over){
        if(i.q1 == i.q2) next
        
        # if intersection
        if(max(t.base[i.q1, V2], t.base[i.q2, V2]) < min(t.base[i.q1, V3], t.base[i.q2, V3]) ){
          # print(c(i.q1, i.q2))
          # print('yes')
          l1 = t.base[i.q1, V3] - t.base[i.q1, V2]
          l2 = t.base[i.q2, V3] - t.base[i.q2, V2]
          if(l1 > l2){
            idx.remove = c(idx.remove, i.q2)
          } else if (l1 < l2){
            idx.remove = c(idx.remove, i.q1)
          } else {
            idx.remove = c(idx.remove, c(i.q2, i.q1))
          }
        }
      }
    }
  }
  
  if(!is.null(idx.remove)){
    idx.remove = unique(idx.remove)
    print('remove')
    print(idx.remove)
    t = t[-idx.remove,]
  }
  
  return(t)
}

hardCleaning <- function(t, base.len, len.min = 1000){
  t = t[t$V7 >= len.min,]
  
  irow.delete <- c()
  for(irow in 1:nrow(t)){
    s1 = strsplit(t$V8[irow], '')[[1]]
    s2 = strsplit(t$V9[irow], '')[[1]]
    if(sum(s1 != '-') == 0) irow.delete <- c(irow.delete, irow)
    if(sum(s2 != '-') == 0) irow.delete <- c(irow.delete, irow)
    if(sum(s2 != s1)/length(s1) > 0.3) irow.delete <- c(irow.delete, irow)
  }
  if(length(irow.delete) > 0) t <- t[-unique(irow.delete),]
  
  t <- cleanTails(t, echo=F)
  
  #Remove duplicated
  pos.dupl <- unique(t[duplicated(t[,c('V2', 'V3')]), c('V2', 'V3'), drop=F])
  # print(pos.dupl)
  for(idupl in 1:nrow(pos.dupl)){
    if(idupl > nrow(pos.dupl)) break
    t <- t[!((t$V2 == pos.dupl[idupl,1]) & ((t$V3 == pos.dupl[idupl,2]))),]
  }
  
  pos.dupl <- unique(t[duplicated(t[,c('V4', 'V5')]), c('V4', 'V5'), drop=F])
  # print(pos.dupl)
  for(idupl in 1:nrow(pos.dupl)){
    if(idupl > nrow(pos.dupl)) break
    t <- t[!((t$V4 == pos.dupl[idupl,1]) & ((t$V5 == pos.dupl[idupl,2]))),]
  }
  
  
  t = removeCompleteOverlaps(t)
  
  return(t)
  
}


removeCentromeres <- function(t, base.len, cen.start.query, cen.end.query,
                              cen.start.base, cen.end.base, i.chr){
  
  t.base <- getBase(t, base.len)
  t <- t[(t.base$V2 <= cen.start.query[i.chr]) | (t.base$V3 >= cen.end.query[i.chr]),]
  
  t.base <- getBase(t, base.len)
  t <- t[(t.base$V4 <= cen.start.base[i.chr]) | (t.base$V5 >= cen.end.base[i.chr]),]
  return(t)
}


getIncorrectPlace <- function(t, base.len, correct.size=10^5){
  # Majority at the correct places
  correct.place = rep(0, nrow(t))
  idx.length = order(-t$V7)
  correct.place[t$V7 > correct.size] = 1
  
  t.base = getBase(t, base.len)
  sum(correct.place == 1)
  if(sum(correct.place == 0) == 0) return(c())
  for(i.place in (sum(correct.place == 1) + 1):nrow(t)){
    
    correct.idx = which(correct.place == 1)
    idx.place = idx.length[i.place]
    # if(idx.place == 6) break
    # ----
    # Get up and down indexes
    idx.up = c(correct.idx[correct.idx > idx.place])
    idx.down = c(correct.idx[correct.idx < idx.place])
    if(length(idx.up) == 0) {
      idx.up = idx.place
      up.V4 = Inf
    } else {
      idx.up = min(idx.up)
      up.V4 = t.base$V4[idx.up]
    }
    if(length(idx.down) == 0) {
      idx.down = idx.place
      down.V5 = 0
    } else {
      idx.down = max(idx.down)
      down.V5 = t.base$V5[idx.down]
    }
    # ----

    mean.place = (t.base$V4[idx.place] + t.base$V5[idx.place]) / 2
    
    if((down.V5 <= mean.place) & (mean.place <= up.V4)){
      correct.place[idx.place] = 1
    }
    if((up.V4 <= mean.place) & (mean.place <= down.V5)){
      correct.place[idx.place] = 1
    }
    
    if(t.base$dir[idx.place] == 0){
      if(t.base$dir[idx.down] == 1){
        if((mean.place <= up.V4)) correct.place[idx.place] = 1
      }
      if(t.base$dir[idx.up] == 1){
        if(down.V5 <= mean.place) correct.place[idx.place] = 1
      }
    }
    
    if(t.base$dir[idx.place] == 1){
      if(t.base$dir[idx.down] == 0){
        if((up.V4 <= mean.place)) correct.place[idx.place] = 1
      }
      if(t.base$dir[idx.up] == 0){
        if(mean.place <= down.V5) correct.place[idx.place] = 1
      }
    }
    
    
    # if(t.base$dir[idx.place] == 0){
    #   if((down.V5 <= mean.place) & (mean.place <= up.V4)){
    #     correct.place[idx.place] = 1
    #   }
    # } else { # dir of place == 1
    #   if((down.V5 == 1) && (t.base$dir[idx.up] == 1)){
    #     if((down.V5 >= mean.place) & (mean.place >= up.V4)){
    #       correct.place[idx.place] = 1
    #     }
    #   } else if((down.V5 == 1) && (t.base$dir[idx.up] == 0)){
    #     if((down.V5 >= mean.place) & (mean.place <= up.V4)){
    #       correct.place[idx.place] = 1
    #     }
    #   } else if((t.base$dir[idx.down] == 0) && (t.base$dir[idx.up] == 1)){
    #     if((down.V5 <= mean.place) & (mean.place >= up.V4)){
    #       correct.place[idx.place] = 1
    #     }
    #   } else {
    #     if((up.V4 <= mean.place) & (mean.place <= down.V5)){
    #       correct.place[idx.place] = 1
    #     }
    #   }
    # }
    
  }
  incorrect.idx = which(correct.place == 0)
  return(incorrect.idx)
}


getFragmentedRegions <- function(t, frag.thresh = 10000){
  # remove fragmented regions
  
  frag.len = (t$V7 < frag.thresh) * 1
  frag.region = frag.len[1:(nrow(t) - 1)] + frag.len[2:nrow(t)]
  frag.index = which(frag.region == 2)
  frag.index = sort(unique(c(frag.index, frag.index + 1)))
  frag.regions <- c()
  frag.beg = 1
  if(length(frag.index) == 0)  return(NULL)
  for(i in 2:length(frag.index)){
    if(frag.index[i] == (frag.index[i-1] + 1)) {
      frag.end = i
    } else {
      frag.regions <- rbind(frag.regions, c(frag.index[frag.beg], frag.index[i-1]))
      frag.beg = i
    }
  }
  frag.regions <- rbind(frag.regions, c(frag.index[frag.beg], frag.index[length(frag.index)]))
  return(frag.regions)
}

getFragmentedRows <- function(t){
  frag.regions <- getFragmentedRegions(t)
  if(is.null(frag.regions)) return(NULL)
  frag.idx <- c()
  for(i in 1:nrow(frag.regions)){
    frag.idx = c(frag.idx, frag.regions[i,1]:frag.regions[i,2])
  }
  return(frag.idx)
}

removeOverlapdAfterBlast <- function(t){
  unique_names = unique(t$V1)
  t_new = c()
  for(name in unique_names){
    #print(name)
    t_tmp = t[t$V1 == name,]
    if(nrow(t_tmp) > 1){
      t_tmp = t_tmp[order(-t_tmp$V6),]
      cover = t_tmp$V3 - t_tmp$V2 + 1
      t_tmp = t_tmp[order(-cover),]
      cover = t_tmp$V3 - t_tmp$V2 + 1
      
      pos = rep(0, 5000)
      
      idx_remove = c()
      pos[t_tmp$V2[1]:t_tmp$V3[1]] = 1
      for(irow in 2:nrow(t_tmp)){
        if(sum(pos[t_tmp$V2[irow]:t_tmp$V3[irow]]) > 0.85 * cover[irow]){
          idx_remove = c(idx_remove, irow)
        } else {
          pos[t_tmp$V2[irow]:t_tmp$V3[irow]] = 1
        }
      }
      # if(length(idx_remove) > 0) t_tmp = t_tmp[-idx_remove,]
      idx_remain = setdiff(1:nrow(t_tmp), idx_remove)
      for(irow in idx_remove){
        idx = ((t_tmp$V2[irow] %in% t_tmp$V2[idx_remain]) & 
                 (t_tmp$V3[irow] %in% t_tmp$V3[idx_remain]) &
                 (t_tmp$V6[irow] %in% t_tmp$V6[idx_remain]))
        if(sum(idx) > 0) idx_remain = c(idx_remain, irow)
      }
      t_tmp = t_tmp[idx_remain,]
    }
    t_new = rbind(t_new, t_tmp)
  }
  return(t_new)
}

getClosestRowQuery <- function(t, base.len, poses){
  close.beg = c()
  t.base = getBase(t, base.len)
  for(pos in poses){
    d = abs(t.base[,c('V2', 'V3')] - pos)
    if(sum(d == 1) > 1) stop('something is wrong')
    if(min(d) != 1){
      close.beg = c(close.beg, 0)
    } else {
      close.beg = c(close.beg, which(rowSums(d == 1) == 1))
    }
    print(min(d))
  }
  return(close.beg)
}

getClosestRowBase <- function(t, base.len, poses){
  close.beg = c()
  t.base = getBase(t, base.len)
  for(pos in poses){
    d = abs(t.base[,c('V4', 'V5')] - pos)
    if(sum(d == 1) > 1) stop('something is wrong')
    if(min(d) != 1){
      close.beg = c(close.beg, 0)
    } else {
      close.beg = c(close.beg, which(rowSums(d == 1) == 1))
    }
  }
  return(close.beg)
}

getSeqByName <- function(f, name){
  if(name %in% names(f)) return(toupper(paste0(f[[name]], sep = '', collapse = '')))
  return(NULL)
}


getMobCovered <- function(b, sim.cutoff = 0.85, f.save = NULL, echo=F, to.get.length=T){

  if(to.get.length){
    if(echo) print('Get Lengths')
    if(echo) print(nrow(b))
    b$len1 = as.numeric(sapply(strsplit(b$V1, "_"), "[", 3))
    b$len2 = as.numeric(sapply(strsplit(b$V8, "_"), "[", 3))    
  }
  b = b[b$V1 != b$V8,]  
  b$comb = paste(b$V1, b$V8, sep = '|')
  
  
  if(!is.null(f.save)) saveRDS(b, f.save)
  
  idx.change = b$V5 < b$V4
  tmp = b$V4[idx.change]
  b$V4[idx.change] = b$V5[idx.change]
  b$V5[idx.change] = tmp
  
  
  # Remove singletons
  if(echo) print('Count combinations')
  comb.counts = table(b$comb)
  comb.counts = comb.counts[comb.counts >= 2]
  b.single = b[!(b$comb %in% names(comb.counts)),]
  b.fragments = b[b$comb %in% names(comb.counts),]
  
  
  # Sure combinations
  b.single$coverage1 = b.single$V3 - b.single$V2 + 1
  b.single$coverage2 = b.single$V5 - b.single$V4 + 1
  # b.sure = b.single[(b.single$coverage1/b.single$len1 >= sim.cutoff) | 
  #                     (b.single$coverage2/b.single$len2 >= sim.cutoff), ]
  
  
  # -- Sequence 1
  # Sorting
  if(echo) print('Fragments1')
  b.fragments = b.fragments[order(b.fragments$V3),]
  b.fragments = b.fragments[order(b.fragments$V2),]
  b.fragments = b.fragments[order(b.fragments$comb),]
  
  
  b.fragments$idx.same = c(F, b.fragments$comb[2:nrow(b.fragments)] == b.fragments$comb[1:(nrow(b.fragments)-1)])
  b.fragments$prev.end = c(0, b.fragments$V3[1:(nrow(b.fragments)-1)])
  
  idx.change = b.fragments$idx.same & (b.fragments$prev.end > b.fragments$V3)
  b.fragments$V3[idx.change] = b.fragments$prev.end[idx.change]
  
  b.fragments$len.aln = b.fragments$V3 - b.fragments$V2 + 1
  b.fragments$gap = c(b.fragments$V2[2:nrow(b.fragments)] - b.fragments$V3[1:(nrow(b.fragments)-1)] - 1, 0)
  b.fragments$gap[b.fragments$gap > 0] = 0 
  idx.comb.diff = c(b.fragments$comb[2:nrow(b.fragments)] != b.fragments$comb[1:(nrow(b.fragments)-1)], T)
  b.fragments$diff = idx.comb.diff
  b.fragments$gap[b.fragments$diff] = 0
  
  coverage1 = tapply(b.fragments$gap, b.fragments$comb, sum) + tapply(b.fragments$len.aln, b.fragments$comb, sum)
  
  b.fragments$coverage1 = coverage1[b.fragments$comb]
  
  
  # -- Sequence 2
  # Sorting
  if(echo) print('Fragments2')
  idx.change = b.fragments$V5 < b.fragments$V4
  tmp = b.fragments$V4[idx.change]
  b.fragments$V4[idx.change] = b.fragments$V5[idx.change]
  b.fragments$V5[idx.change] = tmp
  
  b.fragments = b.fragments[order(b.fragments$V5),]
  b.fragments = b.fragments[order(b.fragments$V4),]
  b.fragments = b.fragments[order(b.fragments$comb),]
  
  
  b.fragments$idx.same = c(F, b.fragments$comb[2:nrow(b.fragments)] == b.fragments$comb[1:(nrow(b.fragments)-1)])
  b.fragments$prev.end = c(0, b.fragments$V5[1:(nrow(b.fragments)-1)])
  
  idx.change = b.fragments$idx.same & (b.fragments$prev.end > b.fragments$V5)
  b.fragments$V5[idx.change] = b.fragments$prev.end[idx.change]
  
  b.fragments$len.aln = b.fragments$V5 - b.fragments$V4 + 1
  b.fragments$gap = c(b.fragments$V4[2:nrow(b.fragments)] - b.fragments$V5[1:(nrow(b.fragments)-1)] - 1, 0)
  b.fragments$gap[b.fragments$gap > 0] = 0 
  idx.comb.diff = c(b.fragments$comb[2:nrow(b.fragments)] != b.fragments$comb[1:(nrow(b.fragments)-1)], T)
  b.fragments$diff = idx.comb.diff
  b.fragments$gap[b.fragments$diff] = 0
  
  coverage2 = tapply(b.fragments$gap, b.fragments$comb, sum) + tapply(b.fragments$len.aln, b.fragments$comb, sum)
  
  b.fragments$coverage2 = coverage2[b.fragments$comb]
  
  # -- All together
  
  if(echo) print('All together')
  b.uni = rbind(b.single[, c('V1', 'V8', 'comb', 'len1', 'len2', 'coverage1', 'coverage2')],
                b.fragments[!duplicated(b.fragments$comb),
                            c('V1', 'V8', 'comb', 'len1', 'len2', 'coverage1', 'coverage2')])
  
  b.uni.cover = b.uni[(b.uni$coverage1/b.uni$len1 >= sim.cutoff) | 
                        (b.uni$coverage2/b.uni$len2 >= sim.cutoff), ]
  
  b.uni.cover$p1 = b.uni.cover$coverage1/b.uni.cover$len1
  b.uni.cover$p2 = b.uni.cover$coverage2/b.uni.cover$len2
  
  
  
  idx = b.uni.cover$p2 > b.uni.cover$p1
  tmp = b.uni.cover[idx, c('V1', 'len1', 'coverage1', 'p1')]
  b.uni.cover[idx, c('V1', 'len1', 'coverage1', 'p1')] = b.uni.cover[idx, c('V8', 'len2', 'coverage2', 'p2')]
  b.uni.cover[idx, c('V8', 'len2', 'coverage2', 'p2')] = tmp
  
  return(b.uni.cover)
}



getMobCovered1 <- function(b, sim.cutoff = 0.85, f.save = NULL, echo=F, to.get.length=T){
  
  if(to.get.length){
    if(echo) print('Get Lengths')
    b = b[b$V1 != b$V8,]  
    if(echo) print(nrow(b))
    b$len1 = as.numeric(sapply(strsplit(b$V1, "_"), "[", 3))
    b$len2 = as.numeric(sapply(strsplit(b$V8, "_"), "[", 3))  
  }
  b$comb = paste(b$V1, b$V8, sep = '|')
  
  
  
  if(!is.null(f.save)) saveRDS(b, f.save)
  
  idx.change = b$V5 < b$V4
  tmp = b$V4[idx.change]
  b$V4[idx.change] = b$V5[idx.change]
  b$V5[idx.change] = tmp
  
  
  # Remove singletons
  if(echo) print('Count combinations')
  comb.counts = table(b$comb)
  comb.counts = comb.counts[comb.counts >= 2]
  b.single = b[!(b$comb %in% names(comb.counts)),]
  b.fragments = b[b$comb %in% names(comb.counts),]
  
  
  # Sure combinations
  b.single$coverage1 = b.single$V3 - b.single$V2 + 1
  # b.sure = b.single[(b.single$coverage1/b.single$len1 >= sim.cutoff) | 
  #                     (b.single$coverage2/b.single$len2 >= sim.cutoff), ]
  
  
  # -- Sequence 1
  # Sorting
  if(echo) print('Fragments1')
  b.fragments = b.fragments[order(b.fragments$V3),]
  b.fragments = b.fragments[order(b.fragments$V2),]
  b.fragments = b.fragments[order(b.fragments$comb),]
  
  
  b.fragments$idx.same = c(F, b.fragments$comb[2:nrow(b.fragments)] == b.fragments$comb[1:(nrow(b.fragments)-1)])
  b.fragments$prev.end = c(0, b.fragments$V3[1:(nrow(b.fragments)-1)])
  
  idx.change = b.fragments$idx.same & (b.fragments$prev.end > b.fragments$V3)
  b.fragments$V3[idx.change] = b.fragments$prev.end[idx.change]
  
  b.fragments$len.aln = b.fragments$V3 - b.fragments$V2 + 1
  b.fragments$gap = c(b.fragments$V2[2:nrow(b.fragments)] - b.fragments$V3[1:(nrow(b.fragments)-1)] - 1, 0)
  b.fragments$gap[b.fragments$gap > 0] = 0 
  idx.comb.diff = c(b.fragments$comb[2:nrow(b.fragments)] != b.fragments$comb[1:(nrow(b.fragments)-1)], T)
  b.fragments$diff = idx.comb.diff
  b.fragments$gap[b.fragments$diff] = 0
  
  coverage1 = tapply(b.fragments$gap, b.fragments$comb, sum) + tapply(b.fragments$len.aln, b.fragments$comb, sum)
  
  b.fragments$coverage1 = coverage1[b.fragments$comb]
  
  
  # -- All together
  
  if(echo) print('All together')
  b.uni = rbind(b.single[, c('V1', 'V8', 'comb', 'len1', 'coverage1')],
                b.fragments[!duplicated(b.fragments$comb),
                            c('V1', 'V8', 'comb', 'len1', 'coverage1')])
  
  b.uni.cover = b.uni[(b.uni$coverage1/b.uni$len1 >= sim.cutoff), ]
  
  return(b.uni.cover)
}

cleanJumps <- function(t, base.len, echo=F){
  f.change = T
  while(f.change){
    f.change = F
    t = t[order(t$V2),]
    rownames(t) <- NULL
    t.base = getBase(t, base.len)
    t.base = t.base[order(t.base$V4),]
    
    idx.names = c(0, as.numeric(rownames(t.base)), (nrow(t.base) +1))
    idx.dist = abs(idx.names[-1] - idx.names[-length(idx.names)]) - 1
    
    d = cbind(idx.dist[-length(idx.dist)], idx.dist[-1])
    idx.bad = which((rowSums(d == 0) == 0))
    # print(d[idx.bad,])
    # print(sort(idx.names[idx.bad] + 1))
    # print(length(idx.bad))
    if(length(idx.bad) == 0) next
    d.sum = rowSums(d)
    idx.bad = idx.names[idx.bad[d.sum[idx.bad] == max(d.sum[idx.bad])] + 1]
    # print(nrow(t))
    # if(idx.bad == 713) stop()
    if(echo) print(idx.bad)
    if(echo) print(t[idx.bad,1:7])
    t = t[-idx.bad,]
    f.change = T
  }
  rownames(t) <- NULL
  return(t)
}

cleanJumps2 <- function(t, base.len, echo=F, min.len = 20000){
  f.change = T
  
  while(f.change){
    f.change = F
    t = t[order(t$V2),]
    
    rownames(t) <- NULL
    t.base = getBase(t, base.len)
    t.base = t.base[order(t.base$V4),]
    idx.ok = which(t.base$V7 > min.len)
    
    idx.names = c(0, as.numeric(rownames(t.base)), (nrow(t.base) +1))
    idx.dist = abs(idx.names[-1] - idx.names[-length(idx.names)]) - 1
    
    d = cbind(idx.dist[-length(idx.dist)], idx.dist[-1])
    idx.bad = which((rowSums(d == 0) == 0))
    idx.bad = setdiff(idx.bad, idx.ok)
    # print(d[idx.bad,])
    # print(sort(idx.names[idx.bad] + 1))
    # print(length(idx.bad))
    if(length(idx.bad) == 0) next
    d.sum = rowSums(d)
    idx.bad = idx.names[idx.bad[d.sum[idx.bad] == max(d.sum[idx.bad])] + 1]
    # print(nrow(t))
    # if(idx.bad == 713) stop()
    if(echo) print(idx.bad)
    if(echo) print(t[idx.bad,1:7])
    t = t[-idx.bad,]
    f.change = T
  }
  rownames(t) <- NULL
  return(t)
}

