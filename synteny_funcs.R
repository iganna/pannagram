#' Clean All Overlaps in a Data Frame
#'
#' This function cleans both big and small overlaps in the provided data frame.
#' It applies the cleaning functions to query and base sequences.
#'
#' @param x.df A data frame containing the alignment with overlaps
#' @param rm.threshold A numeric value between 0 and 1 indicating the threshold for 
#'                     determining significant overlaps. Default is 0.5.
#'
#' @return A data frame cleaned of both big and small overlaps.
#'
cleanOverlaps <- function(x.df, rm.threshold = 0.5){
  
  # Remove short overlaps: twice, because from "both sides"
  for(i.tmp in 1:2){
    x.df = cleanBigOverlaps(x.df, rm.threshold = rm.threshold)
    x.df = cutSmallOverlaps(x.df)
  }
  
  for(i.tmp in 1:2){
    x.df = cleanBigOverlapsQuery(x.df, rm.threshold = rm.threshold)
    x.df = cutSmallOverlapsQuery(x.df)
  }
  
  return(x.df)
  
}

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
defineOverlappsQuery <- function(x.df){
  
  # if not sorted - SORT!
  if(is.unsorted(x.df$V2)){
    # Sort by the position in the reference
    x.df = x.df[order(-x.df$V7),]
    x.df = x.df[order(x.df$V2),]
  }
  
  x.df$rm.len = 0
  idx.overlap = which(x.df$V2[-1] <= x.df$V3[-nrow(x.df)])
  
  # x.df[c(irow, irow+1),]
  
  for(irow in idx.overlap){
    # We cut either irow or [irow+1]
    
    # Which row to cut
    icut = ifelse(x.df$V7[irow] > x.df$V7[irow + 1], irow + 1, irow)
    ibig = ifelse(x.df$V7[irow] > x.df$V7[irow + 1], irow, irow + 1)
    
    # Find left(value = 1) or right(value = -1) tails 
    # of irow(i.e., tail.irow) or [irow+1](i.e., tail.next), 
    # which are involved in the overlap
    tail.iсut = ifelse((x.df$V2[icut] >= x.df$V2[ibig]) & 
                         (x.df$V2[icut] <= x.df$V3[ibig]), 
                       1, -1)
    
    # how much to cut
    ncut = length(intersect(x.df$V2[irow]:x.df$V3[irow],
                            x.df$V2[irow+1]:x.df$V3[irow+1]))
    # Remember the cut
    x.df$rm.len[icut] = ncut * tail.iсut
  }
  
  # if(!sort.flaf){
  #   x.df[order(x.df[,new.col.name]),]
  #   x.df = x.df[,colnames(x.df) != new.col.name]
  # }
  
  
  return(x.df)
}


#' Clean Big Overlaps in Query Sequences in the alignment
#'
#' This function iteratively cleans rows in a data frame based on overlap criteria.
#' It uses the `defineOverlappsQuery` function to determine overlaps 
#' and removes rows based on the `rm.threshold`.
#' If the ovelaps is higher than 50% of length - remove.
#'
#' @param x.df A data frame containing overlap or alignment information.
#' @param rm.threshold A numeric value between 0 and 1 indicating the threshold for removing rows based on overlap criteria. Default is 0.5.
#'
#' @return A data frame after cleaning rows based on overlap criteria.
#'
cleanBigOverlapsQuery <- function(x.df, rm.threshold = 0.5){
  
  # Check if rm.threshold is between 0 and 1
  if (rm.threshold < 0 || rm.threshold > 1) {
    stop("Error: rm.threshold must be between 0 and 1.")
  }
  
  while(T){
    n.row = nrow(x.df)
    x.df = defineOverlappsQuery(x.df)
    x.df =  x.df[((abs(x.df$rm.len) / x.df$V7) <= rm.threshold) | (x.df$rm.len == 0),]
    if(n.row == nrow(x.df)) break
    pokaz(n.row)
  }
  
  return(x.df)
}

#' Remove Small Overlaps in Query sequences
#'
#' The function identifies and processes peviously defined overlaps.
#'
#' @param x.sk A data frame with the alignment
#'
#' @return A modified data frame without overlaps.
#' 
cutSmallOverlapsQuery <- function(x.df){
  
  # if(!isSorted(x.df$p.beg)) pokazAttention('2!!')
  
  idx.cut = which(x.df$rm.len != 0)
  for(irow in idx.cut){
    
    # print(info(x.df[irow,]))
    # if(x.df$dir[irow] == 1) stop('dir')
    
    adjustment <- x.df$rm.len[irow]  # sing says about "begin" or "end"
    # adjustment.dir <- ifelse(x.df$dir[irow] == 0, adjustment, (-1) * adjustment)
    n.symbols = abs(adjustment)
    
    if(adjustment > 0){
      # stop('3')
      # From the befinning!
      
      # Adjust strings in the alignment
      seq = seq2nt(x.df$V8[irow])
      aln.adjust = which(seq != '-')[n.symbols]
      x.df$V8[irow] = substr(x.df$V8[irow], (aln.adjust+1), nchar(x.df$V8[irow]))
      x.df$V9[irow] = substr(x.df$V9[irow], (aln.adjust+1), nchar(x.df$V9[irow]))

      
      # Adjust positions - query
      s.q.cut = seq2nt(x.df$V8[irow])
      adjustment.q = sum(s.q.cut != '-') - 1
      x.df$V2[irow] = x.df$V3[irow] - adjustment.q
      
      # Adjust positions - base
      s.b.cut = seq2nt(x.df$V9[irow])
      adjustment.b = sum(s.b.cut != '-') - 1
      x.df$V4[irow] = x.df$V5[irow] - adjustment.b * sign(0.5 - x.df$dir[irow])
      
      
    }
    if(adjustment < 0){
      # stop('4')
      # From the ending!
      
      # Adjust strings in the alignment
      seq = seq2nt(x.df$V8[irow])
      aln.adjust = tail(which(seq != '-'), n=n.symbols)[1]
      x.df$V8[irow] = substr(x.df$V8[irow], 1, (aln.adjust-1))
      x.df$V9[irow] = substr(x.df$V9[irow], 1, (aln.adjust-1))
      
      
      
      # Adjust positions
      
      s.q.cut = seq2nt(x.df$V8[irow])
      adjustment.q = sum(s.q.cut != '-') - 1
      x.df$V3[irow] = x.df$V2[irow] + adjustment.q
      
      
      # if(x.df$dir[irow] == 1) stop('dir')
      # Adjust positions - base
      s.b.cut = seq2nt(x.df$V9[irow])
      adjustment.b = sum(s.b.cut != '-') - 1
      x.df$V5[irow] = x.df$V4[irow] + adjustment.b * sign(0.5 - x.df$dir[irow])
    }
    
    # print(info(x.df[irow,]))
    
    x.df$V7[irow] = nchar(x.df$V8[irow])  # doesn't matter 8 or 9 here
    
    # print(x.df[irow, -c(8,9)])
    # if(x.df$V1[irow] == 'acc_10001|chr_1|part_157|780001') stop()
  }
  
  
  # Undate begin-eng positions
  x.df$p.beg <- ifelse(x.df$V4 < x.df$V5, x.df$V4, x.df$V5)
  x.df$p.end <- ifelse(x.df$V4 < x.df$V5, x.df$V5, x.df$V4)
  
  return(x.df)
}



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


#' Clean Big Overlaps in Base Sequences in the alignment
#'
#' This function iteratively cleans rows in a data frame based on overlap criteria.
#' It uses the `defineOverlappsQuery` function to determine overlaps 
#' and removes rows based on the `rm.threshold`.
#' If the ovelaps is higher than 50% of length - remove.
#'
#' @param x.df A data frame containing overlap or alignment information.
#' @param rm.threshold A numeric value between 0 and 1 indicating the threshold for removing rows based on overlap criteria. Default is 0.5.
#'
#' @return A data frame after cleaning rows based on overlap criteria.
#'
cleanBigOverlaps <- function(x.df, rm.threshold = 0.5){
  
  # Check if rm.threshold is between 0 and 1
  if (rm.threshold < 0 || rm.threshold > 1) {
    stop("Error: rm.threshold must be between 0 and 1.")
  }
  
  while(T){
    n.row = nrow(x.df)
    x.df = defineOverlapps(x.df)
    x.df =  x.df[((abs(x.df$rm.len) / x.df$V7) <= rm.threshold) | (x.df$rm.len == 0),]
    if(n.row == nrow(x.df)) break
    pokaz(n.row)
  }
  
  return(x.df)
}

#' Remove Small Overlaps in Base sequences
#'
#' The function identifies and processes peviously defined overlaps.
#'
#' @param x.sk A data frame with the alignment
#'
#' @return A modified data frame without overlaps.
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


#' Check Correspondence to Genome
#'
#' This function examines the alignment data and compares the end positions of synteny blocks
#' to the original genomes. It ensures that all positions remain accurate and were not distorted
#' during any insinuations in the alignment process.
#'
#' @param t A data frame containing alignment information.
#' @param base.fas.fw A character vector representing the forward base genome sequence.
#' @param base.fas.bw A character vector representing the backward base genome sequence.
#' @param query.fas A character vector representing the query genome sequence.
#' @param k An integer determining the length of the end positions to check. Default is 10.
#'
#' @return This function does not return a value. Instead, it stops with an error message 
#'         if any discrepancies are found in the end positions of the alignment.
#'
checkCorrespToGenome <- function(x, base.fas.fw, base.fas.bw, query.fas, k = 10){
  for(irow in 1:nrow(x)) {
    s1 = toupper(remainLastN(x[irow, 'V8'], k))
    s1 = gsub("\\-","",s1)
    
    if(nchar(s1) == 0){
      s2 = s1
    } else {
      s2 = toupper(paste0(query.fas[(-(nchar(s1)-1):0) + x[irow, 'V3']], collapse = ''))
    }
    
    
    if(s1 != s2 ){
      pokaz('Row', irow)
      pokaz(s1)
      pokaz(s2)
      stop('Checking the query not passed')
    }
    
    s1 = toupper(remainLastN(x[irow, 'V9'], k))
    s1 = gsub("\\-","",s1)
    if(x$dir[irow] == 0) {
      base.fas = base.fas.fw
    } else {
      base.fas = base.fas.bw
    }
    
    if(nchar(s1) == 0){
      s2 = s1
    } else {
      s2 = toupper(paste0(base.fas[(-(nchar(s1)-1):0) + x[irow, 'V5']], collapse = ''))
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

#' Reconstruct Path from Predecessor Vector
#'
#' This function reconstructs a path based on a given predecessor vector.
#' The path is determined by traversing the predecessor vector from end to start.
#'
#' @param predecessor A numeric vector containing the predecessor index for each position.
#'
#' @return A numeric vector indicating the reconstructed path.
#'
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



#' Identify the Optimal Alignment Path in the Dencrease-Pos-Query (Down) Direction from the Synteny Block Begin
#'
#' This function identifies the optimal alignment path when base and query directions are the same (Plus).
#' 
#' It first flips the base and query positions and then uses the `pathUpPlus` function 
#' to determine the best alignment path.
#'
#' @param x.tmp A data frame containing alignment information.
#' 
#' @return An integer vector indicating the rows of x.tmp, which form the best alignment path.
#'
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

#' Identify the Optimal Alignment Path in the Dencrease-Pos-Query (Down) Direction from the Synteny Block Begin
#'
#' This function identifies the optimal alignment path when base and query directions are opposite (Minus).
#' It first flips the query positions and then uses the `pathUpPlus` function 
#' to determine the best alignment path.
#'
#' @param x.tmp A data frame containing alignment information.
#'
#' @return An integer vector indicating the rows of x.tmp, which form the best alignment path.
#'
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


#' Identify the Optimal Alignment Path in the Increase-Pos-Query (Up) Direction from the Synteny Block End
#'
#' This function identifies the optimal alignment path when base and query directions are opposite (Minus). 
#' It first inverts the base positions and then utilizes the `pathUpPlus` function to determine the 
#' best alignment path.
#'  
#' @param x.tmp A data frame containing alignment information.
#'
#' @return An integer vector indicating the rows of x.tmp, which form the best alignment path.
#'
pathUpMinus <- function(x.tmp){
  
  # Flip base
  b.max = max(c(x.tmp$V4, x.tmp$V5)) + 1
  x.tmp$V4 = b.max - x.tmp$V4
  x.tmp$V5 = b.max - x.tmp$V5
  
  # DON'T SORT
  # PLEASE NO SORTING AND NO CHANGES OF V4-V5
  
  idx.visit = pathUpPlus(x.tmp)
  return(idx.visit)
}

#' Identify the Optimal Alignment Path in the Increase-Pos-Query (Up) Direction from the Synteny Block End
#'
#' This function identifies the optimal alignment path when base and query directions are the same (Plus).
#' 
#' @param x.tmp A data frame containing alignment information.
#'
#' @return An integer vector indicating the rows of x.tmp, which form the best alignment path.
#'
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



#' Get Correspondence of Query coordinates to Base Coordinates
#'
#' This function determines the correspondence of positions from the alignment.
#'
#' @param x A data frame or matrix containing alignment information.
#' @param base.len An integer representing the length of the base sequence.
#'
#' @return A numeric vector indicating the correspondence of positions between the sequences.
#' Positive values indicate direct correspondence, while negative values indicate reverse correspondence.
#'
getCorresp2BaseSign <- function(x, base.len){
  
  pos.corresp = rep(0, base.len)
  for(irow in 1:nrow(x)){
    
    aln.len = nchar(x$V8[irow])
    
    # Occupied positions in query
    positions.q = rep(0, aln.len)
    positions.q[seq2nt(x$V8[irow]) != '-'] = x$V2[irow]:x$V3[irow]
    
    # Occupied positions in base
    positions.b = rep(0, aln.len)
    if(x$dir[irow] == 0){
      positions.b[seq2nt(x$V9[irow]) != '-'] = t.base$V4[irow]:t.base$V5[irow]  
    } else {
      positions.b[seq2nt(x$V9[irow]) != '-'] = -(t.base$V5[irow]:t.base$V4[irow])
    }
    
    positions.q = positions.q[positions.b != 0]
    positions.b = positions.b[positions.b != 0]
    pos.corresp[abs(positions.b)] = positions.q * sign(positions.b)
    
  }
  return(pos.corresp)
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


