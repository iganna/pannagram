#' Previous and Next values in the vector
#'
#' If some elements are equal to 0, then this fucntion will find the previous and next non-zero element,
#' if this 0 is in one block with the flanking non-zero elements
#'
#' @param x.acc A numeric vector with positions.
#' @return A list containing two elements: 'prev' and 'next', each is a vector of calculated values.
getPrevNext <- function(x.acc){
  ### ---- Find prev and next ----
  n = length(x.acc)
  for(i.tmp in 1:2){
    v.acc = x.acc
    if(i.tmp == 2) {
      v.acc = rev(v.acc)
    }
    v.rank = rank(v.acc)
    v.rank[v.acc == 0] = 0
    
    v.zero.beg = which((v.acc[-1] == 0) & (v.acc[-n] != 0)) + 1
    v.zero.end = which((v.acc[-1] != 0) & (v.acc[-n] == 0))
    if(v.acc[1] == 0) v.zero.end = v.zero.end[-1]
    if(v.acc[n] == 0) v.zero.end = c(v.zero.end, n)
    
    # ..... WITHIN ONE STRATCH BLOCK .....
    idx = which(abs(v.rank[v.zero.beg-1] - v.rank[v.zero.end+1]) != 1)
    v.zero.beg = v.zero.beg[-idx]
    v.zero.end = v.zero.end[-idx]
    # .....
    
    tmp = rep(0, n)
    tmp[v.zero.beg] = v.acc[v.zero.beg-1]
    tmp[v.zero.end] = -v.acc[v.zero.beg-1]
    tmp[v.zero.end[v.zero.beg == v.zero.end]] = 0
    tmp = cumsum(tmp)
    tmp[v.zero.end] = v.acc[v.zero.beg-1]
    
    if(i.tmp == 1){
      v.prev = x.acc + tmp
    } else {
      tmp = rev(tmp)
      v.next = x.acc + tmp
    }
  }
  
  return(list(v.prev = v.prev, 
              v.next = v.next))
}
