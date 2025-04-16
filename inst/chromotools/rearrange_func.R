#' Select Best Correspondence Regions Among Chromosomes
#'
#' Finds the best matching regions between reference and query chromosomes based on a BLAST-scoring matrix.
#'
#' @param pos Numeric matrix with scores, rows correspond to queries, columns to positions in reference.
#' @param i.chr.ref.len Integer, length of the reference chromosome.
#' @param i.chr.ref Integer or string, identifier of the ID of the reference chromosome.
#' @param i.chr.corresp Vector of identifiers for query chromosomes.
#' @param min.overlap.fragment Numeric, minimal fraction of reference length for a valid region.
#'
# @export
findBestChromosome <- function(pos, 
                       i.chr.ref.len, 
                       i.chr.ref, 
                       i.chr.corresp, 
                       min.len, to.extend = T) {
  
  # Find maxumum by chromosomes
  pos.max <- pos[1, ]
  pos.i <- as.integer(pos.max != 0)
  for (i in 2:nrow(pos)) {
    idx.i <- which(pos[i, ] > pos.max)
    pos.max[idx.i] <- pos[i, idx.i]
    pos.i[idx.i] <- i
  }
  rm('pos')
  
  pos.idx = 1:i.chr.ref.len
  idx.remain = pos.i != 0
  pos.idx = pos.idx[idx.remain]
  pos.i = pos.i[idx.remain]
  
  # Ger regions
  df.all = c()
  for(i in 1:length(i.chr.corresp)){
    df = findOnes((pos.i == i) * 1)
    if(nrow(df) == 0) next
    df$len = df$end - df$beg + 1
    
    df$i.ref = i.chr.ref
    df$i.acc =  i.chr.corresp[i]
    
    df.all = rbind(df.all, df)
  }
  df.all = df.all[order(df.all$beg),]
  
  # Gradually remove the smallest
  while(min(df.all$len) < min.len){
    idx = which.min(df.all$len)[1]
    df.all = df.all[-idx,,drop=F]
    idx.merge = which(diff(df.all$i.acc) == 0)
    if(length(idx.merge) == 0) next
    idx.merge = idx.merge[1]
    df.all$end[idx.merge] = df.all$end[idx.merge + 1]
    df.all$len[idx.merge] = df.all$end[idx.merge] - df.all$beg[idx.merge] + 1
    df.all = df.all[-(idx.merge + 1),,drop=F]
    if(nrow(df.all) == 0) stop("no correspondence is left")
    # print(df.all)
  }
  df.all$beg = pos.idx[df.all$beg]
  df.all$end = pos.idx[df.all$end]
  df.all$len = df.all$end - df.all$beg + 1
  
  
  
  # print(df.all)
  if(to.extend){
    
    df.all = df.all[order(df.all$beg),]
    
    for (i in 1:(nrow(df.all)-1)) {
      gap <- df.all$beg[i+1] - df.all$end[i] - 1
      shift <- floor(gap / 2)
      
      df.all$end[i]   <- df.all$end[i] + shift
      df.all$beg[i+1] <- df.all$end[i] + 1
    }
    
    df.all$beg[1] = 1
    df.all$end[nrow(df.all)] = i.chr.ref.len
    df.all$len = df.all$end - df.all$beg + 1
  }

  return(df.all)
}
