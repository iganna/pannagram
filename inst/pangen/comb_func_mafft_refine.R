
refineMafft <- function(mx, n.flank = 30){
  
  
  sim.cutoff = 0.2
  s.nt.fake = '!'
  
  
  # Create the mayrix with posisionts
  mx.pos = matrix(0, nrow = nrow(mx), ncol = ncol(mx), 
                  dimnames = list(rownames(mx), NULL))
  for(irow in 1:nrow(mx.pos)){
    pos.row = as.numeric(strsplit(rownames(mx)[irow], '\\|')[[1]][3:4])
    pos.row = pos.row[1]:pos.row[2] 
    pos.row = 1:length(pos.row)
    
    pos.filled = which(mx[irow,] != '-')
    pos.filled = pos.filled[-(1:n.flank)]
    pos.filled <- pos.filled[-((length(pos.filled) - n.flank + 1):length(pos.filled))]
    
    mx.pos[irow, pos.filled] = pos.row
  }
  
  # Remove flanking regions
  mx[mx.pos == 0] = '-'
  idx.col.remain = (colSums(mx != '-') != 0)
  mx = mx[,idx.col.remain]
  mx.pos = mx.pos[,idx.col.remain]
  
  # Save variables
  mx.init = mx
  mx.pos.init = mx.pos
  
  # Plot
  # p0 = msaplot(mx)
  # p0
  
  # --- --- --- --- --- --- --- ---
  # ---- Split ----
  
  # Final alignment
  mx.ok = c()
  mx.pos.ok = c()
  
  # To arrange a cycle
  mx.rest = mx.init
  mx.pos.rest = mx.pos.init
  
  for(n.round in 1:5){
    # pokaz('Disentangle Round', n.round)
    
    mx = mx.rest           # The matrix you work with and change values
    mx.pos = mx.pos.rest   # The matrix you work with and change values
    
    # Diversity by each position
    aln.len = ncol(mx)
    mx = toupper(mx)
    s.nts = c('A', 'C', 'G', 'T')
    pos.profile = matrix(0, nrow = 4, ncol = aln.len, dimnames = list(c(s.nts, NULL)))
    for(s.nt in s.nts){
      pos.profile[s.nt,] = colSums(mx == s.nt)
    }
    pos.variation = (colSums(pos.profile == 0) != 3) * 1
    
    # Define blocks, were the alignment non well
    blocks.all = c()
    for(i.seq in 1:nrow(mx)){
      s = mx[i.seq,]
      blocks = findOnes((s != '-') *1)
      blocks$len = blocks$end - blocks$beg + 1
      if(nrow(blocks) == 0) next
      
      # Estimate diversity within each block
      blocks$pi = 0
      blocks$acc = rownames(mx)[i.seq]
      for(irow in 1:nrow(blocks)){
        idx.block = blocks$beg[irow]:blocks$end[irow]
        blocks$pi[irow] = sum(pos.variation[idx.block]) / blocks$len[irow]
      }
      
      blocks.all = rbind(blocks.all, blocks)
    }
    if(nrow(blocks) == 0){
      # pokaz('Exit 1')
      mx.ok = cbind(mx.ok, mx)
      mx.pos.ok = cbind(mx.pos.ok, mx.pos)
      break
    } 
    
    blocks = blocks.all
    rm(blocks.all)
    blocks = blocks[blocks$pi>sim.cutoff,,drop = F]
    blocks = blocks[order(-blocks$pi),,drop = F]
    
    if(nrow(blocks) == 0){
      # pokaz('Exit 2')
      mx.ok = cbind(mx.ok, mx)
      mx.pos.ok = cbind(mx.pos.ok, mx.pos)
      break
    }
    
    blocks$remove = 0
    for(irow in 1:nrow(blocks)){
      idx.block = blocks$beg[irow]:blocks$end[irow]
      blocks$pi[irow] = sum(pos.variation[idx.block]) / blocks$len[irow]
      if(blocks$pi[irow] <= sim.cutoff){
        # stop()
        # pokaz(irow)
        next
      } 
      
      seq.tmp = mx[blocks$acc[irow], idx.block]  
      for(s.nt in s.nts){
        pos.profile[s.nt,idx.block] = pos.profile[s.nt,idx.block] - 1 * (seq.tmp == s.nt)
      }
      
      blocks$remove[irow] = 1
      pos.variation[idx.block] = (colSums(pos.profile[,idx.block,drop = F] == 0) != 3) * 1
    }
    
    
    # # Optional
    # for(irow in 1:nrow(blocks)){
    #   idx.block = blocks$beg[irow]:blocks$end[irow]
    #   loneliness = max(pos.profile[,idx.block,drop=F])
    #   if(loneliness == 1){
    #     # pokaz(irow)
    #     # stop()
    #     blocks$remove[irow] = 1
    #   } 
    # }
    
    
    
    # Remove blocks with very high diversity
    for(irow in 1:nrow(blocks)){
      if(blocks$remove[irow] == 0) next
      idx.block =  blocks$beg[irow]:blocks$end[irow]
      mx[blocks$acc[irow],idx.block] = s.nt.fake
    }
    
    
    # Well alignmed positions
    mask.ok = (mx != s.nt.fake)
    mx[!mask.ok] = '-'
    mx.pos[!mask.ok] = 0
    
    mx.ok = cbind(mx.ok, mx)
    mx.pos.ok = cbind(mx.pos.ok, mx.pos)
    
    # The rest for the next round
    
    mx.rest[mask.ok] = '-'
    mx.pos.rest[mask.ok] = 0
    
    cols.rest = colSums(mx.rest != '-') > 0
    mx.rest = mx.rest[, cols.rest, drop = F]
    mx.pos.rest = mx.pos.rest[, cols.rest, drop = F]
    
    if(ncol(mx.rest) == 0){
      # pokaz('Exit 3')
      break
    }
    
  }
  if(sum(rowSums(mx.pos.ok != 0) != rowSums(mx.pos.init != 0)) != 0) {
    pokazAttention('Not enough rounds') 
    return(list(mx = mx.init, pos = mx.pos.init))
  }
  if(sum(dim(mx.ok) != dim(mx.pos.ok)) != 0) stop('Sizes of matrices are different')
  
  # p2 = msaplot(mx.ok)
  # p2
  
  
  # ---- Ordering ----
  mx = mx.ok
  mx.pos = mx.pos.ok
  # 
  # # Remove total gaps
  # idx.remain = colSums(mx != '-') > 0
  # mx = mx[,idx.remain]
  # mx.pos = mx.pos[,idx.remain]
  # mx.len = ncol(mx)
  # 
  # 
  # # blocks of connetedness
  # b.beg = 1
  # b.end = c()
  # mx.beg = c()
  # mx.end = c()
  # 
  # val.beg = mx.pos[,b.beg]
  # val.end = mx.pos[,b.beg]
  # for(icol in 2:mx.len){
  #   val.nex = mx.pos[,icol]
  #   
  #   d = val.nex - val.end
  #   d = d[(val.end != 0) & (val.nex != 0)]
  #   if(length(d) == 0) d = 0
  #   
  #   if(sum(d != 1) == 0){
  #     val.beg[val.beg == 0] = val.nex[val.beg == 0]
  #     val.end[val.nex != 0] = val.nex[val.nex != 0]
  #   } else {
  #     # stop()
  #     b.end = c(b.end, icol-1)
  #     b.beg = c(b.beg, icol)
  #     
  #     mx.beg = cbind(mx.beg, val.beg)
  #     mx.end = cbind(mx.end, val.end)
  #     
  #     if(sum(val.beg == 0) != sum(val.beg * val.end == 0)) stop()
  #     
  #     val.beg = val.nex
  #     val.end = val.nex
  #   }
  # }
  # b.end = c(b.end, mx.len)  # last ending
  # mx.beg = cbind(mx.beg, val.beg)
  # mx.end = cbind(mx.end, val.end)
  # 
  # # Common data.frame
  # blocks = data.frame(beg = b.beg, end = b.end)
  # blocks$len = blocks$end - blocks$beg + 1
  # 
  # # Sorting
  # n = nrow(blocks)
  # index <- 1:n 
  # i <- 2
  # while (i <= n) {
  #   
  #   beg.i <- mx.beg[,index[i]]
  #   end.i <- mx.end[,index[i]]
  #   
  #   key.i <- index[i]
  #   k <- 0  # position which fit the
  #   for (j in (i - 1):1) {
  #     beg.j = mx.beg[,index[j]]
  #     end.j = mx.end[,index[j]]
  #     
  #     idx.common = (beg.j * beg.i) != 0
  #     
  #     if(sum(idx.common) == 0){
  #       next
  #     } else if(max(abs(beg.j[idx.common] - end.i[idx.common])) == 1) {
  #       k = j
  #       break
  #     } else if(max(abs(beg.i[idx.common] - end.j[idx.common])) == 1) {
  #       k = j+1
  #       break
  #     }else if(sum(beg.i[idx.common] > end.j[idx.common])) {
  #       break  
  #     }
  #   }
  #   # stop()
  #   if(k == 0){
  #     # stop('skip')
  #     tmp = index[i]
  #     index[i:(n-1)] = index[(i+1):n]
  #     index[n] = tmp
  #     next # not to increase the i index
  #   } else if (k != i) {
  #     # stop('middle fit')
  #     tmp = index[i]
  #     index[(k+1):i] = index[k:(i-1)]
  #     index[k] = tmp
  #     
  #   } else {
  #     # stop('stay')
  #   }
  #   i <- i + 1  
  #   
  # }
  # 
  # pos = c()
  # for(i in index){
  #   pos = c(pos, blocks$beg[i]:blocks$end[i])
  # }
  # 
  # mx = mx[,pos]
  # mx.pos = mx.pos[,pos]
  # 
  # # Check that every sequence is sorted
  # for(irow in 1:nrow(mx.pos)){
  #   p = mx.pos[irow,]
  #   p = p[p != 0]
  #   if(is.unsorted(p)) print(irow)
  #   
  #   # which( mx.pos[irow,] %in% which(diff(p) !=1))
  # }
  # 
  # # p4 = msaplot(mx)
  # # p4
  
  
  
  mx.pos.real = mx.pos * 0
  for(irow in 1:nrow(mx.pos)){
    tmp = as.numeric(strsplit(rownames(mx.pos)[irow], '\\|')[[1]][3:4])
    idx.nogap = which(mx.pos[irow] != 0)
    mx.pos.real[irow, idx.nogap] = tmp[mx.pos[irow, idx.nogap]]
  }
  
  return(list(mx = mx, pos = mx.pos.real))
  
}


#' Mask Unaligned Blocks in The Alignment Matrix
#'
#' This function identifies a masks unaligned or poorly aligned blocks in
#' a matrix of nucleotides. It evaluates the diversity at each position and
#' masks blocks based on a diversity cutoff.
#'
#' @param mx A character matrix where rows represent sequences and columns
#' represent aligned positions. Sequences should contain only 'A', 'C', 'G',
#' 'T', and '-' for gaps.
#' @param sim.cutoff A numeric threshold for masking blocks based on diversity.
#' Blocks with diversity below this cutoff are considered poorly aligned and masked.
#'
#' @return A logical matrix of the same dimensions as `mx`, where `TRUE` indicates
#' that the position is well aligned and `FALSE` indicates it should be masked.
#'
#' @export
maskUnaligned <- function(mx, sim.cutoff = 0.2){
  
  s.nts = c('A', 'C', 'G', 'T')
  
  # Diversity by each position
  pos.profile = getProfile(mx)
  pos.variation = (colSums(pos.profile == 0) != 3) * 1
  
  # Define blocks, were the alignment non well
  blocks.all = c()
  for(i.seq in 1:nrow(mx)){
    s = mx[i.seq,]
    blocks = findOnes((s != '-') *1)
    blocks$len = blocks$end - blocks$beg + 1
    if(nrow(blocks) == 0) next
    
    # Estimate diversity within each block
    blocks$pi = 0
    blocks$acc = rownames(mx)[i.seq]
    for(irow in 1:nrow(blocks)){
      idx.block = blocks$beg[irow]:blocks$end[irow]
      blocks$pi[irow] = sum(pos.variation[idx.block]) / blocks$len[irow]
    }
    
    blocks.all = rbind(blocks.all, blocks)
  }
  
  
  if(nrow(blocks) == 0){
    pokaz('Exit 1')
    return(NULL)
  } 
  
  blocks = blocks.all
  rm(blocks.all)
  blocks = blocks[blocks$pi>sim.cutoff,,drop = F]
  blocks = blocks[order(-blocks$pi),,drop = F]
  
  if(nrow(blocks) == 0){
    pokaz('Exit 2')
    return(NULL)
  }
  
  blocks$remove = 0
  for(irow in 1:nrow(blocks)){
    idx.block = blocks$beg[irow]:blocks$end[irow]
    blocks$pi[irow] = sum(pos.variation[idx.block]) / blocks$len[irow]
    if(blocks$pi[irow] <= sim.cutoff){
      # stop()
      # pokaz(irow)
      next
    } 
    
    seq.tmp = mx[blocks$acc[irow], idx.block]  
    for(s.nt in s.nts){
      pos.profile[s.nt,idx.block] = pos.profile[s.nt,idx.block] - 1 * (seq.tmp == s.nt)
    }
    
    blocks$remove[irow] = 1
    pos.variation[idx.block] = (colSums(pos.profile[,idx.block,drop = F] == 0) != 3) * 1
  }
  
  # Create a mask for blocks with very high diversity
  mask.ok = !(mx == 'X')  # a way to make a matrix with all True
  for(irow in 1:nrow(blocks)){
    if(blocks$remove[irow] == 0) next
    idx.block =  blocks$beg[irow]:blocks$end[irow]
    mask.ok[blocks$acc[irow],idx.block] = F
  }
  
  return(mask.ok)
}



#' Calculate Nucleotide Profile for Each Position in The Alignment Matrix
#'
#' This function computes a profile matrix showing the count of each nucleotide
#' ('A', 'C', 'G', 'T') at every position in a given matrix of DNA sequences.
#'
#' @param mx Alignment matrix where each row represents a sequence and each
#' column represents a nucleotide position in alignment.
#'
#' @return A numeric matrix with 4 rows, each representing one of the nucleotides
#' ('A', 'C', 'G', 'T'). The columns correspond to the positions in the input alignment,
#' and the values represent the count of each nucleotide at each position.
#'
#' @export
getProfile <- function(...) {
  pokazAttention('Please replace `getProfile` with `mx2profile`')
  result <- mx2profile(...)
  return(result)
}















