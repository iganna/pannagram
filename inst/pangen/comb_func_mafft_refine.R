
refineMafft <- function(mx, n.flank = 30){
  
  mx = toupper(mx)
  
  sim.cutoff = 0.2
  s.nt.fake = '!'
  
  
  # Create the matrix with positions
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


# path.work = '/Volumes/Samsung_T5/vienn/test/pannagram_test/mafft/'
refineAlignment <- function(seqs, path.work, n.flank = 30){
  
  seqs.mx = aln2mx(seqs)
  n.seqs = length(seqs)
  
  # ---- Distance matrix ----
  
  dist.mx = calсDistanceMatrix(seqs.mx)
  
  # ---- Clustering ----
  hc = hclust(as.dist(dist.mx))
  
  clusters <- cutree(hc, h = 0.1)
  
  # ---- Sequences from clusters ----
  
  seqs.clean = sapply(seqs, 
                      function(seq) substr(gsub('-', '', seq), 
                                           n.flank+1, 
                                           nchar(seq) - n.flank))

  writeFastaMy(seqs.clean, paste0(path.work, 'seqs_clean.fasta'))
  
  seqs.cl = c()
  
  positions = list()
  alignments = list()
  for(i.cl in 1:max(clusters)){
    pokaz(i.cl)
    seqs.tmp = seqs.clean[names(clusters)[clusters == i.cl]]
    
    if(length(seqs.tmp) == 1) {
      seqs.cl = c(seqs.cl, seqs.tmp)
      positions[[i.cl]] = matrix(1:nchar(seqs.tmp), nrow = 1, 
                                 dimnames = list(names(seqs.tmp), NULL))
      alignments[[i.cl]] = matrix(seq2nt(seqs.tmp), nrow = 1, 
                                  dimnames = list(names(seqs.tmp), NULL))
      next
    }
    
    seqs.cl.fasta = paste0(path.work, 'seqs_', i.cl, '.fasta')
    aln.fasta = paste0(path.work, 'aln_', i.cl, '.fasta')
    
    writeFastaMy(seqs.tmp, seqs.cl.fasta)
    
    system(paste('mafft --op 3 --ep 0.1 --quiet --maxiterate 100 ', seqs.cl.fasta, '>', aln.fasta,  sep = ' '))
    
    seqs.cl.aln = readFastaMy(aln.fasta)
    seqs.cl.mx = aln2mx(seqs.cl.aln)
    
    pos.cl.mx = matrix(0,nrow = nrow(seqs.cl.mx), ncol = ncol(seqs.cl.mx))
    for(irow in 1:nrow(seqs.cl.mx)){
      pos.cl.mx[irow,seqs.cl.mx[irow,] != '-'] = 1:nchar(seqs.tmp[irow])
    }
    rownames(pos.cl.mx) = names(seqs.tmp)
    positions[[i.cl]] = pos.cl.mx
    alignments[[i.cl]] = seqs.cl.mx
    
    
    # msadiff(seqs.cl.mx)
    seq.cons = mx2cons(seqs.cl.mx)
    
    seqs.cl = c(seqs.cl, nt2seq(seq.cons))
  }
  names(seqs.cl) = paste0('clust_', 1:length(seqs.cl))
  
  
  # ---- Hierarchical clustering ----
  df.merge = as.data.frame(hc$merge)
  df.merge$id1 = ifelse(hc$merge[,1] < 0, clusters[abs(hc$merge[,1])], NA)
  df.merge$id2 = ifelse(hc$merge[,2] < 0, clusters[abs(hc$merge[,2])], NA)
  
  df.merge$cl = rep(0, nrow(hc$merge))
  n.cl = max(clusters) + 1
  
  for (i in 1:nrow(df.merge)) {
    df.merge$id1[i] = ifelse(is.na(df.merge$id1[i]), df.merge$cl[abs(hc$merge[i,1])], df.merge$id1[i])
    df.merge$id2[i] = ifelse(is.na(df.merge$id2[i]), df.merge$cl[abs(hc$merge[i,2])], df.merge$id2[i])
    
    df.merge$cl[i] = ifelse(df.merge$id1[i] == df.merge$id2[i], df.merge$id1[i], n.cl)
    
    if (df.merge$cl[i] == n.cl) {
      n.cl = n.cl + 1
    }
  }
  
  # ---- Merge clusters ----
  mx.pos.init = aln2pos(seqs)
  
  for(i.merge in which(df.merge$cl > max(clusters))){
    pokaz(i.merge)
    
    # if(i.merge == 23) stop()
    
    i.cl1 = df.merge$id1[i.merge]
    i.cl2 = df.merge$id2[i.merge]
    s1 = seqs.cl[i.cl1] 
    s2 = seqs.cl[i.cl2] 
    
    # dotplot.s(s1, s2, 15, 14)

    ## ---- Mafft add alignments ----
    # seq1 = mx2aln(alignments[[i.cl1]])
    # seq2 = mx2aln(alignments[[i.cl2]])
    
    # # Reduce the number of sequences
    seq1 = mx2aln(mx2cons(alignments[[i.cl1]], amount = 3))
    seq2 = mx2aln(mx2cons(alignments[[i.cl2]], amount = 3))
    
    pokaz('before', nrow(alignments[[i.cl1]]), nrow(alignments[[i.cl2]]))
    pokaz('after', length(seq1), length(seq2))
    
    mafft.res = mafftAdd(seq1, seq2, path.work)
    
    result = mafft.res$result
    df = mafft.res$df
    pos.mx = mafft.res$pos.mx
    
    
    ## ---- Check all the synteny chunks ----
    mafft.mx = matrix('-', nrow = 2, ncol = ncol(pos.mx))
    mafft.mx[1, pos.mx[1,] != 0] = toupper(seq2nt(s1))
    mafft.mx[2, pos.mx[2,] != 0] = toupper(seq2nt(s2))
    
    irow_support = c()
    cnunk.min.len = 50
    sim.score = 0.9
    for(irow in 1:nrow(result)){
      if(result$len[irow] < cnunk.min.len) next
      idx = result$beg[irow]:result$end[irow]
      score = sum(mafft.mx[1,idx] == mafft.mx[2,idx]) / result$len[irow]
      if(score >= sim.score){
        irow_support = c(irow_support, irow)
      }
    }

    ## ---- BLAST----
  
    x = data.frame(tmp=numeric())   
    if(length(irow_support) != nrow(result)){
      x = blastTwoSeqs(s1, s2, path.work)  
    }
    
    # x = x[x$V7 > 30,]
    if(nrow(x) > 0){
      
      x$idx = 1:nrow(x)
      x$dir = (x$V4 > x$V5) * 1
      x = x[x$dir == 0,]   # TODO: keep inversions
      x = glueZero(x)
      
      # test df lines
      for (irow in setdiff(1:nrow(df), irow_support)) {
        
        start_df <- df$V2[irow]
        end_df <- df$V3[irow]
        
        length_df <- end_df - start_df
        
        overlap_start <- pmax(start_df, x$V2) 
        overlap_end <- pmin(end_df, x$V3)     
        
        overlap_length <- overlap_end - overlap_start
        
        valid_overlap <- overlap_length >= 0.9 * length_df & overlap_length > 0
        
        if(sum(valid_overlap) == 0) next
        
        df.tmp = x[valid_overlap,]
        
        df.tmp[,c('V4', 'V5')] = df.tmp[,c('V4', 'V5')] - df$V4[irow] + df$V2[irow]
        
        valid_overlap2 = pmin(df.tmp$V3, df.tmp$V5) - pmax(df.tmp$V2, df.tmp$V4)
        valid_overlap2 = valid_overlap2 / df.tmp$V7
        
        if(sum(valid_overlap2 > 0.9) > 0) {
          irow_support = c(irow_support, irow)
        }
        
      }
    }
    
    if(length(irow_support) != 0){
      irow_support = sort(irow_support)
      
      n.sup = length(irow_support)
      
      result.sup = result[irow_support,,drop=F]
      result.sup$type = 0
      result.sup$id = 1:n.sup
      result.sup$len = result.sup$end - result.sup$beg + 1
      
      df.sup = df[irow_support,,drop=F]
      
      # stop()
      
      # ---- Add the gaps from the first sequence ----
      len1 = sum(pos.mx[1,] != 0)
      result.sup.tmp = data.frame(beg = c(1, df.sup$V3+1),
                                  end = c(df.sup$V2-1, len1),
                                  # id = (0:(n.sup)) + 0.3),
                                  id = (0:(n.sup)),
                                  type=1)
      
      # Lengths of all blocks
      result.sup.tmp$len = result.sup.tmp$end - result.sup.tmp$beg + 1
      result.sup.tmp = result.sup.tmp[result.sup.tmp$len > 0,,drop = F]
      
      # define ID by the order in the initial alignment
      if(nrow(result.sup.tmp) > 0){
        len.mafft.aln = ncol(pos.mx)
        for(irow in 1:nrow(result.sup.tmp)){
          tmp = which((result.sup.tmp$beg[irow] <= pos.mx[1,]) & (pos.mx[1,] <= result.sup.tmp$end[irow]))
          result.sup.tmp$id[irow] = result.sup.tmp$id[irow] + (mean(tmp)) / len.mafft.aln
        }
        result.sup = rbind(result.sup, result.sup.tmp)        
      }
      
      
      # ---- Add the gaps from the second sequence ----
      len2 = sum(pos.mx[2,] != 0)
      result.sup.tmp =  data.frame(beg = c(1, df.sup$V5+1),
                                   end = c(df.sup$V4-1, len2),
                                   # id = (0:(n.sup)) + 0.5),
                                   id = (0:(n.sup)),
                                   type=2)
      
      # Lengths of all blocks
      result.sup.tmp$len = result.sup.tmp$end - result.sup.tmp$beg + 1
      result.sup.tmp = result.sup.tmp[result.sup.tmp$len > 0,,drop = F]
      
      # define ID by the order in the initial alignment
      if(nrow(result.sup.tmp) > 0){
        len.mafft.aln = ncol(pos.mx)
        for(irow in 1:nrow(result.sup.tmp)){
          tmp = which((result.sup.tmp$beg[irow] <= pos.mx[2,]) & (pos.mx[2,] <= result.sup.tmp$end[irow]))
          result.sup.tmp$id[irow] = result.sup.tmp$id[irow] + mean(tmp) / len.mafft.aln
        }
        result.sup = rbind(result.sup, result.sup.tmp)
      }
      
      
      # ---- Order ----
      result.sup = result.sup[order(result.sup$id),]
      
      result.sup$a.beg = cumsum(c(0, result.sup$len[-nrow(result.sup)])) + 1
      result.sup$a.end = result.sup$a.beg + result.sup$len - 1
      
      if(sum((result.sup$a.end - result.sup$a.beg + 1) != result.sup$len) > 0) {
        stop("SOMETHING IS WROTG")
      }
      
      # ---- Consensus positions ----
      
      mx.comb = matrix(0, nrow = 2, ncol = result.sup$a.end[nrow(result.sup)])
      for(irow in 1:nrow(result.sup)){
        if(result.sup$type[irow] == 0){
          mx.comb[,result.sup$a.beg[irow]:result.sup$a.end[irow]] = 
            pos.mx[,result.sup$beg[irow]:result.sup$end[irow]]
        } else {
          
          mx.comb[result.sup$type[irow],result.sup$a.beg[irow]:result.sup$a.end[irow]] = 
            result.sup$beg[irow]:result.sup$end[irow]
        }
      }
      
      rowSums(pos.mx != 0)
      
      if(sum(rowSums(mx.comb != 0) !=  rowSums(pos.mx != 0)) > 0){
        stop('MISSED POSITIONS')
        setdiff(pos.mx[1,], mx.comb[1,])
        setdiff(pos.mx[2,], mx.comb[2,])
      } 
      
      # stop()
      
      for(irow in 1:2){
        pp = mx.comb[irow,]
        pp = pp[pp != 0]
        if(is.unsorted(pp)){
          stop('WRONG SORTING')
        }
        
        if(length(pp) != length(unique(pp))){
          stop('WRONG UNIQUE')
        }
      }
      
    } else {
      
      n1 = nchar(s1)
      n2 = nchar(s2)
      mx.comb = matrix(0, nrow = 2, ncol = n1 + n2)
      mx.comb[1, 1:n1] = 1:n1
      mx.comb[2, n1 + (1:n2)] = 1:n2
    }
    
    # ---- Consensus sequences ----
    
    s1 = seq2nt(s1)
    s2 = seq2nt(s2)
    
    non.zero.indices.1 <- mx.comb[1,] != 0
    non.zero.indices.2 <- mx.comb[2,] != 0
    
    n1 = nrow(positions[[i.cl1]])
    n2 = nrow(positions[[i.cl2]])
    mx.comb.pos = matrix(0, 
                         nrow = n1 + n2,
                         ncol = ncol(mx.comb))
    mx.comb.pos[1:n1, non.zero.indices.1]        = positions[[i.cl1]][,mx.comb[1, non.zero.indices.1]]
    mx.comb.pos[n1 + (1:n2), non.zero.indices.2] = positions[[i.cl2]][,mx.comb[2, non.zero.indices.2]]
    
    mx.comb.seq = matrix('-', 
                         nrow = n1 + n2,
                         ncol = ncol(mx.comb))
    
    mx.comb.seq[1:n1, non.zero.indices.1]        = alignments[[i.cl1]][,mx.comb[1, non.zero.indices.1]]
    mx.comb.seq[n1 + (1:n2), non.zero.indices.2] = alignments[[i.cl2]][,mx.comb[2, non.zero.indices.2]]                     
    
    
    tmp.names = c(rownames(alignments[[i.cl1]]), rownames(alignments[[i.cl2]]))
    rownames(mx.comb.pos) = tmp.names
    rownames(mx.comb.seq) = tmp.names
    
    positions[[df.merge$cl[i.merge]]] = mx.comb.pos
    alignments[[df.merge$cl[i.merge] ]] = mx.comb.seq
    
    # if(i.merge == 21) stop()
    # 
    # msaplot(mx.comb.seq)
    # msadiff(mx.comb.seq)
    # 
    
    seq.cons.comb = mx2cons(mx.comb.seq)
    
    seqs.cl[paste0('clust_', df.merge$cl[i.merge])] = nt2seq(seq.cons.comb)
    
  }
  
  return(list(pos = positions,
              aln = alignments))
  
}

#' Align two Alignments with MAFFT
#'
#' This function performs sequence alignment using the MAFFT tool. It takes two sets of alignments,
#' aligns them using MAFFT, and then identifies common aligned regions.
#'
#' @param seq1 Character vector of the first set of aligned sequences.
#' @param seq2 Character vector of the second set of aligned sequences.
#' @param path.work String representing the working directory path where intermediate files will be stored.
#' @param n.diff Integer specifying the maximum allowable gap between alignment blocks to be merged. Default is 5.
#' @return A data frame with the alignment positions and lengths of the identified blocks.
#' @export
mafftAdd <- function(seq1, seq2, path.work, n.diff = 5){
  
  # Create a tablefile (so-called subMSAtable) for mafft 
  id1 = 1:length(seq1)
  id2 = (length(seq1) + 1):(length(seq1) + length(seq2))
  
  first_line <- paste(id1, collapse = " ")
  second_line <- paste(id2, collapse = " ")
  file.table = paste0(path.work, "tbl.txt")
  writeLines(c(first_line, second_line, ''), file.table)
  
  # Merge sequences into the input file
  file.seqs.merge = paste0(path.work, 'seqs_tmp_merge.fasta')
  writeFastaMy(c(seq1, seq2), file.seqs.merge)
  
  # File for storing MAFFT alignment output
  file.mafft = paste0(path.work, 'seqs_tmp_mafft.fasta')
  
  # Run MAFFT with the specified table and merged sequences
  system(paste('mafft --op 3  --ep 0.1 --quiet --merge', file.table, file.seqs.merge, '>', file.mafft, sep = ' '))
  
  # Read the aligned sequences
  xx = readFastaMy(file.mafft)
  mafft.mx = aln2mx(xx)
  
  # p = msaplot(mafft.mx)
  
  # Initialize position matrix to map alignment positions
  pos.mx = matrix(0, nrow = 2, ncol = ncol(mafft.mx))
  
  # Identify positions with sequence content (no gaps) in the alignment
  pos1 = colSums(mafft.mx[id1,,drop=F] != '-') > 0
  pos2 = colSums(mafft.mx[id2,,drop=F] != '-') > 0
  
  pos.mx[1, pos1] = 1:sum(pos1)
  pos.mx[2, pos2] = 1:sum(pos2)
  
  # Identify aligned blocks in the merged alignment
  pos.x = pos1 * 2 + pos2
  result = findOnes((pos.x == 3) * 1)
  
  # Merge close alignment blocks based on `n.diff` threshold
  i <- 1
  while (i < nrow(result)) {
    if (result$beg[i + 1] - result$end[i] <= n.diff) {
      result$end[i] <- result$end[i + 1]
      result <- result[-(i + 1), ]      
    } else {
      i <- i + 1
    }
  }
  result$len = result$end - result$beg + 1
  
  # Create a data frame with the final alignment block information
  df = data.frame(
    V2 = pos.mx[1, result$beg],    # Start position in seq1
    V3 = pos.mx[1, result$end],    # End position in seq1
    V4 = pos.mx[2, result$beg],    # Start position in seq2
    V5 = pos.mx[2, result$end],    # End position in seq2
    V7 = result$end - result$beg + 1 # Length of the alignment block
  )
  
  return(list(result = result,
              df = df,
              pos.mx = pos.mx))
}


#' Perform BLAST Alignment Between Two Sequences
#'
#' This function performs a BLAST alignment between two input sequences.
#'
#' @param s1 Character string representing the first sequence.
#' @param s2 Character string representing the second sequence.
#' @param path.work String representing the working directory path where intermediate files will be stored.
#' @return A data frame containing BLAST alignment results with columns for sequence identifiers, 
#' alignment positions, percentage identity, and alignment lengths.
#' @export
blastTwoSeqs <- function(s1, s2, path.work){
  
  # Write the sequences to a temporary FASTA file
  file.seqs.tmp = paste0(path.work, 'seqs_tmp.fasta')
  writeFastaMy(c(s1, s2), file.seqs.tmp)
  
  # Define the output file for BLAST results
  file.blast.cons = paste0(file.seqs.tmp, '.out')
  
  # Create a BLAST database from the temporary sequence file
  system(paste0('makeblastdb -in ', file.seqs.tmp,' -dbtype nucl  > /dev/null'))
  
  # Run BLASTn alignment between the sequences in the temporary file
  system(paste0('blastn -db ',file.seqs.tmp,' -query ',file.seqs.tmp,
                ' -num_alignments 50 ',
                ' -out ',file.blast.cons,
                ' -outfmt "7 qseqid qstart qend sstart send pident length sseqid"'))
  
  # Read the BLAST output into a data frame
  x = readBlast(file.blast.cons)
  x = x[x$V1 != x$V8,,drop=F] # Remove self-alignments
  
  # Split data into two groups: alignments starting from sequence 1 and sequence 2
  idx1 = x$V1 %in% names(s1)
  x1 = x[idx1,,drop=F]
  x2 = x[!idx1,,drop=F]
  
  # Swap columns for the second group to match the first group's format
  colnames(x2)[2:5] = c('V4', 'V5', 'V2', 'V3') 
  colnames(x2)[c(1,8)] = colnames(x2)[c(8,1)]
  
  # Merge both groups into a single data frame
  x = rbind(x1, x2)
  x = unique(x) # Remove duplicates
  
  return(x)
}



# s1.blast = s1
# s2.blast = s2
# 
# for(irow in 1:nrow(df)){
#   
#   if(df$V7[irow] < 50) next
#   
#   sx1 = seq2nt(s1)
#   sx1 = sx1[df$V2[irow]:df$V3[irow]]
#   sx1 = nt2seq(sx1)
#   
#   s1.blast[paste0('id_', irow, '_1_', df$V2[irow]-1)] = sx1
#   
#   sx2 = seq2nt(s2)
#   sx2 = sx2[df$V4[irow]:df$V5[irow]]
#   sx2 = nt2seq(sx2)
#   
#   s2.blast[paste0('id_', irow, '_1_', df$V4[irow]-1)] = sx2
#   
# }
# 
# xx = blastTwoSeqs2(s1.blast, s2.blast, path.work)


blastTwoSeqs2 <- function(s1, s2, path.work){
  
  # Write the sequences to a temporary FASTA file
  file.seqs.tmp1 = paste0(path.work, 'seqs_tmp1.fasta')
  writeFastaMy(s1, file.seqs.tmp1)
  
  file.seqs.tmp2 = paste0(path.work, 'seqs_tmp2.fasta')
  writeFastaMy(s2, file.seqs.tmp2)
  
  # Define the output file for BLAST results
  file.blast1 = paste0(file.seqs.tmp1, '.out')
  file.blast2 = paste0(file.seqs.tmp2, '.out')
  
  # Create a BLAST database from the temporary sequence file
  system(paste0('makeblastdb -in ', file.seqs.tmp1,' -dbtype nucl  > /dev/null'))
  system(paste0('makeblastdb -in ', file.seqs.tmp2,' -dbtype nucl  > /dev/null'))
  
  # Run BLASTn alignment between the sequences in the temporary file
  system(paste0('blastn -db ',file.seqs.tmp2,' -query ',file.seqs.tmp1,
                ' -num_alignments 50 ',
                ' -out ',file.blast1,
                ' -outfmt "7 qseqid qstart qend sstart send pident length sseqid"'))
  
  system(paste0('blastn -db ',file.seqs.tmp1,' -query ',file.seqs.tmp2,
                ' -num_alignments 50 ',
                ' -out ',file.blast2,
                ' -outfmt "7 qseqid qstart qend sstart send pident length sseqid"'))
  
  # Read the BLAST output into a data frame
  x1 = readBlast(file.blast1)
  x2 = readBlast(file.blast2)
  
  # Swap columns for the second group to match the first group's format
  colnames(x2)[2:5] = c('V4', 'V5', 'V2', 'V3') 
  colnames(x2)[c(1,8)] = colnames(x2)[c(8,1)]
  
  # Merge both groups into a single data frame
  x = rbind(x1, x2)
  x = unique(x) # Remove duplicates
  
  return(x)
}


#' Calculate Distance Matrix for Sequence Alignment
#'
#' This function calculates a pairwise distance matrix for a given set of aligned sequences.
#'
#' @param seqs.mx A matrix where each row represents a sequence, and each column represents a position in the alignment.
#' Gaps in the alignment are represented as `'-'`.
#'
#' @return A symmetric matrix of pairwise distances between the sequences.
#' @export
#'
calсDistanceMatrix <- function(seqs.mx) {
  
  n.seqs = nrow(seqs.mx)
  
  dist.mx <- matrix(0, nrow = n.seqs, ncol = n.seqs, 
                    dimnames = list(rownames(seqs.mx), rownames(seqs.mx)))
  
  # Precompute non-gap positions
  valid.pos.all <- lapply(1:n.seqs, function(i) seqs.mx[i,] != '-')
  
  for(i in 1:(n.seqs-1)) {
    for(j in (i+1):n.seqs) {
      
      k <- ifelse(sum(valid.pos.all[[i]]) > sum(valid.pos.all[[j]]), i, j)
      
      # Calculate the positions that are valid for comparison (non-gap positions)
      valid.pos <- valid.pos.all[[i]] | valid.pos.all[[j]]
      
      # Distance
      d <- sum((seqs.mx[i,] != seqs.mx[j,]) & valid.pos) / sum(valid.pos)
      
      dist.mx[i, j] <- d
      dist.mx[j, i] <- d
    }
  }
  
  return(dist.mx)
}






