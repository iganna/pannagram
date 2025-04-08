refineAlignment <- function(seqs.clean, path.work, mx.cons1){
  
  n.seqs = length(seqs.clean)
  
  # # ---- Distance matrix ----
  # 
  # dist.mx = calсDistAln(seqs.mx)
  dist.mx = calcDistKmer(seqs.clean)
  # 
  # # ---- Clustering ----
  hc = hclust(as.dist(dist.mx))
  # 
  sim.cutoff = 0.1
  clusters <- cutree(hc, h = sim.cutoff) # 0.1 was before
  
  # dist.mx = calcDistKmer(seqs.clean)
  # clusters = clustDistKmer(dist.mx)
  
  
  # ---- Sequences from clusters ----
  
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
    
    writeFasta(seqs.tmp, seqs.cl.fasta)
    
    #--op 3 --ep 0.1
    system(paste('mafft  --quiet --maxiterate 100 ', seqs.cl.fasta, '>', aln.fasta,  sep = ' '))
    
    seqs.cl.aln = readFasta(aln.fasta)
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
  
  # ---- Extend the first alignment ----
  name1 = names(seqs.clean)[1]
  for(i.aln in 1:length(alignments)){
    if(!(name1 %in% rownames(alignments))) next
    
    new.rownames = paste0('new_tmp_', 1:nrow(mx.cons1))
    
    mx.cons1.aln = matrix('-', nrow = nrow(mx.cons1), ncol = ncol(alignments[[i.aln]]),
                          dimnames = c(new.rownames, NULL))
    pos.i.aln = which(alignments[[i.aln]][name1,] != '-')
    mx.cons1.aln[,pos.i.aln] = mx.cons1
    
    mx.cons1.pos = matrix(0, nrow = nrow(mx.cons1), ncol = ncol(positions[[i.aln]]),
                          dimnames = c(new.rownames, NULL))
    for(irow in 1:nrow(mx.cons1.pos)){
      idx.irow = mx.cons1[irow] != '-'
      mx.cons1.pos[,pos.i.aln[idx.irow]] = 1:sum(idx.irow)
    }
    
    alignments[[i.aln]] = rbind(alignments[[i.aln]], mx.cons1.aln)
    
    positions[[i.aln]] = rbind(positions[[i.aln]], mx.cons1.pos)
    
    idx.remove = which(rownames(alignments[[i.aln]]) == name1)
    
    alignments[[i.aln]] = alignments[[i.aln]][-idx.remove,]
    positions[[i.aln]] = positions[[i.aln]][-idx.remove,]
    
    break
  }
  
  # ---- Merge clusters ----
  
  for(i.merge in which(df.merge$cl > max(clusters))){
    pokaz(i.merge, df.merge$cl[i.merge])
    
    # if(i.merge == 23) stop()
    
    i.cl1 = df.merge$id1[i.merge]
    i.cl2 = df.merge$id2[i.merge]
    s1 = seqs.cl[i.cl1] 
    s2 = seqs.cl[i.cl2] 
    
    # dotplot.s(s1, s2, 15, 14)
    
    ## ---- Mafft add alignments ----
    
    # Reduce the number of sequences
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
    
    # msaplot(mafft.mx)
    
    irow_support = c()
    cnunk.min.len = 25
    sim.score = 0.9
    for(irow in 1:nrow(result)){
      if(result$len[irow] < min(c(nchar(s1)/3, nchar(s2)/3, cnunk.min.len))) next
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