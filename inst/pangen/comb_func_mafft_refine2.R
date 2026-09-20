#' Refine Sequence Alignment
#'
#' This function refines sequence alignments by clustering sequences, removing flanking gaps, 
#' and merging aligned clusters based on similarity and synteny analysis.
#'
#' @param seqs A list of aligned nucleotide sequences.
#' @param path.work A character string specifying the working directory where temporary files 
#' will be saved.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{pos}{A list of matrices representing the positions of nucleotides in the refined alignments for each cluster.}
#'   \item{aln}{A list of matrices representing the refined alignments for each cluster.}
#' }
#'
#' @export
refineAlignment_prev <- function(seqs.clean, path.work, keep.pos = FALSE){

  n.seqs = length(seqs.clean)
  seqs.clean = toupper(seqs.clean)
  
  # ---- Handle N in sequences ----
  # n.n = sapply(gregexpr("N", seqs.clean, ignore.case = TRUE), function(m) sum(m > 0))
  # seqs.names.n = names(seqs.clean)[n.n != 0]
  # if(length(seqs.names.n) > 0){
  #   pos.replace = list()
  #   for(s.name in seqs.names.n){
  #     s = seqs.clean[s.name]
  #     s = seq2nt(s)
  #     pos.replace[[s.name]] = which((s != 'n') & (s != 'N'))
  #     s = s[pos.replace[[s.name]]]
  #     seqs.clean[s.name] = nt2seq(s)
  #   }
  # }
  
  # ---- Distance matrix ----
  # 
  # dist.mx = calcDistAln(seqs.mx)
  dist.mx = calcDistKmer(seqs.clean)
  
  # save(list = ls(), file = "tmp_workspace_refine_aln.RData")
  
  # ---- Clustering ----
  hc = hclust(as.dist(dist.mx))
  
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
      if(keep.pos) positions[[i.cl]] = matrix(1:nchar(seqs.tmp), nrow = 1, 
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
    
    if(keep.pos){
      pos.cl.mx = matrix(0L,nrow = nrow(seqs.cl.mx), ncol = ncol(seqs.cl.mx))
      for(irow in 1:nrow(seqs.cl.mx)){
        pos.cl.mx[irow,seqs.cl.mx[irow,] != '-'] = 1:nchar(seqs.tmp[irow])
      }
      rownames(pos.cl.mx) = names(seqs.tmp)
      positions[[i.cl]] = pos.cl.mx
    }
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
  
  merge.rows = which(df.merge$cl > max(clusters))
  for(i.merge in merge.rows){
    
    # if(df.merge$cl[i.merge] == 17) stop()
    pokaz(i.merge, df.merge$cl[i.merge])
    
    
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
    
    df = mafft.res$df
    pos.mx = mafft.res$pos.mx
    result = mafft.res$result
    result$support = rep(0, nrow(result))  # 0-row result -> no support branch below
    
    len.chunk.merge = 8
    irow = 1
    while(irow < nrow(result)){
      irow.d = result$beg[irow + 1] - result$end[irow] - 1
      if(irow.d <= len.chunk.merge){
        # Merge chunks
        result$end[irow] = result$end[irow + 1]
        result$len[irow] = result$end[irow] - result$beg[irow] + 1
        result = result[-(irow + 1),]
        
        result$support[irow] = result$support[irow] + irow.d
        
        df$V3[irow] = df$V3[irow + 1]
        df$V5[irow] = df$V5[irow + 1]
        df = df[-(irow + 1),]
        df$len = result$len
      }
      irow = irow + 1
    }
    
    
    ## ---- Check all the synteny chunks ----
    mafft.mx = matrix('-', nrow = 2, ncol = ncol(pos.mx))
    mafft.mx[1, pos.mx[1,] != 0] = toupper(seq2nt(s1))
    mafft.mx[2, pos.mx[2,] != 0] = toupper(seq2nt(s2))
    
    # msaplot(mafft.mx)
    
    irow_support = c()
    cnunk.min.len = 25
    sim.score = 0.8
    for(irow in seq_len(nrow(result))){
      if((irow != 1) & (irow != nrow(result))){
        if(result$len[irow] < min(c(nchar(s1)/3, nchar(s2)/3, cnunk.min.len))) next  
      }
      idx = result$beg[irow]:result$end[irow]
      score = sum(mafft.mx[1,idx] == mafft.mx[2,idx]) / (result$len[irow] - result$support[irow])
      if(score >= sim.score){
        irow_support = c(irow_support, irow)
      }
    }
    result$support = NULL
    
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
    
    n1 = nrow(alignments[[i.cl1]])
    n2 = nrow(alignments[[i.cl2]])
    if(keep.pos){
      mx.comb.pos = matrix(0L, 
                           nrow = n1 + n2,
                           ncol = ncol(mx.comb))
      mx.comb.pos[1:n1, non.zero.indices.1]        = positions[[i.cl1]][,mx.comb[1, non.zero.indices.1]]
      mx.comb.pos[n1 + (1:n2), non.zero.indices.2] = positions[[i.cl2]][,mx.comb[2, non.zero.indices.2]]
    }
    
    mx.comb.seq = matrix('-', 
                         nrow = n1 + n2,
                         ncol = ncol(mx.comb))
    
    mx.comb.seq[1:n1, non.zero.indices.1]        = alignments[[i.cl1]][,mx.comb[1, non.zero.indices.1]]
    mx.comb.seq[n1 + (1:n2), non.zero.indices.2] = alignments[[i.cl2]][,mx.comb[2, non.zero.indices.2]]                     
    
    
    tmp.names = c(rownames(alignments[[i.cl1]]), rownames(alignments[[i.cl2]]))
    rownames(mx.comb.seq) = tmp.names
    if(keep.pos){
      rownames(mx.comb.pos) = tmp.names
      positions[[df.merge$cl[i.merge]]] = mx.comb.pos
    }
    alignments[[df.merge$cl[i.merge] ]] = mx.comb.seq
    
    # ---- Release the two parents (see the same note in refineAlignment) ----
    ids.left = merge.rows[merge.rows > i.merge]
    ids.left = c(df.merge$id1[ids.left], df.merge$id2[ids.left])
    ids.drop = setdiff(c(i.cl1, i.cl2), ids.left)
    if(length(ids.drop) > 0){
      alignments[ids.drop] <- list(NULL)
      if(keep.pos) positions[ids.drop] <- list(NULL)
    }
    
    # if(i.merge == 21) stop()
    # 
    # msaplot(mx.comb.seq)
    # msadiff(mx.comb.seq)
    # 
    
    seq.cons.comb = mx2cons(mx.comb.seq)
    
    seqs.cl[paste0('clust_', df.merge$cl[i.merge])] = nt2seq(seq.cons.comb)
    
  }
  
  # # ---- Put N back ----
  # if(length(seqs.names.n) > 0){
  #   for(i.p in 1:length(positions)){
  #     pos = positions[[i.p]]
  #     for(s.name in seqs.names.n){
  #       if(s.name %in% row.names(pos)){
  #         pos.s = pos[s.name,]
  #         pos.s[pos.s != 0] = pos.replace[[s.name]]
  #         # stop()
  #       }
  #     }
  #   }
  # }
  
  return(list(pos = if(keep.pos) positions else NULL,
              aln = alignments))
}


#' Routed refinement. The old MAFFT --merge escalation used the width of the fast
#' alignment ("spread") as the signal, but a wide alignment is usually correct: two
#' different alleles at the same breakpoint belong in their own columns. It is off by
#' default; spread.max (or env PANNAGRAM_SPREAD_MAX) brings it back.
refineAlignmentRouted <- function(seqs.clean, path.work, spread.max = NA, keep.pos = FALSE){
  if(is.na(spread.max)){
    e = Sys.getenv("PANNAGRAM_SPREAD_MAX")
    spread.max = if(nzchar(e)) as.numeric(e) else Inf
  }
  R = refineAlignment(seqs.clean, path.work, keep.pos = keep.pos)
  if(!is.finite(spread.max)) return(R)
  B = R$aln[[length(R$aln)]]
  spread = ncol(B) / max(nchar(seqs.clean))
  if(spread <= spread.max) return(R)
  refineAlignment_prev(seqs.clean, path.work, keep.pos = keep.pos)   # opt-in: old MAFFT merge
}


# ============================================================================
#  FAST refinement (replaces the O(L^2) full MAFFT pairwise merge)
# ============================================================================

#' Align two long sequences via BLAST anchors + local MAFFT of short gaps.
#'
#' Returns a 2-row position matrix pos.mx: pos.mx[1,col] = position in s1 (or 0),
#' pos.mx[2,col] = position in s2 (or 0). Guarantees every position of s1 and s2 is
#' present exactly once and in ascending order (all positions kept, none dropped).
#'
#' A single MAFFT is tried first and accepted only if the identity over the co-aligned
#' bases is high enough; MAFFT is global and otherwise smears a non-homologous allele
#' along its partner. The threshold is stricter for longer stretches.
#'
#' Otherwise the backbone is a collinear chain of BLAST HSPs (chainHsps), the stretches
#' between anchors go to MAFFT under the same test, and a stretch too long for MAFFT is
#' re-anchored inside, up to max.depth levels. Only a stretch with no anchors at all is
#' kept SEPARATED, each side in its own columns.
alignTwoLong <- function(s1, s2, path.work, gap.mafft.max = 8000,
                         ident.min = 0.7, ident.long = 0.85, len.strict = 3000,
                         short.skip = 20, max.depth = 3){
  s1 = toupper(s1); s2 = toupper(s2); n1 = nchar(s1); n2 = nchar(s2)
  nt1 = seq2nt(s1); nt2 = seq2nt(s2)

  # identity a MAFFT-aligned stretch must reach to be accepted
  need.ident <- function(len) if(len > len.strict) ident.long else ident.min

  mafftPairMx <- function(a, b){
    writeFasta(setNames(c(a, b), c('a', 'b')), paste0(path.work, 'seg.fasta'))
    system(paste0('mafft --quiet --op 3 --ep 0.1 ', path.work, 'seg.fasta > ', path.work, 'seg_a.fasta'))
    aln2mx(readFasta(paste0(path.work, 'seg_a.fasta')))[c('a', 'b'), , drop = F]
  }
  mxToPos <- function(mx){
    c1 = rep(0, ncol(mx)); c2 = rep(0, ncol(mx))
    c1[mx[1, ] != '-'] = 1:n1; c2[mx[2, ] != '-'] = 1:n2; rbind(c1, c2)
  }
  # invariant: every position of s1 and s2 present exactly once, in ascending order
  valid <- function(pm){
    if(is.null(pm) || nrow(pm) != 2) return(FALSE)
    a = pm[1, pm[1, ] != 0]; b = pm[2, pm[2, ] != 0]
    length(a) == n1 && !anyDuplicated(a) && !is.unsorted(a) &&
    length(b) == n2 && !anyDuplicated(b) && !is.unsorted(b)
  }
  # identity over the columns where BOTH rows have a base
  identity2 <- function(mx){
    both = (mx[1, ] != '-') & (mx[2, ] != '-')
    if(sum(both) == 0) return(0)
    mean(mx[1, both] == mx[2, both])
  }

  # ---- Collinear chain of HSPs: maximum-weight chaining, O(n log n) ----
  # Greedy picking (longest first) drops most of the homology on long repeat-rich loci.
  # DP over the HSPs sorted by start, with a Fenwick tree (prefix maximum) over the
  # coordinate of the second sequence.
  chainHsps <- function(x){
    x = x[nrow(x) > 0 & (x$V2 < x$V3) & (x$V4 < x$V5), , drop = F]   # forward on BOTH
    if(nrow(x) == 0) return(x)
    x = x[order(x$V2, x$V4), , drop = F]
    n = nrow(x)
    w = x$V7                                  # weight = aligned length of the HSP
    keys = sort(unique(c(x$V4, x$V5)))        # compressed coordinates of sequence 2
    idx.end = match(x$V5, keys)
    idx.qry = match(x$V4, keys) - 1L          # prefix strictly before the HSP start
    m = length(keys)
    tree.val = rep(-Inf, m); tree.arg = rep(0L, m)
    fen.update <- function(i, val, arg){
      while(i <= m){
        if(val > tree.val[i]){ tree.val[i] <<- val; tree.arg[i] <<- arg }
        i = i + bitwAnd(i, -i)
      }
    }
    fen.query <- function(i){
      best = -Inf; arg = 0L
      while(i > 0){
        if(tree.val[i] > best){ best = tree.val[i]; arg = tree.arg[i] }
        i = i - bitwAnd(i, -i)
      }
      c(best, arg)
    }
    score = numeric(n); parent = integer(n)
    ord.end = order(x$V3); ptr = 1L
    for(i in seq_len(n)){
      # everything ending before this HSP starts may precede it
      while(ptr <= n && x$V3[ord.end[ptr]] < x$V2[i]){
        j = ord.end[ptr]; fen.update(idx.end[j], score[j], j); ptr = ptr + 1L
      }
      q = if(idx.qry[i] >= 1) fen.query(idx.qry[i]) else c(-Inf, 0)
      if(is.finite(q[1])){ score[i] = q[1] + w[i]; parent[i] = as.integer(q[2]) }
      else { score[i] = w[i]; parent[i] = 0L }
    }
    i = which.max(score); keep = c()
    while(i > 0){ keep = c(i, keep); i = parent[i] }
    x[keep, , drop = F]
  }

  cols1 = integer(0); cols2 = integer(0)
  em <- function(a, b){ cols1 <<- c(cols1, a); cols2 <<- c(cols2, b) }
  emitBlock <- function(q, s, off1, off2){        # BLAST-aligned strings of an HSP
    qc = seq2nt(q); sc = seq2nt(s); p1 = off1; p2 = off2
    c1 = integer(length(qc)); c2 = integer(length(qc))
    for(j in seq_along(qc)){
      if(qc[j] != '-'){ c1[j] = p1; p1 = p1 + 1 }
      if(sc[j] != '-'){ c2[j] = p2; p2 = p2 + 1 }
    }
    em(c1, c2)
  }
  separate <- function(g1a, g1b, g2a, g2b){       # no homology: own columns for each
    em(g1a:g1b, rep(0, g1b - g1a + 1)); em(rep(0, g2b - g2a + 1), g2a:g2b)
  }

  # Align the stretch s1[g1a..g1b] against s2[g2a..g2b].
  fillGap <- function(g1a, g1b, g2a, g2b, depth){
    l1 = g1b - g1a + 1; l2 = g2b - g2a + 1
    if(l1 <= 0 && l2 <= 0) return(invisible())
    if(l1 <= 0){ em(rep(0, l2), g2a:g2b); return(invisible()) }
    if(l2 <= 0){ em(g1a:g1b, rep(0, l1)); return(invisible()) }

    if(max(l1, l2) <= gap.mafft.max){
      mx = mafftPairMx(nt2seq(nt1[g1a:g1b]), nt2seq(nt2[g2a:g2b]))
      # forcing non-homologous stretches together is worse than leaving them apart
      if(max(l1, l2) >= short.skip && identity2(mx) < need.ident(max(l1, l2))){
        separate(g1a, g1b, g2a, g2b); return(invisible())
      }
      c1 = rep(0, ncol(mx)); c2 = rep(0, ncol(mx))
      c1[mx[1, ] != '-'] = g1a:g1b; c2[mx[2, ] != '-'] = g2a:g2b; em(c1, c2)
      return(invisible())
    }

    # too long for MAFFT: re-anchor inside instead of giving up
    if(depth < max.depth){
      a = tryCatch(chainHsps(blastTwoSeqs(nt2seq(nt1[g1a:g1b]), nt2seq(nt2[g2a:g2b]), path.work)),
                   error = function(e) NULL)
      if(!is.null(a) && nrow(a) > 0){
        cur1 = 1; cur2 = 1
        for(i in 1:nrow(a)){
          fillGap(g1a + cur1 - 1, g1a + a$V2[i] - 2, g2a + cur2 - 1, g2a + a$V4[i] - 2, depth + 1)
          emitBlock(a$V8[i], a$V9[i], g1a + a$V2[i] - 1, g2a + a$V4[i] - 1)
          cur1 = a$V3[i] + 1; cur2 = a$V5[i] + 1
        }
        fillGap(g1a + cur1 - 1, g1b, g2a + cur2 - 1, g2b, depth + 1)
        return(invisible())
      }
    }
    separate(g1a, g1b, g2a, g2b)
  }

  viaAnchors <- function(){
    cols1 <<- integer(0); cols2 <<- integer(0)
    a = chainHsps(blastTwoSeqs(s1, s2, path.work))
    cur1 = 1; cur2 = 1
    if(nrow(a) > 0) for(i in 1:nrow(a)){
      fillGap(cur1, a$V2[i] - 1, cur2, a$V4[i] - 1, 1)
      emitBlock(a$V8[i], a$V9[i], a$V2[i], a$V4[i])
      cur1 = a$V3[i] + 1; cur2 = a$V5[i] + 1
    }
    fillGap(cur1, n1, cur2, n2, 1)
    rbind(cols1, cols2)
  }

  # ---- a single MAFFT, if it is not forced ----
  if(max(n1, n2) <= gap.mafft.max){
    mx = tryCatch(mafftPairMx(s1, s2), error = function(e) NULL)
    if(!is.null(mx) && identity2(mx) >= need.ident(max(n1, n2))){
      pm = tryCatch(mxToPos(mx), error = function(e) NULL)
      if(valid(pm)) return(pm)
    }
  }

  pm = tryCatch(viaAnchors(), error = function(e) NULL)
  if(valid(pm)) return(pm)

  # ---- fallback: full MAFFT of the pair (rare edge cases; keeps every position) ----
  mxToPos(mafftPairMx(s1, s2))
}


#' Fast drop-in for refineAlignment: same clustering + progressive combine, but the
#' expensive pairwise merge (old: full MAFFT --merge + BLAST synteny + heavy R per
#' node) is replaced by alignTwoLong() (BLAST anchors + short-gap MAFFT). Output
#' structure is identical: list(pos, aln) with aln[[last]] the final MSA matrix.
refineAlignment <- function(seqs.clean, path.work, gap.mafft.max = 8000, keep.pos = FALSE){

  n.seqs = length(seqs.clean)
  seqs.clean = toupper(seqs.clean)

  # ---- Cluster ----
  dist.mx = calcDistKmer(seqs.clean)
  hc = hclust(as.dist(dist.mx))
  clusters <- cutree(hc, h = 0.1)

  # ---- Align each cluster (MAFFT on small groups) ----
  seqs.cl = c(); positions = list(); alignments = list()
  for(i.cl in 1:max(clusters)){
    seqs.tmp = seqs.clean[names(clusters)[clusters == i.cl]]
    if(length(seqs.tmp) == 1){
      seqs.cl = c(seqs.cl, seqs.tmp)
      if(keep.pos) positions[[i.cl]] = matrix(1:nchar(seqs.tmp), nrow = 1, dimnames = list(names(seqs.tmp), NULL))
      alignments[[i.cl]] = matrix(seq2nt(seqs.tmp),  nrow = 1, dimnames = list(names(seqs.tmp), NULL))
      next
    }
    seqs.cl.fasta = paste0(path.work, 'seqs_', i.cl, '.fasta')
    aln.fasta     = paste0(path.work, 'aln_',  i.cl, '.fasta')
    writeFasta(seqs.tmp, seqs.cl.fasta)
    system(paste('mafft  --quiet --maxiterate 100 ', seqs.cl.fasta, '>', aln.fasta, sep = ' '))
    seqs.cl.mx = aln2mx(readFasta(aln.fasta))
    if(keep.pos){
      pos.cl.mx = matrix(0L, nrow = nrow(seqs.cl.mx), ncol = ncol(seqs.cl.mx))
      for(irow in 1:nrow(seqs.cl.mx)) pos.cl.mx[irow, seqs.cl.mx[irow, ] != '-'] = 1:nchar(seqs.tmp[irow])
      rownames(pos.cl.mx) = names(seqs.tmp)
      positions[[i.cl]] = pos.cl.mx
    }
    alignments[[i.cl]] = seqs.cl.mx
    seqs.cl = c(seqs.cl, nt2seq(mx2cons(seqs.cl.mx)))
  }
  names(seqs.cl) = paste0('clust_', 1:length(seqs.cl))

  # ---- Hierarchical merge order (same as before) ----
  df.merge = as.data.frame(hc$merge)
  df.merge$id1 = ifelse(hc$merge[,1] < 0, clusters[abs(hc$merge[,1])], NA)
  df.merge$id2 = ifelse(hc$merge[,2] < 0, clusters[abs(hc$merge[,2])], NA)
  df.merge$cl = rep(0, nrow(hc$merge)); n.cl = max(clusters) + 1
  for(i in 1:nrow(df.merge)){
    df.merge$id1[i] = ifelse(is.na(df.merge$id1[i]), df.merge$cl[abs(hc$merge[i,1])], df.merge$id1[i])
    df.merge$id2[i] = ifelse(is.na(df.merge$id2[i]), df.merge$cl[abs(hc$merge[i,2])], df.merge$id2[i])
    df.merge$cl[i]  = ifelse(df.merge$id1[i] == df.merge$id2[i], df.merge$id1[i], n.cl)
    if(df.merge$cl[i] == n.cl) n.cl = n.cl + 1
  }

  # ---- Merge clusters via alignTwoLong (fast) ----
  merge.rows = which(df.merge$cl > max(clusters))
  for(i.merge in merge.rows){
    i.cl1 = df.merge$id1[i.merge]; i.cl2 = df.merge$id2[i.merge]
    s1 = seqs.cl[i.cl1]; s2 = seqs.cl[i.cl2]

    mx.comb = alignTwoLong(s1, s2, path.work, gap.mafft.max = gap.mafft.max)  # 2-row correspondence (all pos kept)

    non.zero.indices.1 = mx.comb[1,] != 0
    non.zero.indices.2 = mx.comb[2,] != 0
    n1 = nrow(alignments[[i.cl1]]); n2 = nrow(alignments[[i.cl2]])

    if(keep.pos){
      mx.comb.pos = matrix(0L, nrow = n1 + n2, ncol = ncol(mx.comb))
      mx.comb.pos[1:n1,        non.zero.indices.1] = positions[[i.cl1]][, mx.comb[1, non.zero.indices.1]]
      mx.comb.pos[n1 + (1:n2), non.zero.indices.2] = positions[[i.cl2]][, mx.comb[2, non.zero.indices.2]]
    }

    mx.comb.seq = matrix('-', nrow = n1 + n2, ncol = ncol(mx.comb))
    mx.comb.seq[1:n1,        non.zero.indices.1] = alignments[[i.cl1]][, mx.comb[1, non.zero.indices.1]]
    mx.comb.seq[n1 + (1:n2), non.zero.indices.2] = alignments[[i.cl2]][, mx.comb[2, non.zero.indices.2]]

    tmp.names = c(rownames(alignments[[i.cl1]]), rownames(alignments[[i.cl2]]))
    rownames(mx.comb.seq) = tmp.names
    if(keep.pos){
      rownames(mx.comb.pos) = tmp.names
      positions[[df.merge$cl[i.merge]]] = mx.comb.pos
    }
    alignments[[df.merge$cl[i.merge]]] = mx.comb.seq

    # ---- Release the two parents ----
    # Each node of the merge tree is consumed by exactly one merge, so once the child
    # matrix exists its parents are dead. Without this the list retained EVERY
    # intermediate matrix until the function returned: peak = sum over all tree nodes
    # (~n.seqs x aln.len x 8 bytes per level x tree depth) instead of the live frontier.
    # `x[i] <- list(NULL)` empties the slot but keeps the indexing, so aln[[last]] is
    # still the root alignment.
    ids.left = merge.rows[merge.rows > i.merge]
    ids.left = c(df.merge$id1[ids.left], df.merge$id2[ids.left])
    ids.drop = setdiff(c(i.cl1, i.cl2), ids.left)
    if(length(ids.drop) > 0){
      alignments[ids.drop] <- list(NULL)
      if(keep.pos) positions[ids.drop] <- list(NULL)
    }

    seqs.cl[paste0('clust_', df.merge$cl[i.merge])] = nt2seq(mx2cons(mx.comb.seq))
  }

  return(list(pos = if(keep.pos) positions else NULL, aln = alignments))
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
mafftAdd <- function(seq1, seq2, path.work, n.diff = 5, n.flank = 0){
  
  # Create a tablefile (so-called subMSAtable) for mafft 
  id1 = 1:length(seq1)
  id2 = (length(seq1) + 1):(length(seq1) + length(seq2))
  
  first_line <- paste(id1, collapse = " ")
  second_line <- paste(id2, collapse = " ")
  file.table = paste0(path.work, "tbl.txt")
  writeLines(c(first_line, second_line, ''), file.table)
  
  # Merge sequences into the input file
  file.seqs.merge = paste0(path.work, 'seqs_tmp_merge.fasta')
  
  if(n.flank == 0){
    writeFasta(c(seq1, seq2), file.seqs.merge)  
  } else {
    s.flank.beg = nt2seq(rep('A', n.flank))
    s.flank.end = nt2seq(rep('T', n.flank))
    
    seqs = c(seq1, seq2)
    for(i in 1:length(seqs)){
      seqs[i] = paste0(s.flank.beg, seqs[i], s.flank.end)
    }
    
    writeFasta(seqs, file.seqs.merge)  
  }
  
  # File for storing MAFFT alignment output
  file.mafft = paste0(path.work, 'seqs_tmp_mafft.fasta')
  
  # Run MAFFT with the specified table and merged sequences
  system(paste('mafft --op 3  --ep 0.1 --quiet --merge', file.table, file.seqs.merge, '>', file.mafft, sep = ' '))
  
  # Read the aligned sequences
  xx = readFasta(file.mafft)
  mafft.mx = aln2mx(xx)
  
  if(n.flank != 0){
    for(i in 1:nrow(mafft.mx)){
      idx.pos = which(mafft.mx[i,] != '-')
      
      n.first <- idx.pos[1:n.flank]
      n.last <- idx.pos[(length(idx.pos) - n.flank + 1):length(idx.pos)]
      mafft.mx[i,n.first] = '-'
      mafft.mx[i,n.last] = '-'
      
    }
    mafft.mx = mafft.mx[,colSums(mafft.mx != '-') > 0, drop = F]
  }
  
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
  
  if(is.null(names(s1))){
    names(s1) = 'xx'
  }
  
  if(is.null(names(s2))){
    names(s2) = 'yy'
  }
  writeFasta(c(s1, s2), file.seqs.tmp)
  
  # Define the output file for BLAST results
  file.blast.cons = paste0(file.seqs.tmp, '.out')
  
  # Create a BLAST database from the temporary sequence file
  system(paste0('makeblastdb -in ', file.seqs.tmp,' -dbtype nucl  > /dev/null'))
  
  # Run BLASTn alignment between the sequences in the temporary file
  system(paste0('blastn -db ',file.seqs.tmp,' -query ',file.seqs.tmp,
                ' -num_alignments 50 ',
                ' -out ',file.blast.cons,
                ' -outfmt "6 qseqid qstart qend sstart send pident length qseq sseq sseqid"'))
  
  # Read the BLAST output into a data frame
  x = readBlast(file.blast.cons)
  if(is.null(x)) {
    return(data.frame(tmp=numeric()))
  }
  x = x[x$V1 != x$V10,,drop=F] # Remove self-alignments
  
  # Split data into two groups: alignments starting from sequence 1 and sequence 2
  idx1 = x$V1 %in% names(s1)
  x1 = x[idx1,,drop=F]
  x2 = x[!idx1,,drop=F]
  
  # Swap columns for the second group to match the first group's format
  colnames(x2)[2:5] = c('V4', 'V5', 'V2', 'V3') 
  colnames(x2)[c(1,10)] = colnames(x2)[c(10,1)]
  
  colnames(x2)[8:9] = c('V9', 'V8') 
  
  
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
  writeFasta(s1, file.seqs.tmp1)
  
  file.seqs.tmp2 = paste0(path.work, 'seqs_tmp2.fasta')
  writeFasta(s2, file.seqs.tmp2)
  
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
                ' -outfmt "6 qseqid qstart qend sstart send pident length sseqid"'))
  
  system(paste0('blastn -db ',file.seqs.tmp1,' -query ',file.seqs.tmp2,
                ' -num_alignments 50 ',
                ' -out ',file.blast2,
                ' -outfmt "6 qseqid qstart qend sstart send pident length sseqid"'))
  
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
calcDistAln <- function(seqs.mx) {
  
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


calcDistKmer <- function(seqs.clean){

  wsize = 7
  seqs.clean <- toupper(seqs.clean)
  n <- length(seqs.clean)

  # Count only the OBSERVED k-mers per sequence (sparse). The old code built the
  # full 4^wsize vocabulary (expand.grid + apply(paste)) and a dense matrix of that
  # many columns every call, then dropped the empties -- pure overhead. The Manhattan
  # distance is identical: never-observed k-mers contribute 0 to every pairwise sum.
  tabs <- lapply(seqs.clean, function(seq){
    kmers <- substring(seq, 1:(nchar(seq) - wsize + 1), wsize:nchar(seq))
    table(kmers)
  })
  all.k <- unique(unlist(lapply(tabs, names), use.names = FALSE))
  M <- matrix(0, nrow = n, ncol = length(all.k), dimnames = list(NULL, all.k))
  for (i in 1:n) M[i, names(tabs[[i]])] <- as.numeric(tabs[[i]])

  dist.mx <- as.matrix(dist(M, method = "manhattan"))

  vec <- nchar(seqs.clean)
  dist.mx <- dist.mx / outer(vec, vec, pmax)

  dimnames(dist.mx) <- list(names(seqs.clean), names(seqs.clean))
  return(dist.mx)
}


clustDistKmer <- function(dist.mx, sim.cutoff = 0.1){
  
  n = nrow(dist.mx)
  dist.names = rownames(dist.mx)
  if(is.null(names)){
    pokazAttention('Names of rows in distance matrix were assigned')
    dist.names = paste0('s', 1:n)
  }
  indices <- which(dist.mx <= sim.cutoff, arr.ind = TRUE)
  if(nrow(indices) == 0){
    return()
  }
  
  # # Define classes
  # classes = 1:n
  # for(irow in 1:nrow(indices)){
  #   cl1 = classes[indices[irow,1]]
  #   cl2 = classes[indices[irow,2]]
  #   if(cl1 != cl2){
  #     classes[classes == cl2] = cl1
  #   }
  # }
  # classes = as.numeric(as.factor(classes))
  
  graph <- igraph::graph_from_data_frame(indices, directed = FALSE)
  comp <- components(graph)
  
  classes = comp$membership = comp$membership[order(as.numeric(names(comp$membership)))]
  names(classes) = dist.names
  
  return(classes)
  
}


# 
# 
# xx = seqs.mx[as.numeric(names(comp$membership)),]
# 
# msadiff(xx)
# 
# msadiff(seqs.mx[as.numeric(names(comp$membership)),])
# 
# msadiff(seqs.mx)
# 
# end <- Sys.time()
# 
# execution_time <- end - start
# print(execution_time)


