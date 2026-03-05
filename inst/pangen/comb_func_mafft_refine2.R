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
refineAlignment <- function(seqs.clean, path.work){
  
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
  
  # ---- Merge clusters ----
  
  for(i.merge in which(df.merge$cl > max(clusters))){
    
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
    result$support = 0
    
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
    for(irow in 1:nrow(result)){
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
  
  nts <- c('A', 'C', 'G', 'T')
  combinations <- expand.grid(rep(list(nts), wsize))
  combinations = apply(combinations, 1, paste0, collapse = '')
  
  start <- Sys.time()
  
  df <- setNames(data.frame(matrix(0, nrow = length(seqs.clean), ncol = length(combinations))), combinations)
  
  for (i in 1:length(seqs.clean)) {
    seq = seqs.clean[i]
    kmers <- substring(seq, 1:(nchar(seq) - wsize + 1), wsize:nchar(seq))
    k.mx.cnt <- table(kmers)
    df[i, names(k.mx.cnt)] <- as.numeric(k.mx.cnt)
  }
  df = df[,combinations]
  df <- df[, colSums(df) > 0]
  
  dist.mx <- as.matrix(dist(df, method = "manhattan"))
  
  vec = nchar(seqs.clean)
  max_matrix <- outer(vec, vec, pmax)
  
  dist.mx = dist.mx / max_matrix
  
  colnames(dist.mx) = names(seqs.clean)
  rownames(dist.mx) = names(seqs.clean)
  
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


