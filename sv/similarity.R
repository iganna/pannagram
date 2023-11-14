# All function to find similarities

"Length (len1) should be defined before"
findHitsInRef <- function(v, sim.cutoff = 0.9){
  
  s.tmp.comb = '___'
  
  idx.include = (v$V7 / v$len1 > sim.cutoff)
  v.include = v[idx.include,]
  idx.non.include = !idx.include
  
  v.sim = v.include
  
  message('Work with partial genes')
  idx.strand = (v$V4 > v$V5) * 1
  for(i.strand in 0:1){
    message(paste('Strand', i.strand))
    v.rest = v[(idx.non.include) & (idx.strand == i.strand),]
    if(i.strand == 1){
      tmp = v.rest$V4
      v.rest$V4 = v.rest$V5
      v.rest$V5 = tmp
    }
    
    # if only one record - delete
    v.rest = v.rest[order(v.rest$V8),]
    v.rest = v.rest[order(v.rest$V1),]
    n.rest = nrow(v.rest)
    idx.one = which((v.rest$V8[-1] == v.rest$V8[-n.rest]) & (v.rest$V1[-1] == v.rest$V1[-n.rest]))
    idx.one = sort(unique(c(idx.one,idx.one+1)))
    v.rest = v.rest[idx.one,]
    
    v.rest = v.rest[order(-v.rest$V5),]
    v.rest = v.rest[order(v.rest$V4),]
    v.rest = v.rest[order(v.rest$V8),]
    v.rest = v.rest[order(v.rest$V1),]
    # remove nestedness
    idx.nested = 1
    while(length(idx.nested) > 0){
      idx.nested = which((v.rest$V2[-1] >= v.rest$V2[-nrow(v.rest)]) & 
                           (v.rest$V3[-1] <=v.rest$V3[-nrow(v.rest)]) & 
                           (v.rest$V4[-1] >= v.rest$V4[-nrow(v.rest)]) & 
                           (v.rest$V5[-1] <=v.rest$V5[-nrow(v.rest)]) & 
                           (v.rest$V1[-1] == v.rest$V1[-nrow(v.rest)]) &
                           (v.rest$V8[-1] == v.rest$V8[-nrow(v.rest)])) + 1
      # print(length(idx.nested))
      if(length(idx.nested) == 0) next
      v.rest = v.rest[-idx.nested,]  
    }
    
    v.rest$cover = v.rest$V3 - v.rest$V2 + 1
    v.rest$ref.overlap1 = c(v.rest$V4[-1] - v.rest$V5[-nrow(v.rest)] - 1, 0)
    v.rest$allowedoverlap1 = v.rest$len1 * (1-sim.cutoff)
    suffixname = (v.rest$ref.overlap1 > v.rest$allowedoverlap1) * 1
    v.rest$suffixname = c(1, suffixname[-length(suffixname)])
    v.rest$suffixname[1 + which(v.rest$V8[-1] != v.rest$V8[-nrow(v.rest)])] = 1
    v.rest$suffixname = cumsum(v.rest$suffixname)
    v.rest$V8 = paste(v.rest$V8, v.rest$suffixname, sep = '|id')
    
    # if only one record - delete
    v.rest = v.rest[order(v.rest$V8),]
    v.rest = v.rest[order(v.rest$V1),]
    n.rest = nrow(v.rest)
    idx.one = which((v.rest$V8[-1] == v.rest$V8[-n.rest]) & (v.rest$V1[-1] == v.rest$V1[-n.rest]))
    idx.one = sort(unique(c(idx.one,idx.one+1)))
    v.rest = v.rest[idx.one,]
    
    v.rest = v.rest[order(-v.rest$V3),]
    v.rest = v.rest[order(v.rest$V2),]
    v.rest = v.rest[order(v.rest$V8),]
    v.rest = v.rest[order(v.rest$V1),]
    # remove nestedness
    idx.nested = 1
    while(length(idx.nested) > 0){
      idx.nested = which((v.rest$V2[-1] >= v.rest$V2[-nrow(v.rest)]) & 
                           (v.rest$V3[-1] <=v.rest$V3[-nrow(v.rest)]) & 
                           (v.rest$V1[-1] == v.rest$V1[-nrow(v.rest)]) &
                           (v.rest$V8[-1] == v.rest$V8[-nrow(v.rest)])) + 1
      # print(length(idx.nested))
      if(length(idx.nested) == 0) next
      v.rest = v.rest[-idx.nested,]  
    }
    
    v.rest$cover = v.rest$V3 - v.rest$V2 + 1
    v.rest$overlap1 = c(v.rest$V2[-1] - v.rest$V3[-nrow(v.rest)] - 1, 0)
    v.rest$overlap1[v.rest$overlap1 > 0] = 0
    idx.diff = which(v.rest$V8[-1] != v.rest$V8[-nrow(v.rest)])
    v.rest$overlap1[idx.diff] = 0
    idx.diff = which(v.rest$V1[-1] != v.rest$V1[-nrow(v.rest)])
    v.rest$overlap1[idx.diff] = 0
    v.rest$cover = v.rest$cover + v.rest$overlap1
    
    v.rest$comb = paste(v.rest$V1, v.rest$V8, sep = s.tmp.comb)
    df.cover = data.frame(V1 = tapply(v.rest$V1, v.rest$comb, unique),
                          V2 = tapply(v.rest$V2, v.rest$comb, min),
                          V3 = tapply(v.rest$V3, v.rest$comb, max),
                          V4 = tapply(v.rest$V4, v.rest$comb, min),
                          V5 = tapply(v.rest$V5, v.rest$comb, max),
                          V6 = tapply(v.rest$V6, v.rest$comb, mean),
                          V7 = tapply(v.rest$cover, v.rest$comb, sum), 
                          V8 = tapply(v.rest$V8, v.rest$comb, unique),
                          # comb = tapply(v.rest$comb, v.rest$comb, unique),
                          len1 = tapply(v.rest$len1, v.rest$comb, unique))
    # df.cover$dir = i.strand
    df.cover$ref.cover = df.cover$V5 - df.cover$V4 + 1
    rownames(df.cover) = NULL
    if(i.strand == 1){
      tmp = df.cover$V4
      df.cover$V4 = df.cover$V5
      df.cover$V5 = tmp
    }
    df.cover$V8 = sapply(df.cover$V8, function(s) strsplit(s, '\\|')[[1]][1])
    
    idx.include = (df.cover$V7 / df.cover$len1 > sim.cutoff) & (df.cover$V6 > sim.cutoff * 100) & 
      (df.cover$ref.cover / df.cover$len1 > sim.cutoff) & (df.cover$len1 / df.cover$ref.cover > sim.cutoff)
    v.sim = rbind(v.sim, df.cover[idx.include, colnames(v.sim)])
  }
  return(v.sim)
}

findMutualSimilarity_old <- function(v, sim.cutoff = 0.9, pos.len1 = NULL, pos.len2 = NULL){
  
  # sim.cutoff = 0.9
  # pos.len1 = 2
  # pos.len2 = 5
  
  s.tmp.comb = '___'
  
  v.strand = list(v[v$V4 < v$V5,], v[v$V4 > v$V5,])
  v.all = c()
  gc()
  # Remove nestedness
  for(i.strand in 1:2){
    message(paste('Strand', i.strand))
    v = v.strand[[i.strand]]
    if(i.strand == 2){
      tmp = v$V4
      v$V4 = v$V5
      v$V5 = tmp
    }
    #Remove duplicates
    v = v[order(-v$V5),]
    v = v[order(v$V4),]
    v = v[order(-v$V3),]
    v = v[order(v$V2),]
    v = v[order(v$V8),]
    v = v[order(v$V1),]
    n = 0
    while(n != nrow(v)){  # another strand strand
      n = nrow(v)
      print(n)
      idx = which((v$V1[-nrow(v)] == v$V1[-1]) & 
                    (v$V2[-nrow(v)] <= v$V2[-1]) & 
                    (v$V3[-nrow(v)] >= v$V3[-1]) &
                    (v$V4[-nrow(v)] <= v$V4[-1]) & 
                    (v$V5[-nrow(v)] >= v$V5[-1]) &
                    (v$V8[-nrow(v)] == v$V8[-1])) + 1
      if(length(idx) == 0) next
      v = v[-idx,]
    }
    if(i.strand == 2){
      tmp = v$V4
      v$V4 = v$V5
      v$V5 = tmp
    }
    v.all = rbind(v.all, v)
    print(nrow(v.all))
  }
  v = v.all
  # rm(v.all)
  gc()
  
  
  # Define lengths - > to a function
  if(!('len1' %in% colnames(v)) && (is.null(pos.len1))) stop('Information about length is not provided')
  if(!('len2' %in% colnames(v)) && (is.null(pos.len2))) stop('Information about length is not provided')
  if(!('len1' %in% colnames(v))){
    v$len1 = as.numeric(sapply(v$V1, function(s) strsplit(s,'\\|')[[1]][pos.len1]))
  }
  if(!('len2' %in% colnames(v))){
    v$len2 = as.numeric(sapply(v$V8, function(s) strsplit(s,'\\|')[[1]][pos.len2]))
  }
  
  idx.include = (v$V7 / v$len1 > sim.cutoff) & (v$V7 / v$len2 > sim.cutoff)
  v.sim = v[idx.include,]
  
  # - - - - - - - - - - - - - - - - - - - - - - - - 
  # Similarity one direction
  # Coverage on query
  v.rest = v[!idx.include,]
  
  # if only one record - delete
  v.rest$comb = paste(v.rest$V1, v.rest$V8, sep = s.tmp.comb)
  v.rest = v.rest[order(-v.rest$V3),]
  v.rest = v.rest[order(v.rest$V2),]
  v.rest = v.rest[order(v.rest$comb),]
  idx.one = which(v.rest$comb[-1] == v.rest$comb[-nrow(v.rest)])
  idx.one = sort(unique(c(idx.one,idx.one+1)))
  v.rest = v.rest[idx.one,]
  v.rest.saved = v.rest
  
  for(i.mode in 1:2){
    
    v.rest = v.rest[v.rest$len1 * sim.cutoff < v.rest$len2]
    # remove nestedness
    idx.nested = 1
    while(length(idx.nested) > 0){
      idx.nested = which((v.rest$V2[-1] >= v.rest$V2[-nrow(v.rest)]) & 
                           (v.rest$V3[-1] <=v.rest$V3[-nrow(v.rest)]) & 
                           (v.rest$comb[-1] == v.rest$comb[-nrow(v.rest)])) + 1
      # print(length(idx.nested))
      if(length(idx.nested) == 0) next
      v.rest = v.rest[-idx.nested,]  
    }
    
    # Coverage 
    v.rest$cover1 = v.rest$V3 - v.rest$V2 + 1
    v.rest$overlap1 = c(v.rest$V2[-1] - v.rest$V3[-nrow(v.rest)] - 1, 0)
    v.rest$overlap1[v.rest$overlap1 > 0] = 0
    idx.diff = which(v.rest$comb[-1] != v.rest$comb[-nrow(v.rest)])
    v.rest$overlap1[idx.diff] = 0
    v.rest$cover1 = v.rest$cover1 + v.rest$overlap1
    
    df.cover = data.frame(V7 = tapply(v.rest$cover1, v.rest$comb, sum),  # coverage
                          len1 = tapply(v.rest$len1, v.rest$comb, unique))
    df.cover = df.cover[df.cover$V7 >= df.cover$len1 * sim.cutoff,]
    
    if(i.mode == 1){
      comb.query = rownames(df.cover)
    } else {
      comb.ref = rownames(df.cover)
      break
    }
    
    # Prepare data for the next mode: oerlap on the reference
    v.rest = v.rest.saved
    # sort the reference coordinates
    idx = v.rest$V4 > v.rest$V5
    tmp = v.rest$V5[idx]
    v.rest$V5[idx] = v.rest$V4[idx]
    v.rest$V4[idx] = tmp
    # put the reference coordinates to the query place
    v.rest$V2 = v.rest$V4
    v.rest$V3 = v.rest$V5
    
    v.rest$V4 = NULL
    v.rest$V5 = NULL
    tmp = v.rest$len1
    v.rest$len1 = v.rest$len2
    v.rest$len2 = tmp
    
    v.rest = v.rest[order(-v.rest$V3),]
    v.rest = v.rest[order(v.rest$V2),]
    v.rest = v.rest[order(v.rest$comb),]
  }
  
  comb.mutual = intersect(comb.query, comb.ref)
  comb.query = setdiff(comb.query, comb.mutual)
  comb.ref = setdiff(comb.ref, comb.mutual)
  
  
  comb.mutual = matrix(unlist(strsplit(comb.mutual, s.tmp.comb)), ncol = 2, byrow = TRUE)
  comb.query = matrix(unlist(strsplit(comb.query, s.tmp.comb)), ncol = 2, byrow = TRUE)
  comb.ref = matrix(unlist(strsplit(comb.ref, s.tmp.comb)), ncol = 2, byrow = TRUE)
  
  comb.ref = comb.ref[,c(2,1)]
  
  
  comb.mutual = rbind(comb.mutual, as.matrix(v.sim[,c('V1', 'V8')]))
  
  return(list(mutual = comb.mutual, covered.query = comb.query, covered.ref = comb.ref))
}



getOneSideCoverage <- function(v.rest, side = 0){
  # it's assumes that everything is presorted before: V4 < V5
  
  # print(head(v.rest))
  if(side == 1){

    
    idx.tmp = v.rest$V4 > v.rest$V5
    if(sum(idx.tmp) > 0){
      tmp = v.rest$V4[idx.tmp]
      v.rest$V4[idx.tmp] = v.rest$V5[idx.tmp]
      v.rest$V5[idx.tmp] = tmp
    }
    
    v.rest$V1 = v.rest$V8
    v.rest$V2 = v.rest$V4
    v.rest$V3 = v.rest$V5
    
  }
  v.rest = v.rest[, 1:3]
  # print(head(v.rest))
  
  # - - - - - - - - - - - - - - - - - - - - - - - - 

  v.rest = v.rest[order(-v.rest$V3),]
  v.rest = v.rest[order(v.rest$V2),]
  v.rest$V1 = as.factor(v.rest$V1)
  v.rest = v.rest[order(v.rest$V1),]
  
  # remove nestedness
  idx.nested = 1
  while(length(idx.nested) > 0){
    idx.nested = which((v.rest$V2[-1] >= v.rest$V2[-nrow(v.rest)]) & 
                         (v.rest$V3[-1] <=v.rest$V3[-nrow(v.rest)]) & 
                         (v.rest$V1[-1] == v.rest$V1[-nrow(v.rest)])) + 1
    print(length(idx.nested))
    if(length(idx.nested) == 0) next
    v.rest = v.rest[-idx.nested,]  
  }
    
  v.rest$cover = v.rest$V3 - v.rest$V2 + 1
  v.rest$overlap = c(v.rest$V2[-1] - v.rest$V3[-nrow(v.rest)] - 1, 0)
  v.rest$overlap[v.rest$overlap > 0] = 0
  idx.diff = which(v.rest$V1[-1] != v.rest$V1[-nrow(v.rest)])
  v.rest$overlap[idx.diff] = 0
  v.rest$cover = v.rest$cover + v.rest$overlap
  
  coverage = tapply(v.rest$cover, v.rest$V1, sum)
  return(coverage)
}


findNestedness <- function(v.res, use.strand = T){
  
  idx.strand = v.res$V4 > v.res$V5
  tmp = v.res$V4[idx.strand]
  v.res$V4[idx.strand] = v.res$V5[idx.strand]
  v.res$V5[idx.strand] = tmp
  
  
  # Make unique names for strands
  if(use.strand == T){
    str.strand = c('+', '-')
    dir.strand = str.strand[idx.strand * 1 + 1]
    v.res$V8 = paste(v.res$V8, dir.strand, sep = '|')
  } 
  v.res$comb = paste(v.res$V1, v.res$V8, sep = '||')
  
  v.res$V1 = v.res$comb
  v.res$V8 = v.res$comb
  cover1 = getOneSideCoverage(v.res)
  cover8 = getOneSideCoverage(v.res, side = 1)
  
  
  v.cover = data.frame(C1 = cover1)
  v.cover$C8 = cover8[rownames(v.cover)]
  v.cover[,c('V1', 'V8')] = stringr::str_split_fixed(rownames(v.cover), "\\|\\|", 2)
  rownames(v.cover) = NULL
  
  if(use.strand){
    v.cover$dir = substr(v.cover$V8, nchar(v.cover$V8), nchar(v.cover$V8))
    v.cover$V8 = substr(v.cover$V8, 1, (nchar(v.cover$V8) - 2))
  } else {
    v.cover$dir = '.'
  }
  
  return(v.cover)
  # v.nest = v.cover[(v.cover$p1 >= sim.cutoff) | (v.cover$p8 >= sim.cutoff),]
  
}













