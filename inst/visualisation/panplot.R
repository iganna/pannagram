# This file contains functions for visualising synteny between genomes

#' @export
getBlocks <- function(v, f.split = T){
  
  len.min = 20000
  v.init = v
  v.idx = 1:length(v)
  
  v.idx = v.idx[v != 0]
  v = v[v != 0]
  # v.r = rank(abs(v))
  v.r = v
  v.r[v < 0] = v.r[v < 0] * (-1)
  v.b = findRuns(v.r)
  
  v.b = v.b[v.b$len >= len.min,]
  
  v.b.rank = matrix(rank(v.b[,c('v.beg', 'v.end')]), ncol = 2, dimnames = list(NULL,c('r.beg', 'r.end')))
  v.b = cbind(v.b, v.b.rank)
  
  
  # Remove 
  irow = 2
  while (irow <= nrow(v.b)) {
    if(((v.b$r.beg[irow] - 1) == v.b$r.end[irow - 1]) &
       ((v.b$v.beg[irow] - 1 - v.b$v.end[irow - 1]) < len.min)){
      v.b$r.end[irow - 1] = v.b$r.end[irow]
      v.b$end[irow - 1] = v.b$end[irow]
      v.b$v.end[irow - 1] = v.b$v.end[irow]
      v.b = v.b[-irow,]
    } else {
      irow = irow + 1
    }
  }
  v.b$i.beg = v.idx[v.b$beg]
  v.b$i.end = v.idx[v.b$end]
  
  df = v.b[,c('v.beg', 'v.end', 'i.beg', 'i.end')]
  
  # Add breaks, if they are with a long gap
  if(f.split){
    len.min.split = 50000
    d = diff(v.idx)
    i.d = which(d >= len.min.split)
    for(i.d in which(d >= len.min.split)){
      v.pos = v[i.d]
      i.pos = v.idx[i.d]
      
      i.v = which((df$v.beg < v.pos) & (df$v.end > v.pos))
      i.i = which((df$i.beg < i.pos) & (df$i.end > i.pos))
      
      if(length(i.v) > 0){
        if(length(i.i) == 0) stop('Something is wrong with i.i')
        if(i.v != i.i) stop('Idxs do not match')
        
        df.tmp = df[i.v,]
        
        df.tmp$v.beg = abs(v[i.d + 1])
        df.tmp$i.beg = v.idx[i.d + 1]
        
        df[i.v,]$v.end = v.pos
        df[i.v,]$i.end = i.pos
        
        df = rbind(df, df.tmp)
      }
    }    
  }
  

  
  df = df[order(df$i.beg),]
  rownames(df) = NULL
  
  colnames(df) <- c('own.b', 'own.e', 'pan.b', 'pan.e')
  df$dir = (df$own.b > df$own.e) * 1
  return(df)
}

extractBlocks <- function(path.aln, n.acc, len.break = 50000, len.break.gap = 10000, echo=T){
  
  idx.synteny = c()
  for(i.chr in 1:5){
    if(echo) message(paste('Chromosome', i.chr))
    v = readRDS(paste0(path.aln, 'val_common_chr_', i.chr, '_ref_add.rds'))
    n.v = nrow(v)
    accessions = colnames(v)[1:n.acc]
    
    n.acc.per.site = rowSums(v[1:n.acc] != 0, )
    is.core = (n.acc.per.site == n.acc) * 1
    
    zeros.beg <- which(is.core == 0 & c(0, is.core[-n.v]) != 0)
    zeros.end <- which(is.core == 0 & c(is.core[-1], 0) != 0)
    if(is.core[1] == 0) zeros.beg = c(1, zeros.beg)
    if(is.core[n.v] == 0) zeros.end = c(zeros.end, n.v)
    
    zeros.len = zeros.end - zeros.beg + 1
    idx.long = which(zeros.len >= len.break)
    
    pos.long = c()
    for(i in idx.long){
      pan.b = zeros.beg[i] - 1
      pan.e = zeros.end[i] + 1
      if(pan.e > n.v) next
      if(pan.b == 0 ) next
      for(acc in accessions){
        acc.b = v[pan.b, acc] - 1
        acc.e = v[pan.e, acc]
        if((acc.e - acc.b < len.break.gap)) next
        pos.long = rbind(pos.long, data.frame(pos = c(acc.b,acc.e), acc = acc))
      }
    }
    
    if(echo) cat('Accessions: ')
    # gaps per accession
    for(acc in accessions){
      if(echo) cat(acc)
      if(echo) cat(' ')
      
      v.acc = v[,acc]
      v.acc = cbind(v.acc, 1:nrow(v))
      
      v.acc = v.acc[v.acc[,1] != 0,] 
      v.acc = cbind(v.acc, 1:nrow(v.acc))
      v.acc = v.acc[order(v.acc[,1]),]
      
      diffs <- diff(v.acc[,3])
      idx.dir <- which(abs(diffs) != 1)
      
      diffs <- diff(v.acc[,1])
      idx.dir <- c(0, sort(c(idx.dir, which(abs(diffs) >= len.break))), nrow(v.acc))
      
      pos.acc = pos.long[pos.long[,2] == acc,1]
      
      idx.dir = sort(c(idx.dir, which(v.acc[,1] %in% pos.acc)))
      
      for(i in 1:(length(idx.dir)-1)){
        pan.b = v.acc[idx.dir[i]+1,2]
        pan.e = v.acc[idx.dir[i+1],2]
        acc.b = v.acc[idx.dir[i]+1,1]
        acc.e = v.acc[idx.dir[i+1],1]
        
        if((acc.e - acc.b + 1) < len.break.gap) next
        idx.synteny = rbind(idx.synteny,
                            data.frame(
                              pan.b = pan.b,
                              pan.e = pan.e,
                              own.b = acc.b,
                              own.e = acc.e,
                              acc = acc,
                              chr = i.chr,
                              # pan.len = pan.e - pan.b - 1,
                              # acc.len = acc.e - acc.b - 1,
                              dir = (v.acc[idx.dir[i+1],3] < v.acc[idx.dir[i]+1,3]) * 1
                            ))
      }
      rm(v.acc)
    }
    
    if(echo) cat('\n')
    
    vars_to_keep <- c("idx.synteny", "i.chr", 'len.break', 'len.break.gap', 'path.aln', 'n.acc')
    all_vars <- ls()
    vars_to_remove <- setdiff(all_vars, vars_to_keep)
    rm(list = vars_to_remove)
    gc()
  }
  
  colnames(idx.synteny)[c(3,4)] = c('own.b', 'own.e')
  return(idx.synteny)
}


# Function to extract blocks for every accession
extractBlocks_old <- function(path.aln, n.acc, len.break = 50000, len.break.gap = 10000, echo=T){
  
  idx.synteny = c()
  for(i.chr in 1:5){
    if(echo) message(paste('Chromosome', i.chr))
    v = readRDS(paste0(path.aln, 'val_common_chr_', i.chr, '_ref_add.rds'))
    n.v = nrow(v)
    accessions = colnames(v)[1:n.acc]
    
    n.acc.per.site = rowSums(v[1:n.acc] != 0, )
    is.core = (n.acc.per.site == n.acc) * 1
    
    zeros.beg <- which(is.core == 0 & c(0, is.core[-n.v]) != 0)
    zeros.end <- which(is.core == 0 & c(is.core[-1], 0) != 0)
    if(is.core[1] == 0) zeros.beg = c(1, zeros.beg)
    if(is.core[n.v] == 0) zeros.end = c(zeros.end, n.v)
    
    zeros.len = zeros.end - zeros.beg + 1
    idx.long = which(zeros.len >= len.break)
    
    pos.long = c()
    for(i in idx.long){
      pan.b = zeros.beg[i] - 1
      pan.e = zeros.end[i] + 1
      if(pan.e > n.v) next
      if(pan.b == 0 ) next
      for(acc in accessions){
        acc.b = v[pan.b, acc] - 1
        acc.e = v[pan.e, acc]
        if((acc.e - acc.b < len.break.gap)) next
        pos.long = rbind(pos.long, data.frame(pos = c(acc.b,acc.e), acc = acc))
      }
    }
    
    if(echo) cat('Accessions: ')
    # gaps per accession
    for(acc in accessions){
      if(echo) cat(acc)
      if(echo) cat(' ')
      
      v.acc = v[,acc]
      v.acc = cbind(v.acc, 1:nrow(v))
      
      v.acc = v.acc[v.acc[,1] != 0,] 
      v.acc = cbind(v.acc, 1:nrow(v.acc))
      v.acc = v.acc[order(v.acc[,1]),]
      
      diffs <- diff(v.acc[,3])
      idx.dir <- which(abs(diffs) != 1)
      
      diffs <- diff(v.acc[,1])
      idx.dir <- c(0, sort(c(idx.dir, which(abs(diffs) >= len.break))), nrow(v.acc))
      
      pos.acc = pos.long[pos.long[,2] == acc,1]
      
      idx.dir = sort(c(idx.dir, which(v.acc[,1] %in% pos.acc)))
      
      for(i in 1:(length(idx.dir)-1)){
        pan.b = v.acc[idx.dir[i]+1,2]
        pan.e = v.acc[idx.dir[i+1],2]
        acc.b = v.acc[idx.dir[i]+1,1]
        acc.e = v.acc[idx.dir[i+1],1]
        
        if((acc.e - acc.b + 1) < len.break.gap) next
        idx.synteny = rbind(idx.synteny,
                            data.frame(
                              pan.b = pan.b,
                              pan.e = pan.e,
                              own.b = acc.b,
                              own.e = acc.e,
                              acc = acc,
                              chr = i.chr,
                              # pan.len = pan.e - pan.b - 1,
                              # acc.len = acc.e - acc.b - 1,
                              dir = (v.acc[idx.dir[i+1],3] < v.acc[idx.dir[i]+1,3]) * 1
                            ))
      }
      rm(v.acc)
    }
    
    if(echo) cat('\n')
    
    vars_to_keep <- c("idx.synteny", "i.chr", 'len.break', 'len.break.gap', 'path.aln', 'n.acc')
    all_vars <- ls()
    vars_to_remove <- setdiff(all_vars, vars_to_keep)
    rm(list = vars_to_remove)
    gc()
  }
  
  colnames(idx.synteny)[c(3,4)] = c('own.b', 'own.e')
  return(idx.synteny)
}



segmentSplit <- function(df.tmp, n){
  
  x1 <- df.tmp$x[1]
  y1 <- df.tmp$y.val[1]
  x2 <- df.tmp$x[4]
  y2 <- df.tmp$y.val[4]
  
  # Steps
  step_x <- (x2 - x1) / n
  step_y <- (y2 - y1) / n
  
  # Create points
  x_points <- c(x1, seq(x1 + step_x, x2 - step_x, by = step_x), x2)
  y_points <- c(y1, seq(y1 + step_y, y2 - step_y, by = step_y), y2)
  
  # Save
  points1 <- data.frame(x = x_points, y = y_points)
  
  #  - - - - - - - - - 
  # Second pair
  x1 <- df.tmp$x[2]
  y1 <- df.tmp$y.val[2]
  x2 <- df.tmp$x[3]
  y2 <- df.tmp$y.val[3]
  
  # Steps
  step_x <- (x2 - x1) / n
  step_y <- (y2 - y1) / n
  
  # Create points
  pokaz(x1, step_x, x2)
  pokaz(y1, step_y, y2)
  x_points <- c(x1, seq(x1 + step_x, x2 - step_x/2, by = step_x), x2)
  y_points <- c(y1, seq(y1 + step_y, y2 - step_y/2, by = step_y), y2)
  
  pokaz(x_points)
  pokaz(y_points)
  
  # Save
  points2 <- data.frame(x = x_points, y = y_points)  
  return(list(p1 = points1, p2 = points2))
}


prepareBlocks <- function(idx.break, file.cen.pos=NULL, file.acc.len=NULL, gap.len = 100000){
  
  cen.pos.chr = NULL
  if(!is.null(file.cen.pos) & !is.null(file.acc.len)){
    cen.pos = read.table(file.cen.pos,
                         stringsAsFactors = F, header = 1)
    cen.pos.chr = cen.pos[cen.pos$Chromosome.x == paste0('Chr',i.chr),]
    cen.pos.chr$Accession.x = as.character(cen.pos.chr$Accession.x)
    cen.pos.chr$End2 = 0
    
    acc.len = read.table(file.acc.len, stringsAsFactors = F)
    
    max.len = max(idx.break$own.e[idx.break$chr == i.chr]) + gap.len
    
    for(acc in accessions){
      chr.pos.acc = cen.pos.chr$End[cen.pos.chr$Accession.x == acc]
      
      acc.len.chr = acc.len$V3[(acc.len$V1 == acc) & (acc.len$V2 == i.chr)]
      
      
      idx.tmp = (idx.break$own.b >= chr.pos.acc) & (idx.break$acc == acc) & (idx.break$own.e >= chr.pos.acc)
      idx.break$own.b[idx.tmp] =
        idx.break$own.b[idx.tmp] + max.len - acc.len.chr
      # idx.tmp = (idx.break$own.e >= chr.pos.acc) & (idx.break$acc == acc)
      idx.break$own.e[idx.tmp] =
        idx.break$own.e[idx.tmp] + max.len - acc.len.chr
      
      idx.cen.acc = idx.break[(idx.break$chr == i.chr)&(idx.break$acc == acc),]
      idx.cen.acc = idx.cen.acc[order(idx.cen.acc$own.b),]
      idx.tmp = max(which(idx.cen.acc$own.b <= chr.pos.acc))
      chr.pos.acc2 = idx.cen.acc$own.e[idx.tmp]
      cen.pos.chr$End[cen.pos.chr$Accession.x == acc] = chr.pos.acc2
      cen.pos.chr$End2[cen.pos.chr$Accession.x == acc] = idx.cen.acc$own.b[idx.tmp+1]
      
    }
  }
  
  # Fix direction
  idx.rev = which(idx.break$dir == 1)
  tmp = idx.break$pan.e[idx.rev]
  idx.break$pan.e[idx.rev] = idx.break$pan.b[idx.rev]
  idx.break$pan.b[idx.rev] = tmp
  
  tmp = idx.break$own.e[idx.rev]
  idx.break$own.e[idx.rev] = idx.break$own.b[idx.rev]
  idx.break$own.b[idx.rev] = tmp
  
  
  return(list(idx.break=idx.break, cen.pos.chr=cen.pos.chr))
  
  
}

# Define blocks shared between neighbouring accessions
getBlocksBwNeiAccs <- function(idx.break, i.chr, accessions, i.order){
  
  df.blocks <- c()
  for(k in 2:length(i.order)){
    
    acc1 = accessions[i.order[k-1]]
    acc2 = accessions[i.order[k]]
    idx.break.k = idx.break[(idx.break$chr == i.chr) & 
                              ((idx.break$acc == acc1) |
                                 (idx.break$acc == acc2)),]
    
    idx.break.k = idx.break.k[order(idx.break.k$pan.b),]
    
    # Collect all pangen breakes between these accessions
    
    for(i.rep in 1:2){
      idx1 = (idx.break.k$acc == acc1)
      idx2 = (idx.break.k$acc == acc2)
      
      br1 = idx.break.k[idx1,]
      
      # add beginning breaks
      pan.br.add1 = c(setdiff(idx.break.k$pan.b[idx2], idx.break.k$pan.b[idx1]),
                      setdiff(idx.break.k$pan.e[idx2], idx.break.k$pan.e[idx1]) + 1)
      if(length(pan.br.add1) > 0){
        for(pos.br in pan.br.add1){
          i = which((br1$pan.b < pos.br) & ((br1$pan.e+1) > pos.br))
          if(length(i) == 0){
            # message(paste('NULL', acc1, acc2, pos.br))
            next
          }
          br1 = rbind(br1, br1[i,])
          j = nrow(br1)
          br1$pan.e[i] = pos.br-1
          br1$pan.b[j] = pos.br
          
          # Proportion of length 
          p = (pos.br - br1$pan.b[i]) / (br1$pan.e[j] - br1$pan.b[i])
          pos.br.own = round(br1$own.b[i] + p * (br1$own.e[i] - br1$own.b[i]))
          br1$own.e[i] = pos.br.own - 1
          br1$own.b[j] = pos.br.own
          
        }
      }
      
      br1 = br1[order(br1$pan.b),]
      idx.break.k = idx.break.k[!idx1,]
      idx.break.k = rbind(idx.break.k, br1)
      
      tmp = acc1
      acc1 = acc2
      acc2 = tmp
      
    }
    
    idx.break.k$sum = idx.break.k$pan.b + idx.break.k$pan.e
    sum.dup = idx.break.k$sum[duplicated(idx.break.k$sum)]
    idx.break.k = idx.break.k[idx.break.k$sum %in% sum.dup, ]
    
    if(length(sum.dup) == 0) next
    
    for(s in sum.dup){
      # print(s)
      b1 = idx.break.k[(idx.break.k$sum == s) & (idx.break.k$acc == acc1),]
      b2 = idx.break.k[(idx.break.k$sum == s) & (idx.break.k$acc == acc2),]
      df.block = data.frame(
        pan.b = b1$pan.b,
        pan.e = b1$pan.e,
        own1.b = b1$own.b,
        own1.e = b1$own.e,
        own2.b = b2$own.b,
        own2.e = b2$own.e,
        y1 = k-1,
        y2 = k,
        acc1 = b1$acc,
        acc2 = b2$acc,
        dir1 = b1$dir == 0,
        dir2 = b2$dir == 0
      )
      df.blocks = rbind(df.blocks, df.block)
    }
    
  }
  
  return(df.blocks)

}


splitBlocksByGrid <- function(df.blocks, wnd.size = 1000000){
  
  # message(paste('window size', wnd.size))
  n.bl = nrow(df.blocks)
  for(irow in 1:n.bl){
    d = df.blocks$pan.e[irow] - df.blocks$pan.b[irow] 
    d.own1 = df.blocks$own1.e[irow] - df.blocks$own1.b[irow]
    d.own2 = df.blocks$own2.e[irow] - df.blocks$own2.b[irow]
    if(abs(d) > wnd.size){
      n.repl = ceiling(abs(d) / wnd.size) 
      df.add <- do.call(rbind, replicate(n.repl, df.blocks[irow,], simplify = FALSE))
      df.add$pan.b = round(df.add$pan.b + (0:(n.repl-1))*d/(n.repl))
      df.add$pan.e[1:(n.repl-1)] = df.add$pan.b[2:n.repl] - 1
      
      
      df.add$own1.b = round(df.add$own1.b + (0:(n.repl-1))*d.own1/(n.repl))
      df.add$own1.e[1:(n.repl-1)] = df.add$own1.b[2:n.repl] - 1
      
      df.add$own2.b = round(df.add$own2.b + (0:(n.repl-1))*d.own2/(n.repl))
      df.add$own2.e[1:(n.repl-1)] = df.add$own2.b[2:n.repl] - 1
      df.blocks[irow, 1] = NA
      df.blocks = rbind(df.blocks, df.add)
    }
  }
  df.blocks = df.blocks[!is.na(df.blocks[,1]),]
  
  
  df.plot = c()
  for(irow in 1:nrow(df.blocks)){
    bl = data.frame(x = c(df.blocks$own1.b[irow], df.blocks$own1.e[irow], 
                          df.blocks$own2.e[irow], df.blocks$own2.b[irow]),
                    y = c(df.blocks$acc1[irow], df.blocks$acc1[irow], 
                          df.blocks$acc2[irow], df.blocks$acc2[irow]),
                    t = irow, chr = i.chr,
                    pan.b = df.blocks$pan.b[irow],
                    # pan.e = df.blocks$pan.e[irow],
                    dir1 = df.blocks$dir1[irow], dir2 = df.blocks$dir2[irow]
    )
    df.plot = rbind(df.plot, bl)
  }
  
  return(df.plot)
  
}


getColorPallete <- function(df.plot){
  custom_colors.cold = c('#27374D',rep(c('#526D82', '#9DB2BF', '#DDE6ED', '#9DB2BF', '#526D82', '#27374D'), 6))
  
  # custom_colors.warm = c('#CE1F6A',
  #                        rep(c('#E76161', '#F99B7D', '#FFE4C0', '#F99B7D', '#E76161', '#CE1F6A'), 5))
  
  custom_colors.warm <- c("#CE1F6A")
  grad.cold <- colorRamp(custom_colors.cold)
  grad.warm <- colorRamp(custom_colors.warm)
  
  
  
  df.plot$c.val = df.plot$pan.b
  df.plot$c.val[!(df.plot$dir1 & df.plot$dir2)] = df.plot$c.val[!(df.plot$dir1 & df.plot$dir2)] * (-1)
  df.plot$col = ''
  
  df.plot$col[df.plot$c.val > 0] = rgb(grad.cold(df.plot$c.val[df.plot$c.val > 0]/ max(df.plot$c.val)),
                                       maxColorValue = 255)
  
  df.plot$col[df.plot$c.val < 0] = rgb(grad.warm(df.plot$c.val[df.plot$c.val < 0]/ min(df.plot$c.val)),
                                       maxColorValue = 255)
  
  gr.col = unique(df.plot[,c('t', 'col')])
  gr.col = gr.col[order(gr.col$t),]
  gr.col = gr.col[,2]
  return(gr.col)
}

# Define vertical shades for inversions
splitInversionBlocks <- function(df.plot, idx.break, gr.col, n.split=5){
  
  
  # Slit blocks
  colors.split <- c('#9DB2BF', "#CE1F6A")
  
  grad.split <- colorRamp(colors.split)
  
  t.correct = unique(df.plot$t[df.plot$dir1 & !df.plot$dir2])
  for(i.rep in 1:2){
    for(i in t.correct){
      df.tmp = df.plot[df.plot$t == i,]
      
      points = segmentSplit(df.tmp, n.split)
      
      df.add111 = c()
      for(j in 1:n.split){
        df.add = df.tmp
        df.add$t = max(df.plot$t) + 1
        df.add$x[1] = points$p1$x[j]
        df.add$x[4] = points$p1$x[j+1]
        df.add$x[2] = points$p2$x[j]
        df.add$x[3] = points$p2$x[j+1]
        
        df.add$y.val[1] = points$p1$y[j]
        df.add$y.val[4] = points$p1$y[j+1]
        df.add$y.val[2] = points$p2$y[j]
        df.add$y.val[3] = points$p2$y[j+1]
        
        # col.new = rgb(grad.split((j-1)/(n.split-1)), maxColorValue = 255)
        if(i.rep == 1){
          col.new = rgb(grad.split((j)/(n.split)), maxColorValue = 255)
        } else {
          col.new = rgb(grad.split(1-(j-1)/(n.split)), maxColorValue = 255)
        }
        
        gr.col = c(gr.col, col.new)
        
        # print(df.add)
        df.add111 = rbind(df.add111, df.add)
        df.plot = rbind(df.plot, df.add)
      }
    }
    gr.col[t.correct] = '#FFFFFF'
    
    t.correct = unique(df.plot$t[!df.plot$dir1 & df.plot$dir2])
    
  }
  return(list(df.plot, gr.col))
}

#' Visualization of synteny blocks along genomes
#'
#' This function visualizes synteny blocks along the specified chromosome with optional centromere positions. 
#' It generates a ggplot object showing synteny blocks between accessions, coloring blocks based on their 
#' relative orientation and positions.
#'
#' @param idx.break Data frame with information on the breakpoints between syntenic blocks.
#' @param i.chr The chromosome number or identifier for which synteny blocks will be plotted.
#' @param accessions A character vector with the names or identifiers of the accessions to be included in the plot.
#' @param i.order Optional, an integer vector specifying the order in which accessions will be plotted.
#' @param wnd.size Numeric, window size for splitting synteny blocks for visualization, default is 1000000.
#' @return A ggplot object showing synteny blocks, their positions, and optionally centromeric regions.
#' 
#' @examples
#' # Example usage:
#' ggSynteny(idx.break, i.chr = 1, accessions = c("Acc1", "Acc2"), i.order = c(1, 2))
#'
#' @export
panplot <- function(idx.break, i.chr, accessions=NULL, i.order=NULL, file.cen.pos=NULL, file.acc.len=NULL, 
                      gap.len = 100000, wnd.size = 1000000){
  
  gap.len = min(gap.len, wnd.size / 10)
  if(is.null(accessions)){
    accessions = unique(idx.break$acc)
  }
  if(is.null(i.order)){
    i.order = 1:length(accessions)
  }
  # Change position with respect to centromeric position
  res = prepareBlocks(idx.break, file.cen.pos, file.acc.len)
  idx.break = res[[1]]
  cen.pos.chr = res[[2]]
  
  # Blocks into sub-blocks between neighbouring accessions
  df.blocks = getBlocksBwNeiAccs(idx.break, i.chr, accessions, i.order)
  
  # Blocks by grid
  # wnd.size = 100000
  
  df.plot = splitBlocksByGrid(df.blocks, wnd.size = wnd.size)
  
  # Define color pallete
  gr.col = getColorPallete(df.plot)
  
  # Split vertically inversion blocks
  
  acc.ord = accessions[i.order]
  df.plot$y = factor(df.plot$y, levels = acc.ord)
  idx.break$acc = factor(idx.break$acc, levels = acc.ord)
  df.plot$y.val = as.numeric(factor(df.plot$y, levels = acc.ord))
  idx.break$acc.val = as.numeric(factor(idx.break$acc, levels = acc.ord))
  
  res = splitInversionBlocks(df.plot, idx.break, gr.col)
  
  df.plot = res[[1]]
  gr.col = res[[2]]
  
  # pokaz(head(df.plot))
  
  # p +  theme(legend.position = "none") 
  
  p = ggplot() +
    geom_polygon(data=df.plot[df.plot$chr == i.chr,], 
                 mapping=aes(x=x, y=y.val, group=t, 
                             # fill = as.factor(dir1+dir2) )
                             fill = as.factor(t)), 
                 # color = 'grey',
                 alpha=0.5) + 
    geom_segment(data = idx.break[idx.break$chr == i.chr,], 
                 aes(x=own.b, y=acc.val, xend=own.e, yend=acc.val, 
                     color=as.factor(( dir == 0) ))) + #+ 2 * trans
    # color=pan.b)) + 
    scale_color_manual(values = c('#CE1F6A', 'black')) + 
    scale_fill_manual(values = gr.col) + 
    # scale_fill_gradientn(colors = custom_colors.cold) + 
    # scale_fill_gradientn(colors = custom_colors.cold) + 
    # scale_fill_viridis_c() + 
    theme_minimal() + ylab('') + xlab('own genome position') + 
    theme(legend.position = "none") + theme(panel.grid.minor.x = element_blank(), 
                                            panel.grid.minor.y = element_blank(), 
                                            panel.grid.major.x = element_blank()) + 
    scale_y_continuous(labels =  accessions[i.order], breaks = 1:length(i.order), expand = c(0, 0)) + 
    scale_x_continuous(expand = c(0, 0))
  
  
  if(!is.null(cen.pos.chr)){
    df.cen = cen.pos.chr
    df.cen$Accession.x = factor(df.cen$Accession.x, levels = acc.ord)
    df.cen$y.val = as.numeric(df.cen$Accession.x)
    df.cen$End = df.cen$End #+ gap.len * 2 
    df.cen = df.cen[order(df.cen$y.val),]
    p = p + geom_ribbon(data = df.cen, aes(y = y.val, xmin = End, xmax = End2), fill = "#F79327", orientation='y', alpha = 0.3)
  }
  
  return(p)
  
}

#' Visualization of synteny blocks along genomes
#' (the same as panplot)
#' @export
ggSynteny <- function(...) {
  pokazAttention('Please exchange ggSynteny to panplot')
  panplot(...)
}



