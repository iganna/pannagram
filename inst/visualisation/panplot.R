#' This file contains functions for visualising synteny between genomes
#' len.min = 20000  Le length of which gap to consider as to blue the regions
#' @export
getBlocks <- function(v, f.split = T, len.min = 10000){
  
  # v <- h5read(file.comb.in, paste0(gr.accs.e, acc))
  
  v.init = v
  v = v.init
  v.idx = 1:length(v)
  
  v.idx = v.idx[v != 0]
  v = v[v != 0]
  
  v.r = v.idx
  v.r = v.idx * sign(v)
  
  v.b = findRuns(v.r)
  
  v.b[,'beg'] = v[v.b[,'beg']]
  v.b[,'end'] = v[v.b[,'end']]
  

  # vals <- c(abs(v.b$v.beg), abs(v.b$v.end))
  vals <- c(abs(v.b$beg), abs(v.b$end))
  ranks <- match(vals, sort(unique(vals))) 
  v.b.rank <- matrix(
    ranks,
    ncol = 2,
    byrow = FALSE,
    dimnames = list(NULL, c("r.beg", "r.end")))
  
  v.b = cbind(v.b, v.b.rank)

  # v.b$dir = sign(v.b$v.end - v.b$v.beg)
  v.b$dir = sign(v.b$beg)
  v.b$v.beg = abs(v.b$v.beg)
  v.b$v.end = abs(v.b$v.end)
  
  # Merge blocks
  n <- nrow(v.b)
  cond <-
    (v.b$r.beg[2:n] - 1 == v.b$r.end[1:(n-1)]) &
    (v.b$dir[2:n]     == v.b$dir[1:(n-1)])     &
    ((v.b$v.beg[2:n] - 1 - v.b$v.end[1:(n-1)]) < len.min)
  
  starts <- c(1, which(!cond) + 1)
  ends <- c(starts[-1] - 1, n)
  
  v.b.new = v.b[starts,]
  
  v.b.new$r.end = v.b$r.end[ends]
  v.b.new$end = v.b$end[ends]
  v.b.new$v.end = v.b$v.end[ends]
  v.b.new$len = v.b.new$end - v.b.new$beg + 1
  rownames(v.b.new) = NULL
  
  
  df = v.b.new[,c('beg', 'end', 'v.beg', 'v.end')]
  df = abs(df)
  rownames(df) = NULL
  colnames(df) <- c('own.b', 'own.e', 'pan.b', 'pan.e')
  df$dir = (df$own.b > df$own.e) * 1
  return(df)
  
  # 
  # 
  # df = v.b[,c('v.beg', 'v.end', 'i.beg', 'i.end')]
  # 
  # # Add breaks, if they are with a long gap
  # if(f.split){
  #   len.min.split = 50000
  #   d = diff(v.idx)
  #   i.d = which(d >= len.min.split)
  #   for(i.d in which(d >= len.min.split)){
  #     v.pos = v[i.d]
  #     i.pos = v.idx[i.d]
  #     
  #     i.v = which((df$v.beg < v.pos) & (df$v.end > v.pos))
  #     i.i = which((df$i.beg < i.pos) & (df$i.end > i.pos))
  #     
  #     if(length(i.v) > 0){
  #       if(length(i.i) == 0){
  #         # save(list = ls(), file = "tmp_workspace_blocks.RData")
  #         # stop('Something is wrong with i.i')
  #         pokazAttention('Something is wrong with i.i')
  #         next
  #         # I have found that it's the case, 
  #         # when i.v points to something between two blocks, 
  #         # so there is no block to split
  #       } 
  #       if(i.v != i.i) stop('Idxs do not match')
  #       
  #       df.tmp = df[i.v,]
  #       
  #       df.tmp$v.beg = abs(v[i.d + 1])
  #       df.tmp$i.beg = v.idx[i.d + 1]
  #       
  #       df[i.v,]$v.end = v.pos
  #       df[i.v,]$i.end = i.pos
  #     
  #       df = rbind(df, df.tmp)
  #     }
  #   }    
  # }
  # 
  # df = df[order(df$i.beg),]
  # rownames(df) = NULL
  # 
  # colnames(df) <- c('own.b', 'own.e', 'pan.b', 'pan.e')
  # df$dir = (df$own.b > df$own.e) * 1
  # return(df)
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
  x_points <- c(x1, seq(x1 + step_x, x2 - step_x/2, by = step_x), x2)
  y_points <- c(y1, seq(y1 + step_y, y2 - step_y/2, by = step_y), y2)
  
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
  x_points <- c(x1, seq(x1 + step_x, x2 - step_x/2, by = step_x), x2)
  y_points <- c(y1, seq(y1 + step_y, y2 - step_y/2, by = step_y), y2)
  
  # Save
  points2 <- data.frame(x = x_points, y = y_points)  
  return(list(p1 = points1, p2 = points2))
}


prepareBlocks <- function(idx.break, file.cen.pos=NULL, file.acc.len=NULL){
  
  cen.pos.chr = NULL
  if(!is.null(file.cen.pos) & !is.null(file.acc.len)){
    cen.pos = read.table(file.cen.pos,
                         stringsAsFactors = F, header = 1)
    cen.pos.chr = cen.pos[cen.pos$Chromosome.x == paste0('Chr',i.chr),]
    cen.pos.chr$Accession.x = as.character(cen.pos.chr$Accession.x)
    cen.pos.chr$End2 = 0
    
    acc.len = read.table(file.acc.len, stringsAsFactors = F)
    
    max.len = max(idx.break$own.e[idx.break$chr == i.chr])
    
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
getBlocksBwNeiAccs <- function(idx.break, accessions, i.order){
  
  df.blocks <- c()
  for(k in 2:length(i.order)){
    pokaz(k)
    acc1 = accessions[i.order[k-1]]
    acc2 = accessions[i.order[k]]
    idx.break.k = idx.break[(idx.break$acc == acc1) |
                            (idx.break$acc == acc2),]
    
    idx.break.k = idx.break.k[order(idx.break.k$pan.b),]
    
    # Collect all pangen breaks between these accessions
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
    
    # Extract all rows for acc1 and acc2 that have duplicated sum values
    b1 <- idx.break.k[idx.break.k$acc == acc1 & idx.break.k$sum %in% sum.dup, ]
    b2 <- idx.break.k[idx.break.k$acc == acc2 & idx.break.k$sum %in% sum.dup, ]
    
    # Align rows of b2 so that they match b1 by 'sum'
    m <- match(b1$sum, b2$sum)
    b2 <- b2[m, ]
    
    # Build all block rows at once
    df.block <- data.frame(
      pan.b  = b1$pan.b,
      pan.e  = b1$pan.e,
      own1.b = b1$own.b,
      own1.e = b1$own.e,
      own2.b = b2$own.b,
      own2.e = b2$own.e,
      y1     = k - 1,
      y2     = k,
      acc1   = b1$acc,
      acc2   = b2$acc,
      dir1   = b1$dir == 0,
      dir2   = b2$dir == 0
    )
    
    df.blocks <- rbind(df.blocks, df.block)
    
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
  
  df.plot <- data.frame(
    x = c(df.blocks$own1.b, df.blocks$own1.e, df.blocks$own2.e, df.blocks$own2.b),
    y = c(df.blocks$acc1,    df.blocks$acc1,    df.blocks$acc2,    df.blocks$acc2),
    t = rep(seq_len(nrow(df.blocks)), 4),
    pan.b = rep(df.blocks$pan.b, 4),
    dir1  = rep(df.blocks$dir1,  4),
    dir2  = rep(df.blocks$dir2,  4)
  )
  df.plot = df.plot[order(df.plot$t),]
  rownames(df.plot) = NULL
  
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
splitInversionBlocks <- function(df.plot, idx.break, gr.col, n.split = 5) {
  colors.split <- c("#9DB2BF", "#CE1F6A")
  grad.split   <- colorRamp(colors.split)
  
  # будем создавать новые ряды в списке, а не через rbind в цикле
  new_rows_list <- list()
  new_cols_vec  <- character()
  
  max_t <- max(df.plot$t)
  t_counter <- 0
  
  t.correct <- unique(df.plot$t[df.plot$dir1 & !df.plot$dir2])
  
  for (i.rep in 1:2) {
    # для текущего i.rep собираем новые ряды во временный список
    rep_rows_list <- vector("list", length(t.correct) * n.split)
    rep_cols_vec  <- character(length(t.correct) * n.split)
    k <- 1
    
    for (i in t.correct) {
      df.tmp <- df.plot[df.plot$t == i, ]
      points <- segmentSplit(df.tmp, n.split)
      
      for (j in seq_len(n.split)) {
        df.add <- df.tmp
        
        # задаём новый t без постоянного max(df.plot$t)
        t_counter <- t_counter + 1
        df.add$t <- max_t + t_counter
        
        df.add$x[1]     <- points$p1$x[j]
        df.add$x[4]     <- points$p1$x[j + 1]
        df.add$x[2]     <- points$p2$x[j]
        df.add$x[3]     <- points$p2$x[j + 1]
        
        df.add$y.val[1] <- points$p1$y[j]
        df.add$y.val[4] <- points$p1$y[j + 1]
        df.add$y.val[2] <- points$p2$y[j]
        df.add$y.val[3] <- points$p2$y[j + 1]
        
        if (i.rep == 1) {
          col.new <- rgb(grad.split(j / n.split), maxColorValue = 255)
        } else {
          col.new <- rgb(grad.split(1 - (j - 1) / n.split), maxColorValue = 255)
        }
        
        rep_rows_list[[k]] <- df.add
        rep_cols_vec[k]    <- col.new
        k <- k + 1
      }
    }
    
    # добавляем то, что насобирали за данный i.rep
    new_rows_list <- c(new_rows_list, rep_rows_list)
    new_cols_vec  <- c(new_cols_vec,  rep_cols_vec)
    
    # перекрашиваем "исходные" блоки в белый
    gr.col[t.correct] <- "#FFFFFF"
    
    # обновляем t.correct для следующего прохода
    t.correct <- unique(df.plot$t[!df.plot$dir1 & df.plot$dir2])
  }
  
  # один раз склеиваем всё, что насобирали
  if (length(new_rows_list) > 0) {
    df.new   <- do.call(rbind, new_rows_list)
    df.plot  <- rbind(df.plot, df.new)
    gr.col   <- c(gr.col, new_cols_vec)
  }
  
  return(list(df.plot, gr.col))
}

# splitInversionBlocks <- function(df.plot, idx.break, gr.col, n.split=5){
#   
#   
#   # Slit blocks
#   colors.split <- c('#9DB2BF', "#CE1F6A")
#   
#   grad.split <- colorRamp(colors.split)
#   
#   t.correct = unique(df.plot$t[df.plot$dir1 & !df.plot$dir2])
#   for(i.rep in 1:2){
#     for(i in t.correct){
#       df.tmp = df.plot[df.plot$t == i,]
#       
#       points = segmentSplit(df.tmp, n.split)
#       
#       for(j in 1:n.split){
#         df.add = df.tmp
#         df.add$t = max(df.plot$t) + 1
#         df.add$x[1] = points$p1$x[j]
#         df.add$x[4] = points$p1$x[j+1]
#         df.add$x[2] = points$p2$x[j]
#         df.add$x[3] = points$p2$x[j+1]
#         
#         df.add$y.val[1] = points$p1$y[j]
#         df.add$y.val[4] = points$p1$y[j+1]
#         df.add$y.val[2] = points$p2$y[j]
#         df.add$y.val[3] = points$p2$y[j+1]
#         
#         # col.new = rgb(grad.split((j-1)/(n.split-1)), maxColorValue = 255)
#         if(i.rep == 1){
#           col.new = rgb(grad.split((j)/(n.split)), maxColorValue = 255)
#         } else {
#           col.new = rgb(grad.split(1-(j-1)/(n.split)), maxColorValue = 255)
#         }
#         
#         gr.col = c(gr.col, col.new)
#         
#         # print(df.add)
#         df.plot = rbind(df.plot, df.add)
#       }
#     }
#     gr.col[t.correct] = '#FFFFFF'
#     
#     t.correct = unique(df.plot$t[!df.plot$dir1 & df.plot$dir2])
#     
#   }
#   return(list(df.plot, gr.col))
# }

#' Visualization of synteny blocks along genomes
#'
#' This function visualizes synteny blocks along the specified chromosome with optional centromere positions. 
#' It generates a ggplot object showing synteny blocks between accessions, coloring blocks based on their 
#' relative orientation and positions.
#'
#' @param idx.break Data frame with information on the breakpoints between syntenic blocks.
#' @param accessions A character vector with the names or identifiers of the accessions to be included in the plot.
#' @param i.order Optional, an integer vector specifying the order in which accessions will be plotted.
#' @param wnd.size Numeric, window size for splitting synteny blocks for visualization, default is 1000000.
#' @return A ggplot object showing synteny blocks, their positions, and optionally centromeric regions.
#' 
panplotInner <- function(idx.break, accessions=NULL, i.order=NULL, file.cen.pos=NULL, file.acc.len=NULL, 
                      wnd.size = 1000000){
  
  wnd.size.min = 10000
  wnd.size = max(wnd.size, wnd.size.min)
  
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
  
  # Filter accessions
  idx.break = idx.break[idx.break$acc %in% accessions,]
  
  # Blocks into sub-blocks between neighbouring accessions
  pokaz('getBlocksBwNeiAccs')
  df.blocks = getBlocksBwNeiAccs(idx.break, accessions, i.order)
  
  # Blocks by grid
  pokaz('splitBlocksByGrid')
  df.plot = splitBlocksByGrid(df.blocks, wnd.size = wnd.size)
  
  # Define color pallete
  pokaz('getColorPallete')
  gr.col = getColorPallete(df.plot)
  
  # Split vertically inversion blocks
  
  acc.ord = accessions[i.order]
  df.plot$y = factor(df.plot$y, levels = acc.ord)
  idx.break$acc = factor(idx.break$acc, levels = acc.ord)
  df.plot$y.val = as.numeric(factor(df.plot$y, levels = acc.ord))
  idx.break$acc.val = as.numeric(factor(idx.break$acc, levels = acc.ord))
  
  pokaz('splitInversionBlocks')
  res = splitInversionBlocks(df.plot, idx.break, gr.col)
  
  df.plot = res[[1]]
  gr.col = res[[2]]
  
  # Color of synteny blocks
  idx.break$dir.factor = factor(idx.break$dir * 1, levels = c(0, 1))
  
  pokaz('Plot')
  p = ggplot() +
    geom_polygon(data=df.plot, 
                 mapping=aes(x=x, y=y.val, group=t, 
                             # fill = as.factor(dir1+dir2) )
                             fill = as.factor(t)), 
                 # color = 'grey',
                 alpha=0.5) + 
    geom_segment(data = idx.break, 
                 aes(x=own.b, y=acc.val, xend=own.e, yend=acc.val, 
                     color=dir.factor )) + #+ 2 * trans
    # color=pan.b)) + 
    scale_color_manual(values = c('#CE1F6A','black'), breaks = c(1, 0)) + 
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
    df.cen$End = df.cen$End
    df.cen = df.cen[order(df.cen$y.val),]
    p = p + geom_ribbon(data = df.cen, aes(y = y.val, xmin = End, xmax = End2), fill = "#F79327", orientation='y', alpha = 0.3)
  }
  
  return(p)
  
}


panplot <- function(path.project, i.chr, accessions = NULL, aln.type='pan', ref.acc='',
                       wnd.size = 100000){
  path.inter.msa = paste0(path.project, '/.intermediate/msa/')
  
  ref.suff <- if (ref.acc == '') '' else paste0('_', ref.acc)
  if(ref.suff != '') aln.type='ref'
  
  file.blocks = paste0(path.inter.msa, aln.type, '_syn_blocks', ref.suff,'.rds')
  
  idx.break = readRDS(file.blocks)
  idx.break = idx.break[idx.break$comb == paste0(i.chr, '_', i.chr),]
  
  p = panplotInner(idx.break, wnd.size = wnd.size, accessions=accessions)
  
  return(p)
  
}

