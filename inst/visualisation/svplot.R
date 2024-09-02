svplotMain <- function(b.graph.sub, sv.graph,
                       te.palette, 
                       fam.palette, 
                       indel.pallete,
                       i.target = NULL){
  
  g.sub <- network(b.graph.sub, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
  g.sub.names = network.vertex.names(g.sub)
  
  g.sub %v% "label.size" = as.character(round((sv.graph$nodes.traits[g.sub.names,]$len / 1000),1))
  g.sub %v% "color.type" = as.character(sv.graph$nodes.traits[g.sub.names,]$type)
  
  set.seed(20)
  p.te <- ggnet2(g.sub, 
                 label = "label.size", 
                 label.size = 2,
                 label.color = 'white',
                 edge.color = "black", 
                 # node.size = sv.graph$nodes.traits[g.sub.names,]$cnt,
                 node.size = 6,
                 color =  "color.type",
                 # palette = te.palette,
                 mode = 'kamadakawai',
                 arrow.gap = 0.07, arrow.size = 3
  ) +
    scale_colour_manual(values = te.palette) +
    guides(size = F) +
    theme_void() +
    scale_x_continuous(expand = c(0.1, 0.1)) +
    scale_y_continuous(expand = c(0.1, 0.1)) +
    theme(legend.key.size = unit(0.5, "cm"))  +
    guides(color = guide_legend(ncol = 2 )) + 
    theme(legend.position = "bottom") + 
    labs(caption = "Value: Length, kbp") +
    theme(plot.caption = element_text(hjust = 0.5))
  
  
  g.sub %v% "label.count" = as.character(sv.graph$nodes.traits[g.sub.names,]$cnt)
  g.sub %v% "color.fam" = as.character(sv.graph$nodes.traits[g.sub.names,]$fam)
  g.sub %v% "node.size.count" = sv.graph$nodes.traits[g.sub.names,]$cnt
  
  set.seed(20)
  p.fam <- ggnet2(g.sub, 
                  label = "label.count", 
                  edge.color = "black", 
                  node.size = "node.size.count", 
                  label.size = 2,
                  label.color = 'white',
                  color =  "color.fam",
                  # palette = fam.palette,
                  mode = 'kamadakawai',
                  arrow.gap = 0.06, arrow.size = 3
  ) +
    scale_colour_manual(values = fam.palette) +
    guides(size = F) +
    theme_void() +
    scale_x_continuous(expand = c(0.1, 0.1)) +
    scale_y_continuous(expand = c(0.1, 0.1)) +
    # theme_minimal() + ylab(NULL) + xlab('Number: node size') +
    theme(legend.key.size = unit(0.5, "cm"))  +
    guides(color = guide_legend(ncol = 2 )) + 
    theme(legend.position = "bottom") + 
    labs(caption = "Value: Number of indels") +
    theme(plot.caption = element_text(hjust = 0.5))
  
  
  g.sub %v% "label.freq" = as.character(sv.graph$nodes.traits[g.sub.names,]$freq)
  g.sub %v% "color.indel" = as.character(indel.new.names[sv.graph$nodes.traits[g.sub.names,]$indel])
  
  set.seed(20)
  p.freq <- ggnet2(g.sub, 
                   label = "label.freq",
                   label.color = 'white',
                   label.size = 2,
                   edge.color = "black", 
                   # node.size = sv.graph$nodes.traits[g.sub.names,]$cnt, 
                   node.size = 6,
                   color =  "color.indel",
                   # palette = indel.pallete,
                   mode = 'kamadakawai',
                   arrow.gap = 0.07, arrow.size = 3
  ) +
    scale_colour_manual(values = indel.pallete) +
    guides(size = F) +
    theme_void() +
    scale_x_continuous(expand = c(0.1, 0.1)) +
    scale_y_continuous(expand = c(0.1, 0.1)) +
    # theme_minimal() + ylab(NULL) + xlab('Number: freq of presence') +
    theme(legend.key.size = unit(0.5, "cm"))  +
    guides(color = guide_legend(ncol = 2 )) + 
    theme(legend.position = "bottom")+ 
    labs(caption = "Value: Mean presence freq.") +
    theme(plot.caption = element_text(hjust = 0.5))
  
  if(!is.null(i.target)){
    p.te = p.te + annotate("text", x = -Inf, y = Inf, 
                             label = i.target, hjust = 0, vjust = 1)
  }
  
  pp = invisible(ggarrange(p.te  + labs(color = NULL),
                           p.fam + labs(color = NULL),
                           p.freq+ labs(color = NULL),
                           nrow=1))
  return(pp)
}



svplotTwoNodes <- function(s.nodes.target, sv.graph, g.comp=NULL,
                           n1 = NULL, n2 = NULL,
                           sv.on.te, 
                           sv.seqs, te.seq,
                           fam.palette){
  
  # s.nodes.target = names(g.comp$membership[g.comp$membership == i.target])
  df.nodes.target = sv.graph$nodes.traits[sv.graph$nodes.traits$node %in% s.nodes.target,]
  
  df.nodes.target = df.nodes.target[order(-df.nodes.target$len),]
  df.nodes.target = df.nodes.target[order(-df.nodes.target$cnt),]
  
  if(is.null(n1)){
    n1 = df.nodes.target$node[1]
  }
  
  if(is.null(n2)){
    n2 = df.nodes.target$node[2]
  }
  
  
  # n1 = 'N269'
  # n2 = 'N203'
  
  
  for(i.tmp in 1:10000){
    sv1 = sv.graph$nodes$name[sv.graph$nodes$node == n1][i.tmp]
    if(sum(sv.on.te$V1 == sv1 ) != 0) break
  }
  
  for(i.tmp in 1:10000){
    sv2 = sv.graph$nodes$name[sv.graph$nodes$node == n2][i.tmp]
    if(sum(sv.on.te$V1 == sv2 ) != 0) break
  }
  
  len1 = as.numeric(strsplit(sv1, '\\|')[[1]][2])
  len2 = as.numeric(strsplit(sv2, '\\|')[[1]][2])
  
  # Name plot
  if(is.null(g.comp)){
    p.name = NULL
  } else {
    b.graph.sub = sv.graph$edges[(sv.graph$edges[,1] %in% s.nodes.target) | 
                                   (sv.graph$edges[,2] %in% s.nodes.target),,drop=F]
    
    g.sub <- network(b.graph.sub, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
    g.sub.names = network.vertex.names(g.sub)
    g.sub %v% "label.name" = as.character(sv.graph$nodes.traits[g.sub.names,]$node)
    g.sub %v% "color.visual" = c('hide', 'show')[(g.sub.names %in% c(n1, n2)) * 1 + 1]
    set.seed(20)
    p.name <- ggnet2(g.sub, 
                     label = "label.name",
                     label.color = 'black',
                     label.size = 3,
                     # edge.color = "black", 
                     # node.size = sv.graph$nodes.traits[g.sub.names,]$cnt, 
                     node.size = 6,
                     color =  "color.visual",
                     # palette = indel.pallete,
                     mode = 'kamadakawai',
                     arrow.gap = 0.07, arrow.size = 3) +
      scale_colour_manual(values = c('hide' = 'white', 'show'='grey')) +
      guides(size = F) +
      theme_void() +
      scale_x_continuous(expand = c(0.1, 0.1)) +
      scale_y_continuous(expand = c(0.1, 0.1)) +
      # theme_minimal() + ylab(NULL) + xlab('Number: freq of presence') +
      theme(legend.key.size = unit(0.5, "cm"))  +
      guides(color = guide_legend(ncol = 2 )) + 
      theme(legend.position = "none")
    
  }
  
  # Dot plot
  
  seq1 = seq2nt(sv.seqs[sv1])
  seq2 = seq2nt(sv.seqs[sv2])
  
  p1 = dotplot(seq1, seq2, 15, 13) + ggtitle(paste(n1, 'vs', n2)) + 
    theme(plot.title = element_text(size = 9))
  p2 = dotplot(seq1, seq1, 15, 13) + ggtitle(paste(n1, 'vs', n1)) + 
    theme(plot.title = element_text(size = 9))
  p3 = dotplot(seq2, seq2, 15, 13) + ggtitle(paste(n2, 'vs', n2)) + 
    theme(plot.title = element_text(size = 9))
  
  
  
  df = sv.on.te[sv.on.te$V1 == sv1,c('V2', 'V3', 'V8'),]
  colnames(df) = c('beg', 'end', 'V8')
  df$len = df$end - df$beg
  te1 = df$V8[df$len == max(df$len)][1]
  df$fam = sapply(df$V8, function(s) strsplit(s, '\\|')[[1]][9])
  po1 = orfplot(df, s.color = 'fam', show.legend = T, arrow.size = 0) + 
    scale_colour_manual(values = fam.palette) + 
    xlim(c(0, len1)) + labs(color = NULL) +
    geom_vline(xintercept=c(0, len1), color = 'black', linetype="dashed") +
    ggtitle(paste(n1, ':', sv1)) + theme(legend.position = "bottom") +
    theme(plot.title = element_text(size = 9))
  
  
  df = sv.on.te[sv.on.te$V1 == sv2,c('V2', 'V3', 'V8'),]
  colnames(df) = c('beg', 'end', 'V8')
  df$len = df$end - df$beg
  te2 = df$V8[df$len == max(df$len)][1]
  df$fam = sapply(df$V8, function(s) strsplit(s, '\\|')[[1]][9])
  po2 = orfplot(df, s.color = 'fam', show.legend = T, arrow.size = 0) + 
    scale_colour_manual(values = fam.palette) + 
    xlim(c(0, len2)) + labs(color = NULL) +
    geom_vline(xintercept=c(0, len2), color = 'black', linetype="dashed") +
    ggtitle(paste(n2, ':', sv2)) + theme(legend.position = "bottom") +
    theme(plot.title = element_text(size = 9))
  
  pp.stat = invisible(ggarrange(p.name, po1, po2, p1, p2, p3, ncol = 3, nrow=2))
  
  
  
  
  # Maximum TE
  
  te.info1 = strsplit(te1, '\\|')[[1]]
  te.info2 = strsplit(te2, '\\|')[[1]]
  p1_te = dotplot(seq1, seq2nt(te.seq[te1]), 15, 13) + 
    ggtitle(paste(n1, 'vs TE:', te.info1[c(9, 7)])) + 
    theme(plot.title = element_text(size = 9)) 
  p2_te = dotplot(seq2, seq2nt(te.seq[te2]), 15, 13) +
    ggtitle(paste(n2, 'vs TE:', te.info2[c(9, 7)])) + 
    theme(plot.title = element_text(size = 9)) 

  # ORFs
  
  orf.res1 = orfFinder(nt2seq(seq1))
  orf.res2 = orfFinder(nt2seq(seq2))
  
  p.orf1 = NULL
  p.orf2 = NULL
  
  if(!is.null(orf.res1$pos)){
    p.orf1 =  orfplot(orf.res1$pos) + ylim(c(0, 10)) + 
      geom_vline(xintercept=c(0, len1), color = 'black', linetype="dashed") +
      ggtitle(paste(n1, ':', sv1)) + theme(legend.position = "bottom") +
      theme(plot.title = element_text(size = 9))
  }
  
  if(!is.null(orf.res2$pos)){
    p.orf2 = orfplot(orf.res2$pos) + ylim(c(0, 10)) +
      geom_vline(xintercept=c(0, len2), color = 'black', linetype="dashed") +
      ggtitle(paste(n2, ':', sv2)) + theme(legend.position = "bottom") +
      theme(plot.title = element_text(size = 9))
  }
  
  # Combine
  pp.te.orf = invisible(ggarrange(p.orf1, p.orf2,
                                  p1_te, p2_te,
                                  nrow = 2, ncol = 2))
  
  return(list(pp.stat = pp.stat,
              pp.te.orf = pp.te.orf))
}

