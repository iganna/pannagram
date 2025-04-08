phyplot <- function(){

}


#' Generate a Tanglegram for Comparing Two Phylogenetic Trees
#'
#' Creates a tanglegram visualization to compare the tip order and structure of two phylogenetic trees.
#'
#' @param t1 The First phylogenetic tree of class \code{phylo}.
#' @param t2 The Second phylogenetic tree of class \code{phylo}.
#' @param color_dict Optional. A named vector of colors for tip labels. Default is \code{NULL}, and colors will be generated automatically.
#' @param group_dict Optional. A named list mapping tip labels to groups. Default is \code{NULL}, and groups will be generated automatically.
#' @param title1 Optional. Title for the first tree plot. Default is \code{NULL}.
#' @param title2 Optional. Title for the second tree plot. Default is \code{NULL}.
#' @param x1.lim Horizontal limit for the first tree plot. Default is \code{0.1}.
#' @param x2.lim Horizontal limit for the second tree plot. Default is \code{0.1}.
#' @param y.expand A numeric value to control vertical spacing in the plots. Default is \code{0.02}.
#'
#' @return A combined plot showing both trees and the lines connecting matching tips.
#' @examples
#' # Example usage:
#' tree1 <- ape::rtree(10)
#' tree2 <- ape::rtree(10)
#' plot <- tanglplot(tree1, tree2, title1 = "Tree 1", title2 = "Tree 2")
#' print(plot)
#' @export

tanglplot <- function(t1, t2, color_dict=NULL, group_dict=NULL, title1=NULL, title2=NULL,
                      x1.lim = 0.1, x2.lim = 0.1, y.expand = 0.02, to.order = T, 
                      show.lab = T, show.algn = T,
                      return.parts = F) {
  # # Required libraries
  # if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
  # if (!requireNamespace("dendextend", quietly = TRUE)) inst-titleall.packages("dendextend")
  # if (!requireNamespace("ggtree", quietly = TRUE)) install.packages("ggtree")
  # if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
  #
  # library("ape")
  # library("dendextend")
  # library("ggtree")
  # library(ggplot2)
  # library(cowplot)

  t1$edge.length[t1$edge.length < 0] <- 0
  t2$edge.length[t2$edge.length < 0] <- 0
  
  # Root and reorder trees
  if (!is.rooted(t1)) t1 <- rootOpt(t1)
  if (!is.rooted(t2)) t2 <- rootOpt(t2)
  
  if(to.order){
    t1 <- reorderTree(t1, t2)
    t2 <- reorderTree(t2, t1)  
  }

  # Colors
  if(is.null(group_dict)){

    tip.names = intersect(getTipOrder(t1), getTipOrder(t2))

    group_dict = as.list(tip.names)
    names(group_dict) = tip.names

    n = length(tip.names)
    colors <- c("#D61355", "#FFB433", "#00215E")
    gradient <- colorRampPalette(colors)
    color_dict <- gradient(n)
    names(color_dict) = tip.names
  }

  # Generate tanglegram
  T1 <- ggtree(t1, ladderize = FALSE) +
    theme_minimal() +
    theme_tree2(legend.position = "none", 
                # plot.margin = unit(c(0, 0, 0, 0), "cm"),
                plot.margin = margin(ggplot2::theme_get()$plot.margin[1],
                       unit(0, "pt"),
                       ggplot2::theme_get()$plot.margin[3],
                       ggplot2::theme_get()$plot.margin[4])) +
    geom_tiplab(aes(label = ifelse(show.lab, label, NA)), hjust = 1, 
                align = show.algn, linetype = "solid", alpha = 0.2) +
    ggtree::geom_text2(aes(label = label,
                           subset = !is.na(as.numeric(label)) & as.numeric(label) >= 75 & as.numeric(label) <= 100),
                       size = 2, hjust = 1, vjust = -1.5, color = 'grey30') +
    scale_y_continuous(expand = c(y.expand, y.expand)) +
    scale_x_continuous(expand = c(0, 0))


  T2 <- ggtree(t2, ladderize=F) +
    theme_minimal() +
    theme_tree2(legend.position = "none", 
                # plot.margin = unit(c(0, 0, 0, 0), "cm"),
                plot.margin = margin(ggplot2::theme_get()$plot.margin[1],
                                     unit(0, "pt"),
                                     ggplot2::theme_get()$plot.margin[3],
                                     ggplot2::theme_get()$plot.margin[4])) +
    geom_tiplab(aes(label = ifelse(show.lab, label, NA)), hjust = 1, 
                align = show.algn, linetype = "solid", alpha = 0.2) +
    ggtree::geom_text2(aes(label=label,
                           subset = !is.na(as.numeric(label)) & as.numeric(label) >=75 & as.numeric(label) <=100),
                       size = 2,
                       hjust = 1,
                       vjust = -1.5, color = 'grey30') +
    scale_y_continuous(expand = c(y.expand, y.expand)) +
    scale_x_continuous(expand = c(0, 0))

  #Add colors
  T1 <- ggtree::groupOTU(T1, group_dict) #+ aes(color = group) + scale_color_manual(values = color_dict, na.value = 'black')
  T2 <- ggtree::groupOTU(T2, group_dict) #+ aes(color = group) + scale_color_manual(values = color_dict, na.value = 'black')


  d1 <- T1$data[T1$data$isTip,]
  d1$x[] <- 1
  d2 <- T2$data[T2$data$isTip,]
  d2$x[] <- 2

  TTcon <- rbind(d1, d2)

  L1 <- ggplot(TTcon, aes(x = x, y = y, colour = group, group = label)) +
    geom_line() +
    theme_void() +
    theme(legend.position = "none", plot.margin = unit(c(1, 0, 1, 0), "cm")) +
    scale_y_continuous(expand = c(y.expand, y.expand)) +
    scale_color_manual(values = color_dict)

  
  T1 = T1 + xlim(c(0, x1.lim)) + ggtitle(title1)
  T2 = T2 + scale_x_reverse(limits = c(x2.lim, 0)) + ggtitle(title2)
  

  if(return.parts){
    return(list("T1" = T1,"L1" = L1,"T2" = T2))
  } else {
    final_plot <- cowplot::plot_grid(T1,
                                     L1,
                                     T2,
                                     nrow = 1, align = "hv")  
    return(final_plot)
  }
  
}


#' Optimize Phylogenetic Tree Rooting
#'
#' Re-roots a phylogenetic tree at the node that minimizes the maximum distance
#' to any tip.
#'
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @return A re-rooted tree of class \code{phylo}.
#' @details The function identifies the optimal root based on distance
#' calculations and re-roots the tree using an appropriate outgroup.
#' @importFrom ape root
#' @examples
#' # Example usage:
#' tree <- ape::rtree(10)
#' optimized_tree <- rootOpt(tree)
rootOpt <- function(tree) {
  n <- length(tree$tip.label)
  mx <- dist.nodes(tree)
  mx <- mx[1:n, (n+1):ncol(mx)]
  d.max <- apply(mx, 2, max)
  node.root <- n + which(d.max == min(d.max))

  get_descendants <- function(tree, node) {
    if (node > Ntip(tree)) {
      children <- list()
      descendants <- function(node) {
        if (!is.null(tree$edge)) {
          for (child in tree$edge[tree$edge[, 1] == node, 2]) {
            children <<- c(children, child)
            descendants(child)
          }
        }
      }
      descendants(node)
      return(children)
    } else {
      stop("Provided node is a tip, not an internal node.")
    }
  }

  children <- unlist(get_descendants(tree, node.root))
  children <- sort(children[children <= n])
  if (length(children) > n/2) {
    children <- setdiff(1:n, children)
  }

  tree_root <- ape::root(tree, outgroup = children, resolve.root = TRUE)
  return(tree_root)
}


#' Reorder Phylogenetic Tree
#'
#' Reorders a phylogenetic tree (\code{t1}) to match the tip order of a reference tree (\code{t2}),
#' optimizing for the smallest absolute differences in tip indices.
#'
#' @param t1 A phylogenetic tree of class \code{phylo} to be reordered.
#' @param t2 A reference phylogenetic tree of class \code{phylo}.
#' @return A reordered phylogenetic tree of class \code{phylo}, with a tip order
#'   that closely matches \code{t2}.
#' @details The function iteratively rotates internal nodes in \code{t1} to minimize
#' the distance between tip indices of \code{t1} and \code{t2}, based on their order.
#' Variance is used as a secondary criterion when distances are equal.
#' @importFrom ape rotate
#' @examples
#' # Example usage:
#' tree1 <- ape::rtree(10)
#' tree2 <- ape::rtree(10)
#' reordered_tree <- reorderTree(tree1, tree2)
reorderTree <- function(t1, t2) {
  n <- length(t1$tip.label)

  tip2 <- getTipOrder(t2)
  tip1 <- getTipOrder(t1)

  indices_vec_ini <- match(tip1, tip2)
  a_ini <- abs((1:n) - indices_vec_ini)
  d_ini <- sum(a_ini, na.rm = T)

  for (i.node in (n+1):(2*n - 1)) {
    if (!(i.node %in% t1$edge[,1])) next

    # t1_rot <- ape::rotate(t1, node = i.node)
    t1_rot <- rotateNode(t1, i.node)
    tip1_rot <- getTipOrder(t1_rot)
    
    # ggtree(t1_rot)
    # ggtree(t1_rot, ladderize = F)
    # ggtree(t1, ladderize = F)

    indices_vec_rot <- match(tip1_rot, tip2)
    a_rot <- abs((1:n) - indices_vec_rot)
    d_rot <- sum(a_rot,na.rm = T)

    if (d_rot < d_ini) {
      d_ini <- d_rot
      a_ini <- a_rot
      t1 <- t1_rot
    } else if (d_rot == d_ini) {
      tmp <- which(a_rot != a_ini)
      if(length(tmp) >= 2){
        b.rot <- a_rot[tmp]
        b.ini <- a_ini[tmp]
        if (var(b.rot) < var(b.ini)) {
          d_ini <- d_rot
          a_ini <- a_rot
          t1 <- t1_rot
        }        
      }
    }
  }
  return(t1)
}

#' Get Names of Tips in the Correct Order
#'
#' Extracts the order of tip labels from a phylogenetic tree based on their
#' appearance in the tree's edge matrix.
#'
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @return A character vector of tip labels in the order they appear in the tree's edge matrix.
#' @export
getTipOrder <- function(tree) {
  n <- length(tree$tip.label)
  x <- tree$edge[,2]
  x <- x[x <= n]
  return(tree$tip.label[x])
}


rotateNode <- function(tree, i.node){
  node.queue = i.node
  i.q = 1
  tips.ids = c()
  idx.edges = c()
  while(i.q <= length(node.queue)){
    node.curr.sign = node.queue[i.q] 
    node.curr = abs(node.curr.sign)
    i.q = i.q + 1
    
    id.next = which(tree$edge[,1] == node.curr)
    if(length(id.next) == 0){
      next
    }
    if(length(id.next) == 1) stop()
    
    if(node.curr == i.node){
      node.queue = c(node.queue, tree$edge[id.next,2] * c(-1, 1))  
      idx.edges = c(idx.edges, id.next * c(-1, 1))
    } else {
      node.queue = c(node.queue, tree$edge[id.next,2] * sign(node.curr.sign))  
      idx.edges = c(idx.edges, id.next * sign(node.curr.sign))
    }
  }
  
  edges = tree$edge
  edges = cbind(edges, 1:nrow(edges))
  
  id1 = sort(idx.edges[idx.edges > 0])
  id2 = sort(abs(idx.edges[idx.edges < 0]))
  
  id1.new = min(id2) + (1:length(id1)) / (length(id1) + 1)
  id2.new = min(id1) + (1:length(id2)) / (length(id2) + 1)
  
  edges[id1,3] = id1.new
  edges[id2,3] = id2.new
  
  edge.order = order(edges[,3])
  
  tree$edge = tree$edge[edge.order,]
  tree$edge.length = tree$edge.length[edge.order]
  
  return(tree)
}



