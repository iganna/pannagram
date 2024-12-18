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
                      x1.lim = 0.1, x2.lim = 0.1, y.expand = 0.02) {
  # # Required libraries
  # if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
  # if (!requireNamespace("dendextend", quietly = TRUE)) install.packages("dendextend")
  # if (!requireNamespace("ggtree", quietly = TRUE)) install.packages("ggtree")
  # if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
  # 
  # library("ape")
  # library("dendextend")
  # library("ggtree")
  # library(ggplot2)
  # library(cowplot)

  # Root and reorder trees
  t1 <- rootOpt(t1)
  t2 <- rootOpt(t2)
  t1 <- reorderTree(t1, t2)
  t2 <- reorderTree(t2, t1)
  
  # Colors
  if(is.null(group_dict)){
    
    tip.names = getTipOrder(t1)
    
    group_dict = as.list(tip.names)
    names(group_dict) = tip.names
    
    n = length(tip.names)
    colors <- c("#D61355", "#FFD23F", "#00215E")
    gradient <- colorRampPalette(colors)
    color_dict <- gradient(n)
    names(color_dict) = tip.names
  }
  
  # Generate tanglegram
  T1 <- ggtree(t1, ladderize = FALSE) +
    theme_tree2(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    geom_tiplab(align = TRUE) +
    ggtree::geom_text2(aes(label = label,
                           subset = !is.na(as.numeric(label)) & as.numeric(label) >= 75 & as.numeric(label) <= 100),
                       size = 2, hjust = 1, vjust = -1.5, color = 'grey30') +
    scale_y_continuous(expand = c(y.expand, y.expand))
  
  
  T2 <- ggtree(t2, ladderize=F) +   
    theme_tree2(legend.position='none', plot.margin = unit(c(0,0,0,0),"cm")) +
    geom_tiplab(hjust =1, align = TRUE) +
    ggtree::geom_text2(aes(label=label, 
                           subset = !is.na(as.numeric(label)) & as.numeric(label) >=75 & as.numeric(label) <=100),
                       size = 2,
                       hjust = 1, 
                       vjust = -1.5, color = 'grey30') +
    scale_y_continuous(expand = c(y.expand, y.expand))
  
  # Add colors
  T1 <- ggtree::groupOTU(T1, group_dict) +
    aes(color = group) + scale_color_manual(values = color_dict)
  T2 <- ggtree::groupOTU(T2, group_dict) +
    aes(color = group) + scale_color_manual(values = color_dict)
  
  
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
  
  
  final_plot <- cowplot::plot_grid(T1 + xlim(c(0, x1.lim)) + ggtitle(title1),
                                   L1,
                                   T2 + scale_x_reverse(limits = c(x2.lim, 0)) + ggtitle(title2),
                                   nrow = 1, align = "hv")
  
  return(final_plot)
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
  d_ini <- sum(a_ini)
  
  for (i.node in (n+1):(2*n - 1)) {
    if (!(i.node %in% t1$edge[,1])) next
    
    t1_rot <- ape::rotate(t1, node = i.node)
    tip1_rot <- getTipOrder(t1_rot)
    
    indices_vec_rot <- match(tip1_rot, tip2)
    a_rot <- abs((1:n) - indices_vec_rot)
    d_rot <- sum(a_rot)
    
    if (d_rot < d_ini) {
      d_ini <- d_rot
      a_ini <- a_rot
      t1 <- t1_rot
    } else if (d_rot == d_ini) {
      tmp <- which(a_rot != a_ini)
      b.rot <- a_rot[tmp]
      b.ini <- a_ini[tmp]
      if (var(b.rot) < var(b.ini)) {
        d_ini <- d_rot
        a_ini <- a_rot
        t1 <- t1_rot
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
#' @examples
#' # Example usage:
#' tree <- ape::rtree(5)
#' tip_order <- getTipOrder(tree)
getTipOrder <- function(tree) {
  n <- length(tree$tip.label)
  x <- tree$edge[,2]
  x <- x[x <= n]
  return(tree$tip.label[x])
}