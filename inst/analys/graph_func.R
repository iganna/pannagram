# The script that was previously called as graph_refinement.Rg

#' Create and visualize a network graph from an igraph object
#'
#' @param g An igraph object representing the graph.
#' @param label Logical value for node labeling (default is TRUE).
#'
#' @return A ggplot2 object of the visualized network graph.
#'
ggigraph <- function(g, label = F, size = 5){
  edges.linear = get.edgelist(g)
  g.part <- network(edges.linear, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
  set.seed(20)
  p <- ggnet2(g.part, label = label, edge.color = "black",
              size = size, color = '#468B97',
              arrow.gap = 0.04, arrow.size = 5,
              # mode = "kamadakawai"
  )
  return(p)
}

#' Create a network graph and visualize it using ggnet2
#'
#' @param edges.linear A list of edges in the form of a matrix (n x 2).
#' @param label Logical value indicating whether to label the nodes (default is TRUE).
#'
#' @return A ggplot2 object representing the visualized network graph.
ggedges <- function(edges.linear, label = F, size = 5){
  g.part <- network(edges.linear, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
  set.seed(20)
  p <- ggnet2(g.part, label = label, edge.color = "black",
              size = size, color = '#468B97',
              arrow.gap = 0.04, arrow.size = 5,
              # mode = "kamadakawai"
  )
  return(p)
}


#' Give traits to nodes based on traits of sequences
#'
#' This function measures cumulative traits for nodes based sequences withint the node.
#' If traits are not unique - returns "Mix"
#'
#' @param nodes A data frame containing node information.
#' @param seqs.trait A vector containing the traits for sequences. 
#' Names of sequences sould be correct!
#' @param explain.mix An indicator to extend "Mix" annotation
#' @param mode A character string specifying the mode of operation. Options include:
#' \itemize{
#'   \item \strong{"unique"}: Returns unique trait values. If there are multiple unique values, it will either label as "Mix" or provide detailed mixed traits based on \code{explain.mix}.
#'   \item \strong{"mean"}: Returns the mean of trait sequences.
#'   \item \strong{"max"}: Returns the maximum value from trait sequences.
#'   \item \strong{"min"}: Returns the minimum value from trait sequences.
#' }
#' Other modes can be added as needed.
#' 
#' @return A vector containing the trait for each node.
#' 
#' @export
traitsSeqToNode <- function(nodes, seqs.trait, explain.mix = F, 
                            mode = 'maxunique',
                            mix.word = 'Mix',
                            mix.word.explain = 'Mix:'){
  seqs.trait = seqs.trait[nodes$name]
  if (mode == "unique") {
    nodes.trait <- tapply(seqs.trait, nodes$node, function(x) {
      return(paste(unique(x), collapse = ','))
    })
    if(explain.mix){
      nodes.trait[grepl(",", nodes.trait)] <- paste(mix.word.explain, sort(nodes.trait[grepl(",", nodes.trait)]), sep = '')
    } else {
      nodes.trait[grepl(",", nodes.trait)] <- mix.word
    }
  } else if (mode == "maxunique") {
    nodes.trait <- tapply(as.character(seqs.trait), nodes$node, function(x) {
      x = x[!is.na(x)]
      x = unique(x)
      if(length(x) == 0){
        return(NA)
      }
      else if(length(x) == 1){
        return(x)
      } else {
        x = table(x)
        x = names(x)[x == max(x)]
        return(x[1])
      }
    })
    
  } else if (mode == "mean") {
    nodes.trait <- tapply(seqs.trait, nodes$node, function(x) {
      return(mean(x))
    })
  } else if (mode == "max") {
    nodes.trait <- tapply(seqs.trait, nodes$node, function(x) {
      return(max(x))
    })
  } else if (mode == "min") {
    nodes.trait <- tapply(seqs.trait, nodes$node, function(x) {
      return(min(x))
    })
  } else if (mode == "cnt") {
    nodes.trait <- tapply(seqs.trait, nodes$node, function(x) {
      return(sum(x))
    })
  } 
  else {
    stop(paste("Invalid mode selected:", mode))
  }
  
  return(nodes.trait)
}

#' Get connected components of a directed graph
#'
#' This function takes a set of directed edges and calculates the connected
#' components in the corresponding directed graph.
#'
#' @param edges A matrix representing directed edges.
#'
#' @return A list of connected components, where each component is represented
#' as a vector of vertex indices.
#' 
#' @export
getGraphComponents <- function(edges){
  g <- igraph::simplify(igraph::make_graph(t(edges), directed = T))
  g.comp <- igraph::components(g)
  g.comp$cnt = table(g.comp$membership)
  return(g.comp)
}


#' Get Louvain Partition of a Graph
#'
#' This function takes an edge list, creates an undirected igraph object,
#' performs Louvain community detection, and returns a vector indicating the
#' community membership of each node.
#'
#' @param edges A data frame or matrix with two columns representing the edge list.
#' Each row corresponds to an edge between two nodes.
#' @param directed Logical, whether the input graph should be treated as directed.
#' Default is \code{TRUE}, but the graph will be converted to undirected internally.
#'
#' @return An integer vector of community memberships, where names correspond to node labels.
#' @examples
#' \dontrun{
#' library(igraph)
#' edges <- data.frame(
#'   from = c("A", "A", "B", "C"),
#'   to   = c("B", "C", "C", "D")
#' )
#' partition <- getGraphCommunities(edges)
#' print(partition)
#' }
#' @export
getGraphCommunities <- function(edges) {
  
  # igraph_g <- igraph::graph_from_edgelist(as.matrix(edges), directed = TRUE)
  # igraph_g <- igraph::as_undirected(igraph_g, mode = "collapse")
  
  igraph_g <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)
  
  louvain_result <- igraph::cluster_louvain(igraph_g) # Perform Louvain community detection
  
  return(igraph::membership(louvain_result))
}


#' Generate Graph Edges from Nestedness Results
#'
#' This function extracts directed edges between query and target names
#' based on a similarity cutoff from nestedness analysis results.
#'
#' @param nestedness A list or data frame containing \code{name.query},
#'   \code{name.target}, \code{coverage.query}, and \code{coverage.target}.
#' @param cov.cutoff Numeric value specifying the similarity cutoff.
#'   Values between 0–1 or 10–100 (interpreted as %) are accepted.
#'
#' @return A two-column matrix of unique edges representing connections
#'   that meet the similarity threshold.
#' @export
getGraphFromNestedness <- function(nestedness,
                                   cov.cutoff=NULL){
  if(is.null(cov.cutoff)){
    cov.cutoff = round(100 * min(rowMax(nestedness.major[,c('coverage.query', 'coverage.target')])))
    pokazAttention("Coverage for edges was set to", cov.cutoff, "\n",
                   "If you want to use a larger cov.cutoff, please provide the parameter.")
    # 
    # idx.1.to.2 = nestedness$coverage.query <= nestedness$coverage.target
    # edges = cbind(nestedness$name.query[idx.1.to.2], nestedness$name.target[idx.1.to.2])
    # idx.2.to.1 = nestedness$coverage.query >= nestedness$coverage.target
    # edges = rbind(edges, cbind(nestedness$name.target[idx.2.to.1], nestedness$name.query[idx.2.to.1]))
    # 
    # edges = unique(edges)
    # return(edges) 
  }
  
  cov.cutoff <- normalizeСovСutoff(cov.cutoff)
  
  # idx.1.to.2 = nestedness$coverage.query >= cov.cutoff
  idx.1.to.2 = (nestedness$coverage.query >= cov.cutoff) & 
    (nestedness$length.target / nestedness$length.query >= cov.cutoff)
  edges = cbind(nestedness$name.query[idx.1.to.2], nestedness$name.target[idx.1.to.2])
  # idx.2.to.1 = nestedness$coverage.target >= cov.cutoff
  idx.2.to.1 = (nestedness$coverage.target >= cov.cutoff) & 
    (nestedness$length.query / nestedness$length.target >= cov.cutoff)
  edges = rbind(edges, cbind(nestedness$name.target[idx.2.to.1], nestedness$name.query[idx.2.to.1]))
  
  edges = unique(edges)
  
  return(edges)
}

#' Create a Compact Representation of a Directed Graph
#'
#' This function constructs a compact version of a directed graph by identifying 
#' strongly connected components and grouping their nodes. It returns a simplified 
#' representation with compacted edges, node information, and component details.
#'
#' @param edges A data frame or matrix with two columns representing directed edges 
#'   (source and target nodes). Column 1 should contain the source nodes, and 
#'   column 2 should contain the target nodes.
#' @param echo Logical; if \code{TRUE}, prints progress information (currently unused 
#'   in the function but retained for compatibility).
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{edges}}{A unique matrix of edges between compacted nodes (no self-loops).}
#'   \item{\code{nodes}}{A data frame mapping each original node to its compact node label. 
#'     Nodes in strongly connected components are labeled as \code{"N#"}, 
#'     and all other nodes as \code{"R#"}.}
#'   \item{\code{nodes.list}}{A named list mapping each compact node ID to the original 
#'     node names it contains.}
#'   \item{\code{nodes.size}}{A named numeric vector giving the number of original nodes 
#'     in each compact node.}
#' }
#' @export
getGraphCompact <- function(edges, 
                            echo = T){
  
  g  <- graph_from_data_frame(edges, directed = TRUE)
  compsStrong <- components(g, mode = "strong")
  clustersStrons = which(compsStrong$csize != 1)
  names.mutual = names(compsStrong$membership)[compsStrong$membership %in% clustersStrons]
  names.rest =  names(compsStrong$membership)[!(compsStrong$membership %in% clustersStrons)]
  
  if(length(names.mutual) > 0){
    nodes.mutual = data.frame(node = paste0('N', compsStrong$membership[names.mutual]), 
                              name = names(compsStrong$membership[names.mutual]))  
  } else {
    nodes.mutual <- data.frame(node = character(), name = character(), 
                               stringsAsFactors = FALSE)  
  }
  if(length(names.rest) > 0){
    nodes.rest = data.frame(node = paste0('R', compsStrong$membership[names.rest]), 
                            name = names(compsStrong$membership[names.rest]))  
  } else {
    nodes.rest <- data.frame(node = character(), name = character(), 
                             stringsAsFactors = FALSE)  
  }
  
  nodes.info = rbind(nodes.mutual, nodes.rest)
  rownames(nodes.info) = nodes.info$name
  
  edges.compact = cbind(nodes.info[edges[,1],]$node,
                        nodes.info[edges[,2],]$node)
  edges.compact <- unique(edges.compact)
  edges.compact = edges.compact[edges.compact[,1] != edges.compact[,2],,drop=F]
  
  
  nodes.list <- split(nodes.info$name, nodes.info$node)
  table(nodes.info$node)
  nodes.size <- unlist(lapply(nodes.list, length))
  
  res.list = list(edges = edges.compact, 
                  nodes = nodes.info, 
                  nodes.list = nodes.list,
                  nodes.size = nodes.size)
  return(res.list)
  
}

#' Get indices of bypassed edges
#'
#' Finds edges that can be bypassed (i.e., there exists an alternative 2-step path).
#'
#' @param edges.compact Two-column matrix or data.frame of edges (from, to).
#' @return Integer vector of indices of bypassed edges.
#' @examples
#' edges <- matrix(c("A","B","B","C","A","C"), ncol=2, byrow=TRUE)
#' getEdgesIdxShortcut(edges)
#' @export
getEdgesIdxShortcut <- function(edges.compact) {
  # Get unique node names
  node.names <- unique(c(edges.compact))
  node.n <- length(node.names)
  
  # Map nodes to indices
  id1 <- match(edges.compact[,1], node.names)
  id2 <- match(edges.compact[,2], node.names)
  
  # Build sparse adjacency matrix
  A <- Matrix::sparseMatrix(i = id1, j = id2, x = TRUE, dims = c(node.n, node.n))
  
  # Get paths of length 2
  A2 <- A %*% A
  
  # Find bypassed edges (edges that have an alternative 2-step path)
  # A.remove <- (A2 & A)
  A.remove <- (A2 != 0) & (A != 0)
  # head(A.remove)
  # idx.remove <- which(A.remove, arr.ind = TRUE)
  nz <- Matrix::which(A.remove)                  # positions of TRUE as linear indices
  idx.remove <- arrayInd(nz, .dim = dim(A.remove))
  
  # Vectorized matching using interaction keys
  key.all <- interaction(id1, id2, drop = TRUE)
  key.remove <- interaction(idx.remove[,1], idx.remove[,2], drop = TRUE)
  
  # Find indices of edges that are bypassed
  idx.bypassed <- which(key.all %in% key.remove)
  
  return(idx.bypassed)
}



#' Remove direct edges from a graph. WORKING VERSION!!
#'
#' This function takes a list of edges and removes direct edges if a longer path exists.
#' Example list of edges:
#' (1) -> (2)
#' (2) -> (3)
#' (1) -> (3)
#' 
#' The edge (1) -> (3) should de removed, because path (1) -> (2) -> (3) exists.
#'
#' @param b.graph A list of edges in the form of a matrix (n x 2).
#' @param echo Logical value for printing progress (default is TRUE).
#'
#' @return A list of edges without direct edges
#'
refineDirectEdges <- function(edges.compact, echo = T){
  g = igraph::make_graph(t(edges.compact), directed = TRUE)
  igraph::E(g)$name <- paste(edges.compact[,1], edges.compact[,2], sep = '-')
  
  # Find nodes, which should be kept with all their edges and delete them from the graph
  nodes.keep = c()
  n.nodes.keep = -1
  while(length(nodes.keep) != n.nodes.keep){
    n.nodes.keep = length(nodes.keep)
    # print(n.nodes.keep)
    print(igraph::vcount(g))
    
    deg.in <- igraph::degree(g, mode = "in")
    deg.out <- igraph::degree(g, mode = "out")
    # tails = names(deg.in)[((deg.in + deg.out) == 1) & ((deg.in * deg.out) == 0) | ((deg.in * deg.out) == 1)]
    tails = names(deg.in)[((deg.in + deg.out) == 1) & ((deg.in * deg.out) == 0) ]
    
    nodes.keep = c(nodes.keep, tails)
    g <- igraph::delete_vertices(g, tails)
  }
  
  # Remove connected components with size 1 and keep node names
  g.comp <- igraph::components(g)
  id.single = which(g.comp$csize == 1)
  nodes.single = names(g.comp$membership)[g.comp$membership %in% id.single]
  nodes.keep = c(nodes.keep, nodes.single)
  g <- igraph::delete_vertices(g, nodes.single)
  
  # Refine each connected component
  g.comp <- igraph::components(g)
  edges.polised = c()
  for(i.comp in 1:g.comp$no){
    if(g.comp$no == 0) break
    # message(i.comp)
    
    # Get subgraph
    names.comp = names(g.comp$membership)[g.comp$membership == i.comp]
    g.sub <- induced_subgraph(g, names.comp)
    
    # print(c(vcount(g.sub), ecount(g.sub)))
    # ggigraph(g.sub)
    
    # Find all nodes, which have "out" degree != 0.
    # TODO: change to >= 2 and check.
    deg.out <- igraph::degree(g.sub, mode = "out")
    targets = names(deg.out)[deg.out != 0]  
    edges.sub =  get.edgelist(g.sub)
    
    # Find all destinations of these nodes: (target) -> (out1)
    idx.target = edges.sub[,1] %in% targets
    out1 = edges.sub[idx.target, 2]
    out1.comb = paste(edges.sub[idx.target, 1], edges.sub[idx.target, 2], sep = '_')  # create "target_out1"
    
    # Find second-layer destinations: (out1) -> (out2)
    idx.out1 = which(edges.sub[,1] %in% unique(out1))
    out2 = tapply(edges.sub[idx.out1,2], edges.sub[idx.out1,1], function(x) list(x))
    
    # Connect second-layer destinations to the initial targets: (target) -> (out1) -> (out2)
    out12 = out2[out1]
    names(out12) = paste0(edges.sub[idx.target, 1], '_', out1, '_')
    
    # Remove the intermediate subname: change "target_out1_out2" into "target_out2"
    out12 = unlist(out12)
    out12.pref = sub("_.*", "", names(out12))
    out2.comb = paste(out12.pref, out12, sep = '_')
    
    # Find edges, which exist twice..
    out.common = intersect(out1.comb, out2.comb)  # intersect "target_out1" and "target_out2"
    # .. and remove these edges
    if(length(out.common) != 0){
      edges.delete <- as.matrix(data.frame(do.call(rbind, strsplit(out.common, "_"))))
      
      id.edges.delete <- sapply(1:nrow(edges.delete), function(i) get.edge.ids(g.sub, edges.delete[i,]))
      g.sub <- delete_edges(g.sub, id.edges.delete)
    }
    # print(c(vcount(g.sub), ecount(g.sub)))
    
    # 
    # OLD CODE, not efficient, but still here to remember the fail state.
    # i = 1
    # for(target.vertex in targets){
    #   i = i + 1
    #   if(round(i/100) == i/100) print(i)
    #   out1 <- V(g.sub)[neighbors(g.sub, target.vertex, mode = "out")]$name
    #   out2 <- unlist(sapply(out1, function(s) V(g.sub)[neighbors(g.sub, s, mode = "out")]$name))
    #   v.remove = intersect(out1, out2)
    #   if(length(v.remove) == 0) next
    #   
    #   edge_to_delete <- sapply(v.remove, function(s) get.edge.ids(g.sub, c(target.vertex, s)))
    #   
    #   g.sub <- delete_edges(g.sub, edge_to_delete)  
    # }
    # ggigraph(g.sub, label = F)
    # ggigraph(g.sub)
    
    # Save
    edges.sub = get.edgelist(g.sub)
    edges.polised = rbind(edges.polised, edges.sub)
  }
  
  # Add all edges for nodes, which were kept before
  edges.add = edges.compact[(edges.compact[,1] %in% nodes.keep) | (edges.compact[,2] %in% nodes.keep), ]
  edges.final = rbind(edges.polised, edges.add)
  return(edges.final)
}


#' Filter a coverage matrix by various criteria
#'
#' Filters a coverage matrix (`nestedness`) based on coverage percentage, sequence 
#' length, and minimum copy number.
#'
#' @param nestedness Data frame with coverage and length information. Must contain 
#'   columns: `name.query`, `name.target`, 
#'            `length.query`, `length.target`,
#'            `coverage.query`, `coverage.target` .
#' @param cov.cutoff Integer. Minimum coverage percentage to retain entries (default = NULL).
#' @param min.len Integer. Minimum sequence length to retain entries (default = NULL).
#' @param max.len Integer. Maximum sequence length to retain entries (default = NULL).
#' @param min.copy Integer. Minimum number of copies to retain entries (default = NULL).
#'
#' @return A filtered data frame with the same structure as `nestedness`.
#' @examples
#' filtered <- filterNestedness(nestedness, cov.cutoff = 80, min.len = 500)
#' @export
filterNestedness <- function(nestedness,
                             cov.cutoff = NULL,
                             min.len    = NULL,
                             max.len    = NULL,
                             min.copy   = NULL,
                             names.remove = NULL,
                             show.echo = F) {
  
  if(show.echo){
    pokaz('Initial number of rows:', nrow(nestedness))
  }
  
  # ---- Filter by names to remove ----
  if(!is.null(names.remove)){
    nestedness = nestedness[!(nestedness$name.query %in% names.remove) & 
                              !(nestedness$name.target %in% names.remove),]
    
    if(show.echo){
      pokaz('Number of rows after name filtration:', nrow(nestedness))
    }
  }
  
  # ---- Filter by coverage cutoff (percentage) ----
  if (!is.null(cov.cutoff)) {
    cov.cutoff <- normalizeСovСutoff(cov.cutoff)
    nestedness <- nestedness[(nestedness$coverage.query  >= cov.cutoff) |
                               (nestedness$coverage.target >= cov.cutoff), ]
    
    if(show.echo){
      pokaz('Number of rows after filtration by coverage:', nrow(nestedness))
    }
  }
  
  # ---- Filter by minimum length ----
  if (!is.null(min.len)) {
    nestedness <- nestedness[(nestedness$length.query  >= min.len) &
                               (nestedness$length.target >= min.len), , drop = FALSE]
    
    if(show.echo){
      pokaz('Number of rows after filtration by minimum length:', nrow(nestedness))
    }
  }
  
  # ---- Filter by maximum length ----
  if (!is.null(max.len)) {
    nestedness <- nestedness[(nestedness$length.query  <= max.len) &
                               (nestedness$length.target <= max.len), , drop = FALSE]
    
    if(show.echo){
      pokaz('Number of rows after filtration by maximum length:', nrow(nestedness))
    }
  }
  
  # ---- Filter by minimum number of copies ----
  if (!is.null(min.copy)) {
    
    if (is.null(cov.cutoff)) {
      warning("'cov.cutoff' sould be provided to filter by copy-number")
    }
    
    cov.cutoff <- normalizeСovСutoff(cov.cutoff)
    
    nestedness.equal <- nestedness[(nestedness$coverage.query  >= cov.cutoff) &
                                     (nestedness$coverage.target >= cov.cutoff), , drop = FALSE]
    
    if (nrow(nestedness.equal) > 0) {
      
      if(min.copy == 2){
        sv.passed = c(nestedness.equal$name.query, nestedness.equal$name.target)
      } else {
        # getGraphComponents should return:
        #   $membership — vector of component IDs (names = vertices)
        #   $csize      — vector of component sizes (index = component ID)
        g.comp <- getGraphComponents(nestedness.equal[, c("name.query", "name.target")])
        
        comp.passed <- which(g.comp$csize >= min.copy)
        sv.passed   <- names(g.comp$membership)[g.comp$membership %in% comp.passed]
      }
      
      nestedness <- nestedness[(nestedness$name.query  %in% sv.passed) &
                                 (nestedness$name.target %in% sv.passed), , drop = FALSE]        
      
    } else {
      # Return an empty data frame with the same structure if no edges remain
      nestedness <- nestedness[0, , drop = FALSE]
    }
    if(show.echo){
      pokaz('Number of rows after filtration by copies:', nrow(nestedness))
    }
  }
  
  return(nestedness)
}


#' Filter graph edges by various criteria
#'
#' @description
#' Filters an edge matrix by:
#' * minimum component size,
#' * removing inter-community edges,
#' * removing bridge edges.
#'
#' @param edges Two-column edge list (matrix or data.frame).
#' @param min.comp.size Integer or NULL. Remove edges connected to components smaller than this size.
#' @param remove.intercommunity Logical. If TRUE, keep only edges inside the same community.
#' @param remove.bridges Logical. If TRUE, remove bridge edges.
#'
#' @return Filtered edge list (matrix).
#' @examples
#' \dontrun{
#' edges <- data.frame(from = c("a","a","b","c","d"),
#'                     to   = c("b","c","c","d","e"))
#' filterEdges(edges, min.comp.size = 3,
#'             remove.intercommunity = TRUE,
#'             remove.bridges = TRUE)
#' }
#' @importFrom igraph graph_from_edgelist simplify bridges delete_edges as_data_frame
#' @export
filterEdges <- function(edges,
                        min.comp.size = NULL,
                        remove.intercommunity = FALSE,
                        remove.bridges = FALSE,
                        show.echo = F) {
  
  # ---- Filter by component size ----
  if (!is.null(min.comp.size)) {
    g.comp <- getGraphComponents(edges)
    comp.remove <- which(g.comp$csize < min.comp.size)
    sv.remove <- names(g.comp$membership)[g.comp$membership %in% comp.remove]
    edges <- edges[!(edges[,1] %in% sv.remove) & !(edges[,2] %in% sv.remove), ]
  }
  
  # ---- Filter edges between communities ----
  if (remove.intercommunity) {
    communities <- getGraphCommunities(edges)
    edges <- edges[communities[edges[,1]] == communities[edges[,2]], ]
  }
  
  # ---- Filter bridges ----
  if (remove.bridges) {
    g <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)
    g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    bridge.ids <- igraph::bridges(g)
    g.no.bridges <- igraph::delete_edges(g, bridge.ids)
    edges.no.bridges <- as.matrix(igraph::as_data_frame(g.no.bridges, what = "edges"))
    
    sv.remain = c(edges.no.bridges[,1], edges.no.bridges[,2])
    
    if(show.echo){
      pokaz('Number of edges after bridge filtration (undirected):', nrow(edges.no.bridges))
    }
    
    edges = edges[(edges[,1] %in% sv.remain) & 
                    (edges[,2] %in% sv.remain),, drop= F]
    
    if(show.echo){
      pokaz('Number of edges after bridge filtration:', nrow(edges))
    }
  }
  
  return(edges)
}

#' @export
filterEdgesShortcut <- function(edges.major, show.echo = F){
  
  
  if(nrow(edges.major) == 0) stop('Edge matrix is empty')
  
  graph.compact = getGraphCompact(edges.major)
  edges.compact = graph.compact$edges
  
  if(nrow(edges.compact) == 0){
    if(show.echo) pokazAttention('The compact graph has no edges')
    return(edges.major)
  }
  
  # Components with two nodes are not needed to be filtered, therefore min.comp.size = 3
  edges.compact = filterEdges(edges.compact, min.comp.size = 3)
  
  if(nrow(edges.compact) == 0){
    if(show.echo) pokazAttention('The compact graph has no edges')
    return(edges.major)
  }
  
  idx.bypassed = getEdgesIdxShortcut(edges.compact)
  
  if(length(idx.bypassed) == 0){
    if(show.echo) pokazAttention('No baypassed Edges in the initial graph')
    return(edges.major)
  }
  
  # Remove corresponding edges from edges.major
  idx.edges.remove = c()
  for(i.edge in idx.bypassed){
    nodes.from = graph.compact$nodes.list[[edges.compact[i.edge, 1]]]
    nodes.to = graph.compact$nodes.list[[edges.compact[i.edge, 2]]]
    i.edges.remove = which((edges.major[,1] %in% nodes.from) & (edges.major[,2] %in% nodes.to))
    idx.edges.remove = c(idx.edges.remove, i.edges.remove)
  }
  
  if(show.echo) pokaz('Number of bypass edges is', length(idx.edges.remove))
  if(length(idx.edges.remove) > 0){
    comp.before = getGraphComponents(edges.major)
    edges.major = edges.major[-idx.edges.remove,]
    comp.after = getGraphComponents(edges.major)
    # Check number of components
    if(comp.before$no != comp.after$no){
      stop("Some deleted edges were crucial")
    }
  }
  
  return(edges.major)
  
}

#' Restore missing SV edges to the graph
#'
#' Re-adds edges for structural variants (SVs) that were removed, assigning each
#' missing SV to the component where most of its neighbors belong (if dominance
#' > `dominant.effect`).
#'
#' @param edges Current edge list (matrix or data frame, 2 columns).
#' @param edges.init Initial full edge list.
#' @param echo Logical, print progress if TRUE. Default: FALSE.
#' @param dominant.effect Numeric threshold (0–1) for dominant component membership. Default: 0.7.
#'
#' @return Updated edge list with restored SV edges.
#' @export
putEdgesBack <- function(edges, edges.init, show.echo=F, dominant.effect = 0.7,
                         small.only = F){
  
  if(small.only){
    sv.back = setdiff(c(as.matrix(edges.init[,1])), c(as.matrix(edges)))
  } else {
    sv.back = setdiff(c(as.matrix(edges.init)), c(as.matrix(edges)))  
  }
  
  sv.len = parseStrings(unique(c(as.matrix(edges.init))), 2, numeric = T)
  if(length(sv.back) == 0){
    if(show.echo) pokaz('No SVs to put back', length(sv.back))
    return(edges)
  } else {
    if(show.echo) pokaz('Number of SVs to put back', length(sv.back))  
  }
  
  # Sort by length (from longest to shortest)
  sv.back = sv.back[order(- sv.len[sv.back])]
  
  edges.back = edges.init
  edges.back = edges.back[(edges.back[,1] %in% sv.back) | (edges.back[,2] %in% sv.back),]
  edges.back = as.matrix(edges.back)
  
  g.comp <- getGraphComponents(edges)
  partition <- g.comp$membership
  partition.add = partition
  
  if(!is.null(partition.add)){
    max.part = max(partition.add)  
  } else {
    max.part = 0
  }
  
  # Map of "in which Edges every SV play the role"
  idx.names = rbind(cbind(edges.back[,1], 1:nrow(edges.back)),
                    cbind(edges.back[,2], 1:nrow(edges.back)))
  idx.names = idx.names[idx.names[,1] %in% sv.back,]
  idx.map = split(idx.names[,2], idx.names[,1])
  
  sv.edges.add.list = list() 
  sv.edges.add.counter = 1  
  
  for(s.sv in sv.back){
    
    idx = as.numeric(idx.map[[s.sv]])
    if (is.null(idx)) {
      stop()
      next 
    }
    
    # Get neighbours
    sv.edges.connect = edges.back[idx,,drop=F]
    
    if(nrow(sv.edges.connect) == 0){
      stop('WHY')
      next
    } 
    
    sv.node.connect = unique(c(sv.edges.connect))
    sv.node.connect = sv.node.connect[sv.node.connect != s.sv]
    sv.node.connect = intersect(sv.node.connect, names(partition.add))
    
    if(length(sv.node.connect) == 0) {
      partition.add[s.sv] = max.part + 1
      max.part = max.part + 1
      next
    }
    # stop()
    sv.node.connect.part = partition.add[sv.node.connect]
    sv.node.connect.part = sv.node.connect.part[!is.na(sv.node.connect.part)]
    part.cnt = table(sv.node.connect.part)
    
    i.part = NA  # In case one want to remove the "else" branch in the following if
    
    if(max(part.cnt)/sum(part.cnt) > dominant.effect){
      i.part = as.numeric(names(part.cnt)[which.max(part.cnt)])
    } 
    # else {
    #   sv.node.connect.part.len = sv.len[names(sv.node.connect.part)]
    #   s.sv.len = sv.len[s.sv]
    #   idx.max.len = which.max(sv.node.connect.part.len)
    #   if(sv.node.connect.part.len[idx.max.len] / s.sv.len > 0.7){
    #     i.part = sv.node.connect.part[idx.max.len]  
    #   }
    # }
    
    if(!is.na(i.part)){
      sv.part = names(partition.add)[partition.add == i.part]
      sv.part = c(sv.part, s.sv)
      sv.edges.add = sv.edges.connect[(sv.edges.connect[,1] %in% sv.part) & 
                                        (sv.edges.connect[,2] %in% sv.part), , drop=F]
      
      # sv.edges.add = sv.edges.add[sv.edges.add[,1] == s.sv,,drop=F]  # Only s.sv -> something
      # pokaz(1, nrow(sv.edges.add))
      # sv.edges.add = sv.edges.add[sv.edges.add[,2] %in% names(partition.add),,drop=F]
      # pokaz(2, nrow(sv.edges.add))
      
      # if(s.sv %in% 'SVgr_2_id_249548|16643'){
      #   print(sv.edges.add)
      #   stop()
      # }
      
      if (nrow(sv.edges.add) > 0) {
        
        partition.add[s.sv] = i.part
        sv.edges.add.list[[sv.edges.add.counter]] = sv.edges.add
        sv.edges.add.counter = sv.edges.add.counter + 1
        
      } else {
        partition.add[s.sv] = 0
        #pokaz('No edges to add')
        # stop('Something is wrong')
      }
    } 
  }
  partition.add = partition.add[partition.add != 0]
  
  if(show.echo) pokaz('Number of edges to put back:', sv.edges.add.counter)
  
  if (length(sv.edges.add.list) > 0) {
    edges.add = do.call(rbind, sv.edges.add.list)
  } else{
    return(edges)
  }
  
  edges = rbind(edges, edges.add)
  
  # if(show.echo) pokaz('Remove shortcuts...')
  # edges <- filterEdgesShortcut(edges)
  
  return(rbind(edges, edges.add))
  
}

#' @export
solveForkNodes <- function(edges,
                             seqs,
                             cutoff.remain.edges = 0.7,
                             flank.cover.cutoff = 0.8,
                             show.echo = FALSE)
{
  
  graph.compact <- getGraphCompact(edges)
  edges.compact <- graph.compact$edges
  
  if (nrow(edges.compact) == 0) {
    return(edges)
  }
  
  comp.before <- NULL
  comp.after  <- NULL
  
  # Consider all nodes that have at least two outgoing edges
  stat.neighbours.all <- NULL
  nodes.parasite <- unique(edges.compact[duplicated(edges.compact[, 1]), 1])
  
  for (node.from in nodes.parasite) {
    stat.neighbours <- data.frame(edge.id = which(edges.compact[, 1] == node.from))
    
    stat.neighbours$node.from <- node.from
    stat.neighbours$node.to   <- edges.compact[stat.neighbours$edge.id, 2]
    stat.neighbours$size.from <- graph.compact$nodes.size[node.from]
    stat.neighbours$size.to   <- graph.compact$nodes.size[stat.neighbours$node.to]
    
    edges.nei <- edges[edges[, 1] %in% graph.compact$nodes.list[[node.from]], ]
    
    # Find how many are connected on the "from"-side
    edges.nei.mod <- edges.nei
    edges.nei.mod[, 2] <- graph.compact$nodes[edges.nei[, 2], ]$node
    edges.nei.mod <- unique(edges.nei.mod)
    
    connect.from     <- split(edges.nei.mod[, 1], edges.nei.mod[, 2])
    connect.from.len <- unlist(lapply(connect.from, length))
    # stat.neighbours$n.from <- connect.from.len[stat.neighbours$node.from]
    stat.neighbours$n.from <- connect.from.len[stat.neighbours$node.to]
    
    # Find how many are connected on the "to"-side
    names.to <- unique(edges.nei[, 2])
    nodes.to <- graph.compact$nodes[names.to, ]$node
    
    cnt.nodes.to <- table(nodes.to)
    stat.neighbours$n.to <- cnt.nodes.to[stat.neighbours$node.to]
    
    if (sum(is.na(stat.neighbours)) > 0)
      stop("Wrong names of neighbours.")
    
    # Compute proportions
    stat.neighbours$p.from <- stat.neighbours$n.from / stat.neighbours$size.from
    stat.neighbours$p.to   <- stat.neighbours$n.to   / stat.neighbours$size.to
    
    stat.neighbours$remain <-
      (stat.neighbours$p.to   >= cutoff.remain.edges) &
      (stat.neighbours$p.from >= cutoff.remain.edges)
    
    # Check whether the sequences are intersect by flanking regions
    
    edges.nei.mod <- edges.nei
    edges.nei.mod[, 2] <- graph.compact$nodes[edges.nei[, 2], ]$node
    
    for (irow in which(stat.neighbours$remain)) {
      node.to <- stat.neighbours$node.to[irow]
      
      idx.tmp <- which(edges.nei.mod[, 2] == node.to)[1]
      edge.tmp <- edges.nei[idx.tmp, ]
      
      s1 <- seq2nt(seqs[edge.tmp[1]])
      s2 <- seq2nt(seqs[edge.tmp[2]])
      
      n.cut <- min(round(length(s1) / flank.cover.cutoff), length(s2))
      
      score.tot <- scoreFlankCoverage(s1, s2, n.cut, 15, 12)
      
      if (score.tot < flank.cover.cutoff) {
        stat.neighbours$remain[irow] <- FALSE
      }
    }
    
    # # Remain the shortest?
    # if(sum(stat.neighbours$remain) > 1){
    #   stat.neighbours$len.to[!stat.neighbours$remain] = 0
    #   irow.shortest = which.min(stat.neighbours$len.to)
    #   
    #   stat.neighbours$remain = F
    #   stat.neighbours$remain[irow.shortest] = T
    # }
    
    stat.neighbours.all <- rbind(stat.neighbours.all, stat.neighbours)
  }
  
  if(!is.null(stat.neighbours.all)){
    idx.edge.remove = stat.neighbours.all$edge.id[!stat.neighbours.all$remain]  
  } else {
    idx.edge.remove = c()
  }
  
  # Remove corresponding edges from edges
  if(length(idx.edge.remove) > 0){
    idx.edges.remove = c()
    for(i.edge in idx.edge.remove){
      nodes.from = graph.compact$nodes.list[[edges.compact[i.edge, 1]]]
      nodes.to = graph.compact$nodes.list[[edges.compact[i.edge, 2]]]
      i.edges.remove = which((edges[,1] %in% nodes.from) & (edges[,2] %in% nodes.to))
      idx.edges.remove = c(idx.edges.remove, i.edges.remove)
    }
    
    if(show.echo) pokaz('Number of edges to remove', length(idx.edges.remove))
    if(length(idx.edges.remove) > 0){
      comp.before = getGraphComponents(edges)
      edges = edges[-idx.edges.remove,]
      comp.after = getGraphComponents(edges)
      # Check number of components
    }
  }
  
  return(edges)
}

#' @export
solveUmbrellaNodes <- function(edges,
                                seqs, 
                                coverage.umbrella.children = 100,
                                cutoff.remain.edges = 0.7,
                                flank.cover.cutoff = 0.8,
                                show.echo = FALSE){
  
  graph.compact = getGraphCompact(edges)
  edges.compact = graph.compact$edges
  
  if(nrow(edges.compact) > 0){
    stat.neighbours.all = c()
    nodes.umbrella = setdiff(unique(edges.compact[duplicated(edges.compact[,2]),2]),
                             edges.compact[,1])
    for(node.to in nodes.umbrella){
      stat.neighbours = data.frame(edge.id = which(edges.compact[,2] == node.to))
      
      stat.neighbours$node.from = edges.compact[stat.neighbours$edge.id,1]
      stat.neighbours$node.to = node.to
      stat.neighbours$size.from = graph.compact$nodes.size[stat.neighbours$node.from] 
      stat.neighbours$size.to = graph.compact$nodes.size[stat.neighbours$node.to] 
      
      edges.nei = edges[edges[,2] %in% graph.compact$nodes.list[[node.to]], ]
      
      # Find how many are connected on the "to"-side
      edges.nei.mod = edges.nei
      edges.nei.mod[,1] = graph.compact$nodes[edges.nei[,1],]$node
      edges.nei.mod = unique(edges.nei.mod)
      
      connect.to = split(edges.nei.mod[,2], edges.nei.mod[,1])
      connect.to.len = unlist(lapply(connect.to, length))
      stat.neighbours$n.to = connect.to.len[stat.neighbours$node.from]
      
      # Find how many are connected on the "from"-side
      names.from = unique(edges.nei[,1])
      nodes.from = graph.compact$nodes[names.from,]$node
      
      cnt.nodes.from = table(nodes.from)
      stat.neighbours$n.from = cnt.nodes.from[stat.neighbours$node.from]
      if(sum(is.na(stat.neighbours)) > 0) stop('Wrong names of neighbours.')
      
      # Compute proportions
      stat.neighbours$p.from = stat.neighbours$n.from / stat.neighbours$size.from
      stat.neighbours$p.to = stat.neighbours$n.to / stat.neighbours$size.to
      
      stat.neighbours$remain = (stat.neighbours$p.to >= cutoff.remain.edges) & 
        (stat.neighbours$p.from >= cutoff.remain.edges)
      
      # Check whether the sequences are intersect by flanking regions
      
      edges.nei.mod = edges.nei
      edges.nei.mod[,1] = graph.compact$nodes[edges.nei[,1],]$node
      stat.neighbours$len.from = 0
      stat.neighbours$len.to = 0
      for(irow in which(stat.neighbours$remain)){
        node.from = stat.neighbours$node.from[irow]
        # stop()
        
        idx.tmp = which(edges.nei.mod[,1] == node.from)[1]
        edge.tmp = edges.nei[idx.tmp,]
        
        s1 = seq2nt(seqs[edge.tmp[1]])
        s2 = seq2nt(seqs[edge.tmp[2]])
        
        stat.neighbours$len.from[irow] = length(s1)
        stat.neighbours$len.to[irow] = length(s2)
        
        n.cut = min(round(length(s1) / flank.cover.cutoff), length(s2))
        
        score.tot <- scoreFlankCoverage(s1, s2, n.cut, 15, 12)
        
        if(score.tot < flank.cover.cutoff){
          stat.neighbours$remain[irow] = F
        }
      }
      
      # # Remain the longest and those which match with the longest
      # if(sum(stat.neighbours$remain) > 1){
      #   stat.neighbours$len.from[!stat.neighbours$remain] = 0
      #   irow.longest = which.max(stat.neighbours$len.from)
      #   
      #   stat.neighbours$remain = F
      #   stat.neighbours$remain[irow.longest] = T
      # }
      
      # Keep the results in the common dataframe
      stat.neighbours.all = rbind(stat.neighbours.all, stat.neighbours)
      # stop()
    }
    
    if(!is.null(stat.neighbours.all)){
      idx.edge.remove = stat.neighbours.all$edge.id[!stat.neighbours.all$remain]
    } else {
      idx.edge.remove = c()
    }
    
    # Remove corresponding edges from edges
    if(length(idx.edge.remove) > 0){
      idx.edges.remove = c()
      for(i.edge in idx.edge.remove){
        nodes.from = graph.compact$nodes.list[[edges.compact[i.edge, 1]]]
        nodes.to = graph.compact$nodes.list[[edges.compact[i.edge, 2]]]
        i.edges.remove = which((edges[,1] %in% nodes.from) & (edges[,2] %in% nodes.to))
        idx.edges.remove = c(idx.edges.remove, i.edges.remove)
      }
      
      if(show.echo) pokaz('Number of edges to remove', length(idx.edges.remove))
      if(length(idx.edges.remove) > 0){
        comp.before = getGraphComponents(edges)
        edges = edges[-idx.edges.remove,]
        comp.after = getGraphComponents(edges)
      }
    }
  }
  
  return(edges)
}

#' @export
normalizeСovСutoff <- function(cov.cutoff) {
  if (cov.cutoff >= 10 && cov.cutoff <= 100) {
    cov.cutoff <- cov.cutoff / 100
  } else if (cov.cutoff < 0 || cov.cutoff > 1) {
    stop("cov.cutoff does not make sense")
  }
  return(cov.cutoff)
}

#' @export
getComponentSequences <- function(seqs.names.comp, 
                                  seqs, 
                                  nestedness){
  
  seqs.names.comp = seqs.names.comp[seqs.names.comp %in% names(seqs)]
  if(length(seqs.names.comp) == 0){
    stop('Names of sequences (seqs.names.comp) were not found among the provided sequences (seqs).') 
  }
  
  # Nestedness for the component
  nest.target = nestedness[(nestedness$name.query %in% seqs.names.comp) & (nestedness$name.target %in% seqs.names.comp),]
  
  # Get remained names
  seqs.names.comp = unique(c(nest.target$name.query, nest.target$name.target))
  if(length(seqs.names.comp) == 0){
    stop('Names of sequences (seqs.names.comp) were not found in nestedness matrix (nestedness).') 
  }
  
  # Get sequences and Lengths
  seqs.target = seqs[seqs.names.comp]
  seqs.target.len = as.numeric(parseStrings(seqs.names.comp, 2))
  
  # Sorting by length
  seqs.target = seqs.target[order(-seqs.target.len)]
  seqs.names.comp = names(seqs.target)
  
  # Find orientations of sequences
  orientation.target = rep('.', length(seqs.target))
  names(orientation.target) = seqs.names.comp
  
  # Define the first orientation by the longest ORF in the first sequence
  s = seqs[seqs.names.comp[1]]
  orf.res = orfFinder(s)
  orientation.target[1] = orf.res$pos$strand[1]
  
  dir.seq = c('-', '+')
  
  n.dot = 0
  while (n.dot != sum(orientation.target == '.')) {
    n.dot = sum(orientation.target == '.')
    for(i in 1:(length(orientation.target))){
      if(orientation.target[i] == '.') next
      nest.tmp = data.frame(name = c(nest.target[nest.target$name.query == seqs.names.comp[i],]$name.target,
                                     nest.target[nest.target$name.target == seqs.names.comp[i],]$name.query),
                            strand = c(nest.target[nest.target$name.query == seqs.names.comp[i],]$strand,
                                       nest.target[nest.target$name.target == seqs.names.comp[i],]$strand),
                            coverage = c(nest.target[nest.target$name.query == seqs.names.comp[i],]$coverage.query,
                                         nest.target[nest.target$name.target == seqs.names.comp[i],]$coverage.target))
        
      nest.tmp = nest.tmp[order(-nest.tmp$coverage),]
      nest.tmp = nest.tmp[order(nest.tmp$name),]
      nest.tmp = nest.tmp[!duplicated(nest.tmp$name),]
      
      nest.tmp = nest.tmp[orientation.target[nest.tmp$name] == '.',]
      if(nrow(nest.tmp) == 0) next
      
      for(irow in 1:nrow(nest.tmp)){
        if(nest.tmp$strand[irow] == '+'){
          orientation.target[nest.tmp$name[irow]] = orientation.target[i]
        } else {
          orientation.target[nest.tmp$name[irow]] = setdiff(dir.seq, orientation.target[i])
        }
      }
      if(sum(orientation.target == '.') == 0) break
    }
  }
  
  if(sum(orientation.target == '.') > 0){
    stop('Not all of the sequences were oriended')
  }
  # 
  for(i in which(orientation.target == '-')){
    seqs.target[i] = revComplSeq(seqs.target[i])
  }
  
  # seqs.target.msa = unlist(lapply(seqs.target, function(s) paste0(s, collapse = '')))
  # seqs.target.msa <- DNAStringSet(seqs.target.msa)
  # alignment <- msa(seqs.target.msa)
  
  # Run the alignment
  
  return(seqs.target)
  
}

#' @export
reduceAlnFreq <- function(aln, 
                          min.presence.freq = 3,
                          show.echo = F){
  
  aln.str.flag = F
  
  # Input handling
  if(is.vector(aln) && is.character(aln)){
    aln <- aln2mx(aln)
    aln.str.flag = T
  }
  
  # If input is matrix, ensure correct type
  if(!is.matrix(aln)){
    stop("Input must be either a character vector of aligned sequences or a matrix.")
  }
  
  if(min.presence.freq > nrow(aln)) stop('Number of aligned sequences is less than min.gap.freq')
  
  nongap.freq = colSums(aln != '-')
  
  aln = aln[, nongap.freq >= min.presence.freq, drop=F]
  
  idx.removed = which(rowSums(aln != 0) == 0)
  if(length(idx.removed) > 0){
    if(show.echo) pokazAttention('The amount of removes sequences is', length(idx.removed))
    aln = aln[-idx.removed, , drop=F]
  }
  
  if(aln.str.flag){
    aln = mx2aln(aln)
  }
  
  return(aln)
  
}

#' @export
reduceAlnFlank <- function(aln, 
                           freq.min.presence = 3,
                           prop.same = 0.9,
                           stretch.length = 5,
                           show.echo = F){
  
  aln.str.flag = F
  
  # Input handling
  if(is.vector(aln) && is.character(aln)){
    aln <- aln2mx(aln)
    aln.str.flag = T
  }
  
  # If input is matrix, ensure correct type
  if(!is.matrix(aln)){
    stop("Input must be either a character vector of aligned sequences or a matrix.")
  }
  
  if(freq.min.presence > nrow(aln)) stop('Number of aligned sequences is less than freq.min.presence')
  
  
  profile = mx2profile(aln)
  freq.max = colMax(profile)
  freq.cnt = colSums(profile)
  
  idx.fit = which( (freq.max >= freq.min.presence) & ((freq.max / freq.cnt) >= prop.same) )
  idx.fit.diff = diff(idx.fit)
  idx.fit.diff[idx.fit.diff != 1] = 0
  idx.stretch = findOnes(idx.fit.diff)
  idx.stretch$len = idx.stretch$end - idx.stretch$beg + 1
  idx.stretch = idx.stretch[idx.stretch$len >= stretch.length,]
  
  idx.start = idx.fit[idx.stretch$beg[1]]
  idx.end = idx.fit[idx.stretch$end[nrow(idx.stretch)] + 1]
  
  
  # msaplot(aln[,1:(idx.start+4)])
  # msaplot(aln[,(idx.end - 4):ncol(aln)])
  
  aln = aln[, idx.start:idx.end, drop = F]
  
  idx.removed = which(rowSums(aln != 0) == 0)
  if(length(idx.removed) > 0){
    if(show.echo) pokazAttention('The amount of removes sequences is', length(idx.removed))
    aln = aln[-idx.removed, , drop=F]
  }
  
  if(aln.str.flag){
    aln = mx2aln(aln)
  }
  
  return(aln)
  
}

#' @export
getFmilyEdges <- function(edges, components, i.comp){
  sv.comp = names(components)[components == i.comp]
  if(length(sv.comp) == 0){
    pokazAttention('Component', i.comp, 'does not exist')
    return(NULL)
  }
  
  edges.comp = edges[(edges[,1] %in% sv.comp) & (edges[,2] %in% sv.comp),,drop=F]
  
  return(edges.comp)
  
}

scoreFlankCoverage <- function(s1, s2, n.cut, wsize = 15, nmatch = 12) {
  score.beg.1 <- max(dotcover(s1, s2[1:n.cut], 15, 12, strand = 1), 
                     dotcover(s1, s2[1:n.cut], 15, 12, strand = -1))
  score.end.1 <- max(dotcover(s1, s2[length(s2) - n.cut + (1:n.cut)], 15, 12, strand = 1),
                     dotcover(s1, s2[length(s2) - n.cut + (1:n.cut)], 15, 12, strand = -1))
  
  score.beg.2 <- max(dotcover(s2[1:n.cut], s1, 15, 12, strand = 1), 
                     dotcover(s2[1:n.cut], s1, 15, 12, strand = -1))
  score.end.2 <- max(dotcover(s2[length(s2) - n.cut + (1:n.cut)], s1, 15, 12, strand = 1),
                     dotcover(s2[length(s2) - n.cut + (1:n.cut)], s1, 15, 12, strand = -1))
  
  score.beg = max(score.beg.1, score.beg.2)  # min or max?
  score.end = max(score.end.1, score.end.2)  # min or max?
  
  score.tot <- max(score.beg, score.end)
  
  return(score.tot)
}

