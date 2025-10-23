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



#' Create a graph from BLAST results
#'
#' This function takes BLAST results in tabular format and generates a graph
#'
#' @param bl.res The BLAST results in tabular format, basically a table with columns:
#' V1 and V8 contain names
#' V2-V3 - positions of match in sequence "query"
#' V4-V5 - positions of match in sequence "subject"
#' V6 - similarity from 0 to 100
#' V7 - length
#' 
#' @param res.nest A dataframe resulting from the nestedness detection procedure.
#' This dataframe must contain specific columns, either \code{c("V1", "V8", "p1", 'p8')} 
#' or \code{c("C1", "C8", "V1", "V8", "len1", "len8")}. In case it contains the latter set 
#' of columns and not the former, two new columns \code{p1} and \code{p8} will be computed.
#' 
#' If both \code{bl.res} and 
#' \code{res.nest} are either NULL or not NULL, an error will be thrown.
#' 
#' @param sim.cutoff similarity to establish an edge between sequences
#' this is for both: (1) sequence blast similarity, (2) length similarity.
#' 
#' @param i.len.field index for the field, where you see the length of the sequence in names
#' @param collapse a boolean value indicating whether to collapse sequences with mutual similarity into nodes (default is TRUE).
#' @param refine An indicating whether to refine the graph (default is TRUE).
#' Example list of edges:
#' (1) -> (2)
#' (2) -> (3)
#' (1) -> (3)
#' 
#' The edge (1) -> (3) should de removed, because path (1) -> (2) -> (3) exists.
#'
#' @param min.length The minimum sequence length to be included in the graph (default is 0).
#' @param max.length The maximum sequence length to be included in the graph (default is Inf).
#'
#' @param return.nest An indicating whether to return the initial nestedness to create the graph (default is FALSE).
#' 
#' @param echo An indicating whether to print progress messages steps (default is TRUE).
#'
#' @return The graph represented as the list with the following components:
#'   \itemize{
#'     \item{edges}{A matrix representing the edges of the graph}
#'     \item{nodes}{A matrix reflecting the relationship between BLAST queries and collapsed nodes.}
#'     \item{node.traits}{A data frame representing traits associated with nodes.}
#'     \term{nestedness}{A dataframe resulting from the nestedness detection procedure.}
#'   }
#' @export
getGraphFromBlast <- function(bl.res = NULL, 
                              res.nest = NULL,
                              sim.cutoff = 0.85, 
                              i.len.field = 5, 
                              collapse = T, 
                              refine = T,
                              min.length = 0,
                              max.length = Inf,
                              return.nest = F,
                              echo = T){
  
  if (is.null(bl.res) == is.null(res.nest)) {
    stop("Both parameters cannot be NULL or both not NULL.")
  }
  
  if(is.null(res.nest)){
    
    # nestedness
    bl.res = bl.res[bl.res$V1 != bl.res$V8,]
    bl.res = bl.res[bl.res$V6 >= sim.cutoff * 100,]
    if(nrow(bl.res) == 0) return(NULL)
    
    # Length of sequences (optimised approach)
    res.nest.len = sapply(unique(c(bl.res$V1, bl.res$V8)), function(s) as.numeric(strsplit(s, '\\|')[[1]][i.len.field]))
    if(sum(is.na(res.nest.len)) != 0){
      stop(paste('Be careful with extracting sequence lengths.\n',
                 'It a positionsl information in their names.\n',
                 'Current length slot is ', i.len.field, sep = ''))
    }
    # lengths are different from the default
    arg.defaults <- formals(sys.function(sys.nframe()))
    if((min.length != arg.defaults$min.length) | (max.length != arg.defaults$max.length)){
      message(paste('New lengths', min.length, max.length))
      res.nest.len = res.nest.len[res.nest.len >= min.length]
      res.nest.len = res.nest.len[res.nest.len <= max.length]
      bl.res = bl.res[bl.res$V1 %in% names(res.nest.len),]
      bl.res = bl.res[bl.res$V8 %in% names(res.nest.len),]
    }
    
    # Find nestedness
    cat('Find Nestedness...')
    res.nest = findNestedness(bl.res, use.strand = F)
    cat(' done!\n')
    
  
    # Gen lengths
    res.nest$len1 = res.nest.len[res.nest$V1]
    res.nest$len8 = res.nest.len[res.nest$V8]
    res.nest$p1 = res.nest$C1 / res.nest$len1
    res.nest$p8 = res.nest$C8 / res.nest$len8
  
  } else {
    
    required.cols <- c("V1", "V8", "p1", 'p8')
    if (!all(required.cols %in% names(res.nest))) {
      required.cols <- c("C1", "C8", "V1", "V8", "len1", "len8")
      if (!all(required.cols %in% names(res.nest))) {
        stop(paste0(c("res.nest must have the following columns:", required.cols), collapse = ', '))
      }
      res.nest$p1 = res.nest$C1 / res.nest$len1
      res.nest$p8 = res.nest$C8 / res.nest$len8
    }
  }
  
  
  # Remain only those blast hits, which satisfy sim.cutoff
  res.nest.sim = res.nest[(res.nest$p1 >= sim.cutoff) | 
                            (res.nest$p8 >= sim.cutoff),]
  
  # get edges of nestedness: (1) -> (2), i.e. (1) is covered by (2)
  idx.1.to.2 = res.nest$p1 >= sim.cutoff
  edges = cbind(res.nest$V1[idx.1.to.2], res.nest$V8[idx.1.to.2])
  idx.2.to.1 = res.nest$p8 >= sim.cutoff
  edges = rbind(edges, cbind(res.nest$V8[idx.2.to.1], res.nest$V1[idx.2.to.1]))
  
  if(!collapse){
    res.list = list(edges = edges, nodes = c(), nodes.traits = c())
    if(return.nest){
      res.list[['nestedness']] = res.nest
    }
    return(res.list)
  }
  
  # define nodes
  # if some have mutual arrows - they should go to one node
  # Idea: create the graph on those elements, which have mutual arrows and take the connected components
  idx.mutual = idx.1.to.2 & idx.2.to.1
  edges.mutual = cbind(res.nest$V1[idx.mutual], res.nest$V8[idx.mutual])
  names.rest = setdiff(c(edges), c(edges.mutual))
  
  graph.mutual <- igraph::simplify(igraph::make_graph(t(edges.mutual), directed = T))
  graph.mutual.comp <- igraph::components(graph.mutual)
  
  nodes.mutual = data.frame(node = paste0('N', graph.mutual.comp$membership), 
                            name = names(graph.mutual.comp$membership))
  if(length(names.rest) > 0){
    nodes.rest = data.frame(node = paste('R', (1:length(names.rest)), sep = ''), 
                            name = names.rest)  
  } else {
    nodes.rest = c()
  }
  
  nodes = rbind(nodes.mutual, nodes.rest)
  rownames(nodes) = nodes$name
  
  nodes.traits = data.frame(cnt = c(table(nodes$node)))
  nodes.traits$node = rownames(nodes.traits)
  
  # Redefine edges but with node names
  edges.compact = cbind(nodes[edges[,1], 'node'], nodes[edges[,2], 'node'])
  edges.compact = edges.compact[edges.compact[,1] != edges.compact[,2],]
  edges.compact = unique(edges.compact)
  
  
  nodes.traits$in.graph = nodes.traits$node %in% c(edges.compact)
  
  # Some sequences are the same, but do not form any unmutual nestedness,so no edges in the graph
  table(nodes.traits$cnt[!nodes.traits$in.graph])
  
  if(refine){
    cat('Refine...')
    if(nrow(edges.compact) > 1){
      edges.compact = refineDirectEdges(edges.compact, echo = F)  
    }
    cat(' done!\n')
  }
  
  res.list = list(edges = edges.compact, nodes = nodes, nodes.traits = nodes.traits)
  if(return.nest){
    res.list[['nestedness']] = res.nest
  }
  return(res.list)
  
}

#' Generate Graph Edges from Nestedness Results
#'
#' This function extracts directed edges between query and target names
#' based on a similarity cutoff from nestedness analysis results.
#'
#' @param res.nest A list or data frame containing \code{name.query},
#'   \code{name.target}, \code{coverage.query}, and \code{coverage.target}.
#' @param coverage.cutoff Numeric value specifying the similarity cutoff.
#'   Values between 0–1 or 10–100 (interpreted as %) are accepted.
#'
#' @return A two-column matrix of unique edges representing connections
#'   that meet the similarity threshold.
#' @export
getGraphFromNestedness <- function(res.nest,
                                   coverage.cutoff){
  
  if (coverage.cutoff >= 10 && coverage.cutoff <= 100) {
    coverage.cutoff <- coverage.cutoff / 100
  } else if (coverage.cutoff < 0 || coverage.cutoff > 1) {
    stop("coverage.cutoff does not make sense")
  }
  
  idx.1.to.2 = res.nest$coverage.query >= coverage.cutoff
  edges = cbind(res.nest$name.query[idx.1.to.2], res.nest$name.target[idx.1.to.2])
  idx.2.to.1 = res.nest$coverage.target >= coverage.cutoff
  edges = rbind(edges, cbind(res.nest$name.target[idx.2.to.1], res.nest$name.query[idx.2.to.1]))
  
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
  names.mutual = names(compsStrong$membership)[compsStrong$membership != 1]
  names.rest =  names(compsStrong$membership)[compsStrong$membership == 1]
  
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
#' getBypassedEdgeIdx(edges)
#' @export
getBypassedEdgeIdx <- function(edges.compact) {
  # Get unique node names
  node.names <- unique(c(edges.compact))
  node.n <- length(node.names)
  
  # Map nodes to indices
  id1 <- match(edges.compact[,1], node.names)
  id2 <- match(edges.compact[,2], node.names)
  
  # Build sparse adjacency matrix
  A <- sparseMatrix(i = id1, j = id2, x = TRUE, dims = c(node.n, node.n))
  
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
#' Filters a coverage matrix (`mx.cov`) based on coverage percentage, sequence 
#' length, and minimum copy number.
#'
#' @param mx.cov Data frame with coverage and length information. Must contain 
#'   columns: `name.query`, `name.target`, 
#'            `length.query`, `length.target`,
#'            `coverage.query`, `coverage.target` .
#' @param cov.cutoff Integer. Minimum coverage percentage to retain entries (default = NULL).
#' @param min.len Integer. Minimum sequence length to retain entries (default = NULL).
#' @param max.len Integer. Maximum sequence length to retain entries (default = NULL).
#' @param min.copy Integer. Minimum number of copies to retain entries (default = NULL).
#'
#' @return A filtered data frame with the same structure as `mx.cov`.
#' @examples
#' filtered <- filterCoverageMatrix(mx.cov, cov.cutoff = 80, min.len = 500)
#' @export
filterCoverageMatrix <- function(mx.cov,
                                 cov.cutoff = NULL,
                                 min.len    = NULL,
                                 max.len    = NULL,
                                 min.copy   = NULL,
                                 names.remove = NULL,
                                 echo = F) {
  
  if(echo){
    pokaz('Initial number of rows:', nrow(mx.cov))
  }
  
  # ---- Filter by names to remove ----
  if(!is.null(names.remove)){
    mx.cov = res.sim[!(mx.cov$name.query %in% names.remove) & 
                        !(mx.cov$name.target %in% names.remove),]
    
    if(echo){
      pokaz('Number of rows after name filtration:', nrow(mx.cov))
    }
  }
  
  # ---- Filter by coverage cutoff (percentage) ----
  if (!is.null(cov.cutoff)) {
    thr <- cov.cutoff / 100
    mx.cov2 <- mx.cov[(mx.cov$coverage.query  >= thr) |
                       (mx.cov$coverage.target >= thr), ]
    
    if(echo){
      pokaz('Number of rows after filtration by coverage:', nrow(mx.cov))
    }
  }
  
  # ---- Filter by minimum length ----
  if (!is.null(min.len)) {
    mx.cov <- mx.cov[(mx.cov$length.query  >= min.len) &
                       (mx.cov$length.target >= min.len), , drop = FALSE]
    
    if(echo){
      pokaz('Number of rows after filtration by minimum length:', nrow(mx.cov))
    }
  }
  
  # ---- Filter by maximum length ----
  if (!is.null(max.len)) {
    mx.cov <- mx.cov[(mx.cov$length.query  <= max.len) &
                       (mx.cov$length.target <= max.len), , drop = FALSE]
    
    if(echo){
      pokaz('Number of rows after filtration by maximum length:', nrow(mx.cov))
    }
  }
  
  # ---- Filter by minimum number of copies ----
  if (!is.null(min.copy)) {
    
    if (is.null(cov.cutoff)) {
      warning("'cov.cutoff' sould be provided to filter by copy-number")
    }
    thr <- cov.cutoff / 100
    
    mx.cov.equal <- mx.cov[(mx.cov$coverage.query  >= thr) &
                             (mx.cov$coverage.target >= thr), , drop = FALSE]
    
    if (nrow(mx.cov.equal) > 0) {
      
      if(min.copy == 2){
        sv.passed = c(mx.cov.equal$name.query, mx.cov.equal$name.target)
      } else {
        # getGraphComponents should return:
        #   $membership — vector of component IDs (names = vertices)
        #   $csize      — vector of component sizes (index = component ID)
        g.comp <- getGraphComponents(mx.cov.equal[, c("name.query", "name.target")])
        
        comp.passed <- which(g.comp$csize >= min.copy)
        sv.passed   <- names(g.comp$membership)[g.comp$membership %in% comp.passed]
      }
      
      mx.cov <- mx.cov[(mx.cov$name.query  %in% sv.passed) &
                         (mx.cov$name.target %in% sv.passed), , drop = FALSE]        
      
    } else {
      # Return an empty data frame with the same structure if no edges remain
      mx.cov <- mx.cov[0, , drop = FALSE]
    }
    if(echo){
      pokaz('Number of rows after filtration by copies:', nrow(mx.cov))
    }
  }
  
  return(mx.cov)
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
                        echo = F) {
  
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
    
    if(echo){
      pokaz('Number of edges after bridge filtration (undirected):', nrow(edges.no.bridges))
    }
    
    edges = edges[(edges[,1] %in% sv.remain) & 
                    (edges[,2] %in% sv.remain),, drop= F]
    
    if(echo){
      pokaz('Number of edges after bridge filtration:', nrow(edges))
    }
  }
  
  
  return(edges)
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
putEdgesBack <- function(edges, edges.init, echo=F, dominant.effect = 0.7){
  
  sv.back = setdiff(c(as.matrix(edges.init)), c(as.matrix(edges)))
  if(length(sv.back) == 0){
    if(echo) pokaz('No SVs to put back', length(sv.back))
    return(edges)
  } else {
    if(echo) pokaz('Number of SVs to put back', length(sv.back))  
  }
  
  
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
    sv.node.connect.part = partition.add[sv.node.connect]
    sv.node.connect.part = sv.node.connect.part[!is.na(sv.node.connect.part)]
    part.cnt = table(sv.node.connect.part)
    
    
    i.part = NA  # In case one want to remove the "else" branch in the following if
    
    if(max(part.cnt)/sum(part.cnt) > dominant.effect){
      i.part = as.numeric(names(part.cnt)[which.max(part.cnt)])
    } else {
      sv.node.connect.part.len = as.numeric(sapply(names(sv.node.connect.part), function(s) strsplit(s, '\\|')[[1]][2]))
      s.sv.len = as.numeric(sapply(s.sv, function(s) strsplit(s, '\\|')[[1]][2]))
      idx.max.len = which.max(sv.node.connect.part.len)
      if(sv.node.connect.part.len[idx.max.len] / s.sv.len > 0.7){
        i.part = sv.node.connect.part[idx.max.len]  
      }
    }
    
    if(!is.na(i.part)){
      partition.add[s.sv] = i.part
      sv.part = names(partition.add)[partition.add == i.part]
      sv.edges.add = sv.edges.connect[(sv.edges.connect[,1] %in% sv.part) & 
                                        (sv.edges.connect[,2] %in% sv.part), , drop=F]
      if (nrow(sv.edges.add) > 0) {
        sv.edges.add.list[[sv.edges.add.counter]] = sv.edges.add
        sv.edges.add.counter = sv.edges.add.counter + 1
      } else {
        stop('Something is wrong')
      }
    } 
  }
  
  if(echo) pokaz('Number of edges to put back:', sv.edges.add.counter)
  
  if (length(sv.edges.add.list) > 0) {
    edges.add = do.call(rbind, sv.edges.add.list)
  } else{
    return(edges)
  }
  
  edges = rbind(edges, edges.add)
  return(edges)
  
}




