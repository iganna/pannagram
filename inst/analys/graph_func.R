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
#' 
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
  return(g.comp)
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
  
  nodes.mutual =     data.frame(node = paste0('N', graph.mutual.comp$membership), 
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
    edges.compact = refineDirectEdges(edges.compact, echo = F)
    cat(' done!\n')
  }
  
  res.list = list(edges = edges.compact, nodes = nodes, nodes.traits = nodes.traits)
  if(return.nest){
    res.list[['nestedness']] = res.nest
  }
  return(res.list)
  
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
  E(g)$name <- paste(edges.compact[,1], edges.compact[,2], sep = '-')
  
  # Find nodes, which should be kept with all their edges and delete them from the graph
  nodes.keep = c()
  n.nodes.keep = -1
  while(length(nodes.keep) != n.nodes.keep){
    n.nodes.keep = length(nodes.keep)
    # print(n.nodes.keep)
    print(vcount(g))
    
    deg.in <- igraph::degree(g, mode = "in")
    deg.out <- igraph::degree(g, mode = "out")
    # tails = names(deg.in)[((deg.in + deg.out) == 1) & ((deg.in * deg.out) == 0) | ((deg.in * deg.out) == 1)]
    tails = names(deg.in)[((deg.in + deg.out) == 1) & ((deg.in * deg.out) == 0) ]
    
    # if( 'R4396' %in% tails) stop()
    
    nodes.keep = c(nodes.keep, tails)
    g <- delete_vertices(g, tails)
  }
  
  # Remove connected components with size 1 and keep node names
  g.comp <- igraph::components(g)
  id.single = which(g.comp$csize == 1)
  nodes.single = names(g.comp$membership)[g.comp$membership %in% id.single]
  nodes.keep = c(nodes.keep, nodes.single)
  g <- delete_vertices(g, nodes.single)
  
  # Refine each connected component
  g.comp <- igraph::components(g)
  edges.polised = c()
  for(i.comp in 1:g.comp$no){
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
#' @param sim.cutoff similarity to establish an edge between sequences
#' this is for both: (1) sequence blast similarity, (2) length similarity.
#' 
#' @param i.len.field index for the field, where you see the length of the sequence in names
#'
#' @return A matrix representing the edges of the graph
#' 
old_getGraphFromBlast <- function(bl.res, sim.cutoff = 0.85, i.len.field = 5){
  bl.res = bl.res[bl.res$V1 != bl.res$V8,]
  
  bl.res.init = bl.res
  bl.res = bl.res[bl.res$V6 >= sim.cutoff * 100,]
  
  res.nest = findNestedness(bl.res, use.strand = F)
  
  message('be careful with extracting sequence lengths, it a positionsl information in their names')
  res.nest.len = sapply(unique(c(res.nest$V1, res.nest$V8)), function(s) as.numeric(strsplit(s, '\\|')[[1]][i.len.field]))
  
  res.nest$len1 = res.nest.len[res.nest$V1]
  res.nest$len8 = res.nest.len[res.nest$V8]
  res.nest$p1 = res.nest$C1 / res.nest$len1
  res.nest$p8 = res.nest$C8 / res.nest$len8
  
  res.nest.sim = res.nest[(res.nest$p1 >= sim.cutoff) | 
                            (res.nest$p8 >= sim.cutoff),]
  
  ## Creating the graph
  
  # all edges
  idx = res.nest$p1 >= sim.cutoff
  edges = cbind(res.nest$V1[idx], res.nest$V8[idx])
  idx = res.nest$p8 >= sim.cutoff
  edges = rbind(edges, cbind(res.nest$V8[idx], res.nest$V1[idx]))
  te.enges.names = unique(c(edges[,1], edges[,2]))
  
  # nodes
  idx = (res.nest$p1 >= sim.cutoff) & (res.nest$p8 >= sim.cutoff)
  te.nodes = cbind(res.nest$V1[idx], res.nest$V8[idx])
  
  te.rest = setdiff(te.enges.names, c(te.nodes[,1], te.nodes[,2]))
  
  te.nodes.graph <- igraph::make_graph(t(te.nodes), directed = T)
  te.nodes.graph <- igraph::simplify(te.nodes.graph)
  te.nodes.comp <- igraph::components(te.nodes.graph)
  
  nodes = data.frame(node = paste0('N', te.nodes.comp$membership), 
                     te = names(te.nodes.comp$membership))
  nodes.rest = data.frame(node = paste('R', (1:length(te.rest)), sep = ''), te = te.rest)
  nodes = rbind(nodes, nodes.rest)
  
  rownames(nodes) = nodes$te
  
  
  nodes.cnt = data.frame(cnt = c(table(nodes$node)))
  nodes.cnt$node = rownames(nodes.cnt)
  
  
  # Redefine edges but with node names
  idx.endes = (edges[,1] %in% nodes$te) & (edges[,2] %in% nodes$te)
  b.graph = cbind(nodes[edges[idx.endes,1], 'node'],nodes[edges[idx.endes,2], 'node'])
  b.graph = unique(b.graph)
  # b.graph = b.graph[b.graph[,1] != b.graph[,2],]
  b.graph.uni = b.graph[b.graph[,1] == b.graph[,2],]
  b.graph = b.graph[b.graph[,1] != b.graph[,2],]
  
  length(unique(c(b.graph[,1], b.graph[,2])))
  
  b.graph = refineDirectEdges_old(b.graph)
  
  return(b.graph)
}

#' Remove direct edges from a graph. OLD VERSION!!
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
old_refineDirectEdges <- function(b.graph, echo = T){
  # reduce indirect arrows
  idx.remove = c()
  for(i.edge in 1:nrow(b.graph)){
    
    if(echo){
      if(i.edge %% 1000 == 0) print(i.edge)
    }
    
    tmp.to = b.graph[b.graph[,1] == b.graph[i.edge,1],2]
    tmp.from = b.graph[b.graph[,2] == b.graph[i.edge,2],1]
    if(length(intersect(tmp.to, tmp.from)) > 0) idx.remove = c(idx.remove, i.edge)
  }
  idx.remove = unique(idx.remove)
  b.graph = b.graph[-idx.remove,]
  
  return(b.graph)
}