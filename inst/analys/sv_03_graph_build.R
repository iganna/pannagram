# Get SV positions, GFF files, dencity files and consensys sequences
# Find SVs and create GFF file

suppressMessages({ 
  # library(rhdf5)
  # library('foreach')
  # library(doParallel)
  library("optparse")
  library(pannagram)
  # library(crayon)
  library(ggplot2)
  
  library(igraph)
  # library(ggnet)
  library(network)
  library(GGally)
  library(Matrix)
})

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option("--path.figures", type = "character", default = "",   help = "Path to folder with figures"),
  make_option("--path.sv", type = "character", default = NULL, help = "Path to sv dir"),
  make_option(c("--cores"),     type = "integer",   default = 1, help = "number of cores to use for parallel processing")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# Paths

path.sv <- opt$path.sv
if(!dir.exists(path.sv)) stop(paste0('No SV dir!', path.sv))

path.figures <- opt$path.figures
if(!dir.exists(path.figures)) stop(paste0('No SV figures dir', path.figures))

# ***********************************************************************
# ---- Variables ----
sim.cutoff.true = 85
sim.cutoff = 90  # have no together with the similarity cutoffs for the initial simsearch results
dominant.effect = 0.7

# Filtration of the graph
len.cutoff = 400
min.comp.size = 4
min.len = 200

min.copy = 2

seed.value = 239
plot.size = 8


# Binning
source(system.file("analys/sv_variables.R", package = "pannagram"))

# Variables for plotting
flag.plot = T
flag.save.plot = F
i.plot = 1

# Echo
show.echo = T

# ***********************************************************************
# ---- Reading the data ----
if(show.echo) pokaz('Reading the data...')

res.sim.file = file.path(path.sv, 'seq_sv_large_85_85.txt')
if(!file.exists(res.sim.file)){
  pokazAttention('All SVs are different, can not biuld a graph.')
  quit(save = "no", status = 0)
}
res.sim = read.table(res.sim.file, header = T)

seqs = readFasta(file.path(path.sv, 'seq_sv_large.fasta'))

# ***********************************************************************
# ---- Get initial edges ----
if(show.echo) pokaz('Get initial edges...')

# edges = res.sim[, c("name.query", "name.target")]
# edges = unique(edges)

res.sim = filterCoverageMatrix(res.sim,
                               min.len = min.len,
                               echo = show.echo)

edges = getGraphFromNestedness(res.sim, coverage.cutoff = sim.cutoff.true)

# if(flag.plot){
if(F){
  g <- network(edges, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
  g.names = network.vertex.names(g)

  set.seed(seed.value)
  p <- ggnet2(g, label = F, edge.color = "black",
              node.size = 1,
              color = '#468B97',
              arrow.gap = 0.01, arrow.size = 2,
  )

  name.output = sprintf("graph_%02d_init", i.plot)
  i.plot = i.plot + 1
  savePNG(p, path = path.figures, name = name.output,
          width = plot.size, height = plot.size)
  if(flag.save.plot){
    saveRDS(p, file.path(path.figures, paste0(name.output, '.rds'))) 
  }

  rmSafe(g)
  rmSafe(p)
}

# ***********************************************************************
# ---- Remain those nodes, that have at least two sequences + Length cutoff ----

if(show.echo) pokaz('Remain those node, that have at least two sequences + Length cutoff...')
res.sim.major = filterCoverageMatrix(res.sim,
                                     cov.cutoff = sim.cutoff,
                                     min.copy = min.copy, 
                                     echo = show.echo)

# edges2 = res.sim.major[, c("name.query", "name.target")]
# edges2 = unique(edges2)

edges2 = getGraphFromNestedness(res.sim.major, coverage.cutoff = 85)

print(edges2)

if(flag.plot){
  suppressMessages(suppressWarnings({
    
  g <- network(edges2, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
    
  set.seed(seed.value)
  p <- ggnet2(g, label = F, edge.color = "black",
              node.size = 1,
              color = '#468B97',
              arrow.gap = 0.01, arrow.size = 2,
  )
  
  name.output = sprintf("graph_%02d_no_doubletons", i.plot)
  i.plot = i.plot + 1
  savePNG(p, path = path.figures, name = name.output,
          width = plot.size, height = plot.size)
  if(flag.save.plot){
    saveRDS(p, file.path(path.figures, paste0(name.output, '.rds'))) 
  }
  
  rmSafe(g)
  rmSafe(p)
  
  }))
}


# ***********************************************************************
# ---- Create a compact graph ----

if(show.echo) pokaz('Create a compact graph...')

graph.compact = getGraphCompact(edges2)
edges.compact = graph.compact$edges

if(T){
# if(flag.plot){
  suppressMessages(suppressWarnings({
    
    g <- network(edges.compact, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
    g.names = network.vertex.names(g)
    
    # Add cluster color
    g %v% "size" <- graph.compact$nodes.size[g.names]
    
    set.seed(seed.value)
    p <- ggnet2(g, label = F, edge.color = "black",
                node.size = 'size',
                color = '#468B97',
                arrow.gap = 0.01, arrow.size = 2,
    ) + theme(legend.position = "none")
    p
    
    name.output = sprintf("graph_%02d_compact", i.plot)
    i.plot = i.plot + 1
    savePNG(p, path = path.figures, name = name.output,
            width = plot.size, height = plot.size)
    
    rmSafe(g)
    rmSafe(p)
    
  }))
}

# ***********************************************************************
# Components with two nodes are not needed to be filtered, therefore min.comp.size = 3

edges.compact = filterEdges(edges.compact, min.comp.size = 3)

# ***********************************************************************
# ---- Refine baypass edges ----

idx.bypassed = getBypassedEdgeIdx(edges.compact)

# Remove corresponding edges from edges2
if(length(idx.bypassed) > 0){
  idx.edges.remove = c()
  for(i.edge in idx.bypassed){
    nodes.from = graph.compact$nodes.list[[edges.compact[i.edge, 1]]]
    nodes.to = graph.compact$nodes.list[[edges.compact[i.edge, 2]]]
    i.edges.remove = which((edges2[,1] %in% nodes.from) & (edges2[,2] %in% nodes.to))
    idx.edges.remove = c(idx.edges.remove, i.edges.remove)
  }
  
  if(show.echo) pokaz('Number of bypass edges is', length(idx.edges.remove))
  if(length(idx.edges.remove) > 0){
    comp.before = getGraphComponents(edges2)
    edges2 = edges2[-idx.edges.remove,]
    comp.after = getGraphComponents(edges2)
    # Check number of components
    if(comp.before$no != comp.after$no){
      stop("Some deleted edges were crucial")
    }
  }
}

# ***********************************************************************
# Update the compact structure and visualize again

graph.compact = getGraphCompact(edges2)
edges.compact = graph.compact$edges

if(T){
  # if(flag.plot){
  suppressMessages(suppressWarnings({
    
    g <- network(edges.compact, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
    g.names = network.vertex.names(g)
    
    # Add cluster color
    g %v% "size" <- graph.compact$nodes.size[g.names]
    
    set.seed(seed.value)
    p <- ggnet2(g, label = F, edge.color = "black",
                node.size = 'size',
                color = '#468B97',
                arrow.gap = 0.01, arrow.size = 2,
    ) + theme(legend.position = "none")
    p
    
    name.output = sprintf("graph_%02d_compact_no_baypass", i.plot)
    i.plot = i.plot + 1
    savePNG(p, path = path.figures, name = name.output,
            width = plot.size, height = plot.size)
    
    rmSafe(g)
    rmSafe(p)
    
  }))
}

# ***********************************************************************
# Remove parasitic edges

cutoff.remain.edges = 0.7
flank.cover.cutoff = 0.8

# Consider all nodes that have at least two outgoing edges
stat.neighbours.all = c()
nodes.parasite = unique(edges.compact[duplicated(edges.compact[,1]),1])
for(node.from in nodes.parasite){
  stat.neighbours = data.frame(edge.id = which(edges.compact[,1] == node.from))
  
  stat.neighbours$node.from = node.from
  stat.neighbours$node.to = edges.compact[stat.neighbours$edge.id,2]
  stat.neighbours$size.from = graph.compact$nodes.size[node.from] 
  stat.neighbours$size.to = graph.compact$nodes.size[stat.neighbours$node.to] 
  
  edges.nei = edges2[edges2[,1] %in% graph.compact$nodes.list[[node.from]],]
  
  # Find how many are connected on the "form"-side
  edges.nei.mod = edges.nei
  edges.nei.mod[,2] = graph.compact$nodes[edges.nei[,2],]$node
  edges.nei.mod = unique(edges.nei.mod)
  
  connect.from = split(edges.nei.mod[,1], edges.nei.mod[,2])
  connect.from.len = unlist(lapply(connect.from, length))
  stat.neighbours$n.from = connect.from.len[stat.neighbours$node.from]
  
  # Find how many are connected on the "to"-side
  names.to = unique(edges.nei[,2])
  nodes.to = graph.compact$nodes[names.to,]$node
  
  cnt.nodes.to = table(nodes.to)
  stat.neighbours$n.to = cnt.nodes.to[stat.neighbours$node.to]
  if(sum(is.na(stat.neighbours)) > 0) stop('Wrong names of neighbours.')
  
  # Compute proportions
  stat.neighbours$p.from = stat.neighbours$n.from / stat.neighbours$size.from
  stat.neighbours$p.to = stat.neighbours$n.to / stat.neighbours$size.to
  
  stat.neighbours$remain = (stat.neighbours$p.to >= cutoff.remain.edges) & 
                                      (stat.neighbours$p.from >= cutoff.remain.edges)
  
  # Check whether the sequences are intersect by flanking regions
  
  edges.nei.mod = edges.nei
  edges.nei.mod[,2] = graph.compact$nodes[edges.nei[,2],]$node
  for(irow in which(stat.neighbours$remain)){
    node.to = stat.neighbours$node.to[irow]
    # stop()
    
    idx.tmp = which(edges.nei.mod[,2] == node.to)[1]
    edge.tmp = edges.nei[idx.tmp,]
    
    s1 = seq2nt(seqs[edge.tmp[1]])
    s2 = seq2nt(seqs[edge.tmp[2]])
    
    n.cut = min(round(length(s1) / flank.cover.cutoff), length(s2))
    
    score.beg = dotcover(s1, s2[1:n.cut], 15, 12)
    score.end = dotcover(s1, s2[length(s2) - n.cut + (1:n.cut)], 15, 12)
    score.tot = max(score.beg, score.end)
    
    if(score.tot < flank.cover.cutoff){
      stat.neighbours$remain[irow] = F
    }
  }
  
  stat.neighbours.all = rbind(stat.neighbours.all, stat.neighbours)
  # stop()
}

if(!is.null(stat.neighbours.all)){
  idx.edge.remove = stat.neighbours.all$edge.id[!stat.neighbours.all$remain]  
} else {
  idx.edge.remove = c()
}

# Remove corresponding edges from edges2
if(length(idx.edge.remove) > 0){
  idx.edges.remove = c()
  for(i.edge in idx.edge.remove){
    nodes.from = graph.compact$nodes.list[[edges.compact[i.edge, 1]]]
    nodes.to = graph.compact$nodes.list[[edges.compact[i.edge, 2]]]
    i.edges.remove = which((edges2[,1] %in% nodes.from) & (edges2[,2] %in% nodes.to))
    idx.edges.remove = c(idx.edges.remove, i.edges.remove)
  }
  
  if(show.echo) pokaz('Number of edges to remove', length(idx.edges.remove))
  if(length(idx.edges.remove) > 0){
    comp.before = getGraphComponents(edges2)
    edges2 = edges2[-idx.edges.remove,]
    comp.after = getGraphComponents(edges2)
    # Check number of components
  }
}


# ***********************************************************************
# ---- Update the compact structure and visualize again ----

graph.compact = getGraphCompact(edges2)
edges.compact = graph.compact$edges

if(T){
  # if(flag.plot){
  suppressMessages(suppressWarnings({
    
    g <- network(edges.compact, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
    g.names = network.vertex.names(g)
    
    # Add cluster color
    g %v% "size" <- graph.compact$nodes.size[g.names]
    
    set.seed(seed.value)
    p <- ggnet2(g, label = T, edge.color = "black",
                label.size = 2,
                node.size = 'size',
                color = '#468B97',
                arrow.gap = 0.01, arrow.size = 2,
    ) + theme(legend.position = "none")
    p
    
    name.output = sprintf("graph_%02d_compact_no_repeats", i.plot)
    i.plot = i.plot + 1
    savePNG(p, path = path.figures, name = name.output,
            width = plot.size, height = plot.size)
    
    rmSafe(g)
    rmSafe(p)
    
  }))
}


# ***********************************************************************
# ---- Remove umbrella edges ----

coverage.umbrella.children = 100

stat.neighbours.all = c()
nodes.umbrella = setdiff(unique(edges.compact[duplicated(edges.compact[,2]),2]),
                         edges.compact[,1])
for(node.to in nodes.umbrella){
  stat.neighbours = data.frame(edge.id = which(edges.compact[,2] == node.to))
  
  stat.neighbours$node.from = edges.compact[stat.neighbours$edge.id,1]
  stat.neighbours$node.to = node.to
  stat.neighbours$size.from = graph.compact$nodes.size[stat.neighbours$node.from] 
  stat.neighbours$size.to = graph.compact$nodes.size[stat.neighbours$node.to] 
  
  edges.nei = edges2[edges2[,2] %in% graph.compact$nodes.list[[node.to]], ]
  
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
    
    score.beg = dotcover(s1, s2[1:n.cut], 15, 12)
    score.end = dotcover(s1, s2[length(s2) - n.cut + (1:n.cut)], 15, 12)
    score.tot = max(score.beg, score.end)
    
    if(score.tot < flank.cover.cutoff){
      stat.neighbours$remain[irow] = F
    }
  }
  
  # Remain the longest and those which match with the longest
  if(sum(stat.neighbours$remain) > 1){
    stat.neighbours$len.from[!stat.neighbours$remain] = 0
    irow.longest = which.max(stat.neighbours$len.from)
    
    # TODO
    # for(irow in which(stat.neighbours$remain)){
    #   if(irow == irow.longest) next
    #   
    #   s1 = 
    #   s2 = 
    #   score.tot = dotcover(s1, s2, 15, 12)
    #   if(score.tot < coverage.umbrella.children){
    #     stat.neighbours$remain[irow] = F
    #   }
    # }
    
  }
  


  # Keep the results in the common dataframe
  stat.neighbours.all = rbind(stat.neighbours.all, stat.neighbours)
  # stop()
}

if(!is.null(stat.neighbours.all)){
  idx.edge.remove = stat.neighbours.all$edge.id[!stat.neighbours.all$remain]
} else {
  idx.edge.remove = c()
}

# Remove corresponding edges from edges2
if(length(idx.edge.remove) > 0){
  idx.edges.remove = c()
  for(i.edge in idx.edge.remove){
    nodes.from = graph.compact$nodes.list[[edges.compact[i.edge, 1]]]
    nodes.to = graph.compact$nodes.list[[edges.compact[i.edge, 2]]]
    i.edges.remove = which((edges2[,1] %in% nodes.from) & (edges2[,2] %in% nodes.to))
    idx.edges.remove = c(idx.edges.remove, i.edges.remove)
  }
  
  if(show.echo) pokaz('Number of edges to remove', length(idx.edges.remove))
  if(length(idx.edges.remove) > 0){
    comp.before = getGraphComponents(edges2)
    edges2 = edges2[-idx.edges.remove,]
    comp.after = getGraphComponents(edges2)
  }
}

# ***********************************************************************
# ---- Update the compact structure and visualize again ----

graph.compact = getGraphCompact(edges2)
edges.compact = graph.compact$edges

if(T){
  # if(flag.plot){
  suppressMessages(suppressWarnings({
    
    g <- network(edges.compact, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
    g.names = network.vertex.names(g)
    
    # Add cluster color
    g %v% "size" <- graph.compact$nodes.size[g.names]
    
    set.seed(seed.value)
    p <- ggnet2(g, label = F, edge.color = "black",
                label.size = 2,
                node.size = 'size',
                color = '#468B97',
                arrow.gap = 0.01, arrow.size = 2,
    ) + theme(legend.position = "none")
    p
    
    name.output = sprintf("graph_%02d_compact_no_umbrella", i.plot)
    i.plot = i.plot + 1
    savePNG(p, path = path.figures, name = name.output,
            width = plot.size, height = plot.size)
    
    rmSafe(g)
    rmSafe(p)
    
  }))
}


# ***********************************************************************
# ---- Remove small clusters ----
if(show.echo) pokaz('Remove small clusters...')

edges22 = filterEdges(edges2, min.comp.size = min.comp.size)

# ***********************************************************************
# ---- Partitioning ----

if(show.echo) pokaz('Partitioning....')

partition = getGraphCommunities(edges22) # Louvain partition

if(flag.plot){
  suppressMessages(suppressWarnings({

  g <- network(edges22, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
  g.names = network.vertex.names(g)

  # Add cluster color
  g %v% "partition" <- as.character(partition[g.names])

  set.seed(seed.value)
  p <- ggnet2(
    g,
    label = FALSE,
    edge.color = "grey70",
    node.size = 1,
    color = "partition",
    arrow.gap = 0.01,
    arrow.size = 2
  ) + theme(legend.position = "none") +
    scale_color_viridis_d(name = "partition")

  name.output = sprintf("graph_%02d_partition", i.plot)
  i.plot = i.plot + 1
  savePNG(p, path = path.figures, name = name.output,
          width = plot.size, height = plot.size)
  if(flag.save.plot){
    saveRDS(p, file.path(path.figures, paste0(name.output, '.rds'))) 
  }

  rmSafe(g)
  rmSafe(g.names)
  rmSafe(p)
  }))
}


# ***********************************************************************
# ---- Remove edges, connecting different communities ----

if(show.echo) pokaz('Remove edges, connecting different communities....')

n.edges = 0
while(n.edges != nrow(edges22)){
  n.edges = nrow(edges22)
  edges22 = filterEdges(edges2, min.comp.size = min.comp.size, echo = show.echo)
  edges22 = filterEdges(edges22, remove.intercommunity = T, echo = show.echo)
}

partition = getGraphCommunities(edges22)

if(flag.plot){
  suppressMessages(suppressWarnings({

  g <- network(edges22, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
  g.names = network.vertex.names(g)

  # Add cluster color
  g %v% "partition" <- as.character(partition[g.names])

  set.seed(seed.value)
  p <- ggnet2(
    g,
    label = FALSE,
    edge.color = "grey70",
    node.size = 1,
    color = "partition",
    arrow.gap = 0.01,
    arrow.size = 2
  ) + theme(legend.position = "none") + scale_color_viridis_d(name = "partition")

  name.output = sprintf("graph_%02d_good_edges", i.plot)
  i.plot = i.plot + 1
  savePNG(p, path = path.figures, name = paste0(name.output),
          width = plot.size, height = plot.size)
  if(flag.save.plot){
    saveRDS(p, file.path(path.figures, paste0(name.output, '.rds'))) 
  }

  rmSafe(g)
  rmSafe(g.names)
  rmSafe(p)

  }))
}


# ***********************************************************************
# ---- Remove bridges ----

if(show.echo) pokaz('Remove bridge nodes....')

edges22.no.bridges = filterEdges(edges22, remove.bridges = T, echo = show.echo)

partition = getGraphCommunities(edges22.no.bridges)

if(flag.plot){
  suppressMessages(suppressWarnings({
      
  g <- network(edges22.no.bridges, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
  g.names = network.vertex.names(g)
  
  # Add cluster color
  g %v% "partition" <- as.character(partition[g.names])
  
  set.seed(seed.value)
  p <- ggnet2(
    g,
    label = FALSE,
    edge.color = "grey70",
    node.size = 1,
    color = "partition", 
    arrow.gap = 0.01,
    arrow.size = 2
  ) + theme(legend.position = "none") + scale_color_viridis_d(name = "partition")
  
  name.output = sprintf("graph_%02d_bridge_removal", i.plot)
  i.plot = i.plot + 1
  savePNG(p, path = path.figures, name = name.output,
          width = plot.size, height = plot.size)
  if(flag.save.plot){
    saveRDS(p, file.path(path.figures, paste0(name.output, '.rds'))) 
  }

  rmSafe(g)
  rmSafe(g.names)
  rmSafe(p)
  }))
}


# ***********************************************************************
# ---- Put back all of the SVs ----
if(show.echo)  pokaz('Put all SVs back with dominant effect...')

edges.solved = edges22.no.bridges
n.edges.solved = 0
while(n.edges.solved != nrow(edges.solved)){
  n.edges.solved = nrow(edges.solved)
  edges.solved = putEdgesBack(edges = edges.solved, 
                              edges.init = edges, 
                              echo = show.echo,
                              dominant.effect = dominant.effect)
}


edges.solved = filterEdges(edges.solved, remove.intercommunity = T)
components.solved = getGraphComponents(edges.solved)
partition.solved = getGraphCommunities(edges.solved)

# sv.edges = unique(unlist(edges))
# sv.edges.solved = unique(unlist(edges.solved))
# sv.diff = setdiff(sv.edges, sv.edges.solved)
# 
# pp.max = c()
# for(sv.tmp in sv.diff){
#   sv.tmp.connect = unique(unlist(edges[(edges[,1] == sv.tmp) | (edges[,2] == sv.tmp),]))
#   pp = patrition.solved[sv.tmp.connect]
#   # pp = pp[!is.na(pp)]
#   pp.cnt = table(pp) / length(pp)
#   pp.max = c(pp.max, max(pp.cnt))
#   if(max(pp.cnt) > dominant.effect) stop()
# }


# ***********************************************************************
# ---- Save ----

# Clustering safe
saveRDS(partition.solved, paste0(path.sv, 'sv_partition_solved.rds'))
write.table(as.matrix(partition.solved), paste0(path.sv, 'sv_partition_solved.txt'), 
            row.names = T, col.names = F, sep = '\t', quote = F)

## Save edges
saveRDS(edges, paste0(path.sv, 'edges_solved.rds'))
write.table(edges, paste0(path.sv, 'edges_solved.txt'), quote = F, sep = '\t', row.names = F, col.names = F)

if(!flag.plot){
  quit(save = "no", status = 0)
}

# ***********************************************************************
# ---- Plot ----

if(flag.plot){
  suppressMessages(suppressWarnings({
    
    g <- network(edges.solved, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
    g.names = network.vertex.names(g)
    
    g %v% "partition" <- as.character(partition.solved[g.names])
    
    set.seed(seed.value)
    p <- ggnet2(
      g,
      label = FALSE,
      edge.color = "grey70",
      node.size = 1,
      color = "partition",
    ) + theme(legend.position = "none") + scale_color_viridis_d(name = "Partition")
    
    
    name.output = sprintf("graph_%02d_solved", i.plot)
    i.plot = i.plot + 1
    savePNG(p, path = path.figures, name = name.output,
            width = plot.size, height = plot.size)
    if(flag.save.plot){
      saveRDS(p, file.path(path.figures, paste0(name.output, '.rds'))) 
    }
    
    rmSafe(p)
  }))
}

# ***********************************************************************
# ---- Colored by length ----


if(flag.plot){
  suppressMessages(suppressWarnings({
    
    g.names.len = as.numeric(sapply(g.names, function(s) strsplit(s, '\\|')[[1]][2]))
    g.names.len.bin <- as.character(cut(g.names.len, breaks = len.bins, right = FALSE, labels = len.labels))
    g %v% "colors" = as.character(g.names.len.bin)
    
    set.seed(seed.value)
    p <- ggnet2(
      g,
      label = FALSE,
      edge.color = "grey70",
      node.size = 1,
      color = "colors",
      palette = color.len,
    ) +  guides(size = F) + theme(aspect.ratio = 1) + 
      guides(color = guide_legend(title = "Length, bp"))
    
    name.output = sprintf("graph_%02d_solved_colored", i.plot)
    i.plot = i.plot + 1
    savePNG(p, path = path.figures, name = name.output,
            width = plot.size + 1, height = plot.size)
    if(flag.save.plot){
      saveRDS(p, file.path(path.figures, paste0(name.output, '.rds'))) 
    }
    
    rmSafe(p)
  }))
}

# ***********************************************************************
# ---- With labels ----
 
if(flag.plot){
  suppressMessages(suppressWarnings({
    
    g.label = rep('', length(g.names))
    names(g.label) = g.names
    g.label[names(partition.solved)] = as.character(partition.solved)
    g.label[duplicated(g.label)] = ''
    
    g %v% "labels" = g.label
    
    set.seed(seed.value)
    p <- ggnet2(
      g,
      label = 'labels',
      edge.color = "grey70",
      node.size = 1,
      color = "colors",
      palette = color.len,
    ) +  guides(size = F) + theme(aspect.ratio = 1) + 
      guides(color = guide_legend(title = "Length, bp"))
    
    name.output = sprintf("graph_%02d_labeled", i.plot)
    i.plot = i.plot + 1
    savePNG(p, path = path.figures, name = name.output,
            width = plot.size + 1, height = plot.size)
    if(flag.save.plot){
      saveRDS(p, file.path(path.figures, paste0(name.output, '.rds'))) 
    }
    
    rmSafe(p)
    rmSafe(g)
    rmSafe(g.names)
  }))
}


# 
# # Separate component
# 
# i.comp = 3
# sv.comp = names(partition.solved)[partition.solved == i.comp]
# edges.comp = edges.solved[edges.solved[,1] %in% sv.comp, ]
# 
# g <- network(edges.comp, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
# 
# g.names = network.vertex.names(g)
# g.names.len = as.numeric(sapply(g.names, function(s) strsplit(s, '\\|')[[1]][2]))
# g.names.len.bin <- as.character(cut(g.names.len, breaks = len.bins, right = FALSE, labels = len.labels))
# g %v% "colors" = as.character(g.names.len.bin)
# 
# set.seed(seed.value)
# p <- ggnet2(
#   g,
#   label = FALSE,
#   edge.color = "grey70",
#   node.size = 1,
#   color = "colors",
#   palette = color.len,
# ) +  guides(size = F) + theme(aspect.ratio = 1) + 
#   guides(color = guide_legend(title = "Length, bp"))
# 
# 
# p
# 
# name.output = 'tmp'
# savePNG(p, path = path.figures, name = name.output,
#         width = plot.size + 1, height = plot.size)
# 
# 
# 
# part.tmp = getGraphCommunities(edges.comp)
# g <- network(edges.comp, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
# g.names = network.vertex.names(g)
# 
# g %v% "partition" <- as.character(part.tmp[g.names])
# 
# set.seed(seed.value)
# p <- ggnet2(
#   g,
#   label = FALSE,
#   edge.color = "grey70",
#   node.size = 1,
#   color = "partition",
#   arrow.size = 2,
#   arrow.gap = 0.01
# ) + theme(legend.position = "none") + scale_color_viridis_d(name = "Partition")
# 
# name.output = 'tmp.part'
# savePNG(p, path = path.figures, name = name.output,
#         width = plot.size + 1, height = plot.size)
