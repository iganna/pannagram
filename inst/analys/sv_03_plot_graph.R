# Get SV positions, GFF files, dencity files and consensys sequences
# Find SVs and create GFF file

suppressMessages({ 
  library(rhdf5)
  library('foreach')
  library(doParallel)
  library("optparse")
  library(pannagram)
  library(crayon)
  library(ggplot2)
  
  library(igraph)
  # library(ggnet)
  library(network)
  library(GGally)
})



args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.cons"), type = "character", default = NULL, help = "path to directory with the consensus"),
  make_option(c("--cores"),     type = "integer",   default = 1, help = "number of cores to use for parallel processing")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# ***********************************************************************

flag.plot = F

# print(opt)

# ***********************************************************************
# Paths

if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop(paste0('Consensus folder does nto exist', path.cons))

path.sv = paste0(path.cons, 'sv/')
if (!dir.exists(path.sv)) dir.create(path.sv)
if(!dir.exists(path.sv)) stop(paste0('SV folder does nto exist', path.cons))


path.figures = paste0(path.cons, 'plot_svs/')
if (!dir.exists(path.figures)) dir.create(path.figures)
if(!dir.exists(path.figures)) stop(paste0('Folder for SV figures does nto exist', path.figures))

# ***********************************************************************
# ---- Values ----
len.min = 15 
sim.cutoff = 0.9
part.cutoff = 0.85

# Filtatio  of the graph
len.cutoff = 400
min.cl.size = 4
min.cnt = 2


# Binning
len.bins <- c(0, 100, 200, 400, 800, 1000, 3000, 5000, 7000, 12000, Inf)
len.labels <- c("0-100", "100-200", "200-400", "400-800", "800-1k", "1k-3k", "3k-5k", "5k-7k", "7k-12k", "12k+")

color.len <- c(
  "0-100" = "#1f77b4",
  "100-200" = "#50B498",
  "200-400" = "#2ca02c",
  "400-800" = "#bcbd22",
  "800-1k" = "#ff7f0e",
  "1k-3k" = "#d62728",
  "3k-5k" = "#9467bd",
  "5k-7k" = "#e377c2",
  "7k-12k" = "#8c564b",
  "12k+" = "#7f7f7f"
)

# ***********************************************************************
# ---- Reading the data ----
pokaz('Reading the data...')

file.sv.pos = paste0(path.sv, 'sv_pangen_pos.rds')
if(!file.exists(file.sv.pos)){
  pokazAttention('SVs were not generated.')
  quit(save = "no", status = 0)
}
sv.all = readRDS(file.sv.pos)
sv.all$chr = as.numeric(sv.all$chr)

sv.se = sv.all[sv.all$single == 1,]
sv.se$len.gr =  cut(sv.se$len, breaks = len.bins, right = FALSE, labels = len.labels)

f.max = max(sv.se$freq.max)

res.cover.file = paste0(path.sv, 'seq_sv_big_on_sv_cover.rds')
if(!file.exists(res.cover.file)){
  pokazAttention('All SVs are different, can not biuld a graph.')
  quit(save = "no", status = 0)
}
res.nest = readRDS(res.cover.file)

# ***********************************************************************
# ---- Make a collapsed graph ----
pokaz('Make a collapsed graph...')

file.g.content = paste0(path.sv, 'g_content_sim',round(sim.cutoff * 100),'.rds')

# remove duplicates in orientations
if(!file.exists(file.g.content)){
  idx = which(duplicated(res.nest[,c('V1', 'V8')]))
  if(length(idx) > 0){
    
    res.nest$comb = paste0(res.nest$V1, '|', res.nest$V8)
    res.nest$id = 1:nrow(res.nest)
    res.check = res.nest[res.nest$comb %in% res.nest$comb[idx],]
    res.check = res.check[order(res.check$comb),]
    cnt = table(res.check$comb)
    if(sum(cnt != 2) > 0) stop('Somethig is wrong with combinations')
    
    res.check$score = res.check$C1 + res.check$C8
    score.higher = res.check$score[seq(1,nrow(res.check), 2)] > res.check$score[seq(2,nrow(res.check), 2)]
    
    idx.remove = c(res.check$id[seq(1,nrow(res.check), 2)][!score.higher],
                   res.check$id[seq(2,nrow(res.check), 2)][score.higher])
    
    if(length(idx.remove) > 0){
      res.nest = res.nest[-idx.remove,,drop=F]
    }
  }
  
  g.content = getGraphFromBlast(res.nest = res.nest, sim.cutoff = sim.cutoff, collapse = T)
  
  saveRDS(g.content, file.g.content)
} else {
  g.content = readRDS(file.g.content)
}

g.content.init = g.content

# ***********************************************************************
# ---- Remain only the nodes with nore than 2 sequences ----

pokaz('Remain only the nodes with nore than 2 sequences...')

g.content = g.content.init
node.remain = g.content$nodes.traits$node[g.content$nodes.traits$cnt >=  min.cnt]
sv.remain = g.content$nodes$name[g.content$nodes$node %in% node.remain]

sv.remain.len = sapply(sv.remain, function(s) as.numeric(strsplit(s, '\\|')[[1]][2]))
sv.remain = sv.remain[sv.remain.len >= 100]

res.nest.remain = res.nest[(res.nest$V1 %in% sv.remain) & 
                             (res.nest$V8 %in% sv.remain),]

g.content = getGraphFromBlast(res.nest = res.nest.remain, sim.cutoff = sim.cutoff, collapse = F)

x = unique(c(g.content$edges))
length(x)

if(flag.plot){
  g <- network(g.content$edges, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
  b.graph.names = network.vertex.names(g)
  
  set.seed(239)
  p.refined <- ggnet2(g, label = F, edge.color = "black",
                      # node.size = g.nodes.cnt[b.graph.names],
                      node.size = 1,
                      color = '#468B97',
                      arrow.gap = 0.01, arrow.size = 2,
                      # color = g.nodes.fam[b.graph.names],
                      # palette = fam.palette,
                      # mode = "kamadakawai"
  )
  
  p.refined
  
  savePNG(p.refined, path = path.figures, name = 'graph_01_init', width = 6, height = 6)  
}


# ***********************************************************************
# ---- Filtration 2: Component size and length of nodes ----
pokaz('Filtration 2: Component size and length of nodes...')

# Components
g.comp <- getGraphComponents(g.content$edges)
memb.comp = g.comp$membership
g.content$edges.small = g.content$edges

## Remove some the smallest - this should be BEFORE THE LENGTH
idx.moderate = which(g.comp$csize >= min.cl.size)
sv.name.moderate = names(g.comp$membership)[g.comp$membership %in% idx.moderate]
g.content$edges.small = g.content$edges.small[(g.content$edges.small[,1] %in% sv.name.moderate),]

## Remove all short sequences, make partition of the rest, and then add all sequences "back" if they belong to only one component

edges = g.content$edges.small
sv.names = unique(c(edges))
sv.names.len = sapply(sv.names, function(s) as.numeric(strsplit(s, '\\|')[[1]][2]))

idx.short = (sv.names.len[edges[,1]] < len.cutoff) | (sv.names.len[edges[,2]] < len.cutoff)
edges.short = edges[idx.short,,drop=F]
g.content$edges.small = edges[!idx.short,,drop=F]

# ***********************************************************************
# ---- Construct the reduced graph ----

save(list = ls(), file = "tmp_workspace_graph.RData")

if(nrow(g.content$edges.small) > 0){
  pokaz('Construct the reduced graph...')
  
  ## Collapse small graph
  names.tmp = unique(c(g.content$edges.small))
  res.nest.small = res.nest[(res.nest$V1 %in% names.tmp) & (res.nest$V8 %in% names.tmp),]
  
  sum(!(res.nest.small$V8 %in% names.tmp))
  
  g.content.small = getGraphFromBlast(res.nest = res.nest.small, sim.cutoff = sim.cutoff, collapse = F)
  
  # Plot small components
  
  if(flag.plot){
    g.sm <- network(g.content.small$edges, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
    b.graph.names = network.vertex.names(g.sm)
    set.seed(239)
    p.sm <- ggnet2(g.sm, label = F, edge.color = "black",
                   node.size = 1,
                   color = '#468B97',
                   arrow.gap = 0.01, 
                   arrow.size = 2
    )
    
    p.sm
    
    savePNG(p.sm, path = path.figures, name = paste0('graph_02_refined_', min.cl.size, '_', len.cutoff), width = 6, height = 6)
  }
} else {
  g.content.small = g.content  # Kostyl
  g.content.small$edges = g.content.small$edges.small
}


# ***********************************************************************
# ---- Partition of the reduced graph ----
if(nrow(g.content.small$edges) > 0){
  pokaz('Partition of the reduced graph...')
  
  # I-graph
  edges <- g.content.small$edges 
  igraph_g <- graph_from_edgelist(as.matrix(edges), directed = TRUE)
  igraph_g <- as.undirected(igraph_g, mode = "collapse")
  
  # Louvain partition
  louvain_result <- cluster_louvain(igraph_g)
  
  # Get clusters
  partition <- membership(louvain_result)
  
  if(flag.plot){
    # Construct the graph
    g <- network(edges, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
    g.names = network.vertex.names(g)
    
    # Add cluster color
    g %v% "partition" <- as.character(partition[g.names])
    
    set.seed(239)
    p.partition <- ggnet2(
      g,
      label = FALSE,
      edge.color = "grey70",
      node.size = 1,
      color = "partition", 
      arrow.gap = 0.01,
      arrow.size = 2
    ) + theme(legend.position = "none") + scale_color_viridis_d(name = "partition")
    
    p.partition
    
    savePNG(p.partition, path = path.figures, name = paste0('graph_03_louvain'), width = 6, height = 6)
  }
} else {
  partition = c()
}

# ***********************************************************************
# ---- Filter edges because of the partition ----

if(nrow(edges) > 0){
  edges.keep <- edges[partition[edges[,1]] == partition[edges[,2]],]
  edges = edges.keep  
}



# ***********************************************************************
# ---- Put back all of the SVs ----
pokaz('Put back all of the SVs...')

# Add short sequences back is they are all connect to the same cluster

# sv.names = c(edges)  # Old version
sv.names = unique(g.content.init$nodes$name)

sv.names.short = setdiff(sv.names, c(edges))

sv.names.short.len = sapply(sv.names.short, function(s) as.numeric(strsplit(s, '\\|')[[1]][2]))
sv.names.short = sv.names.short[order(-sv.names.short.len)]

sv.names.short.len = sapply(sv.names.short, function(s) as.numeric(strsplit(s, '\\|')[[1]][2]))
sv.names.short = sv.names.short[sv.names.short.len > 100]

# Edges
edges.short = g.content$edges
edges.short = edges.short[(edges.short[,1] %in% sv.names.short) | 
                            (edges.short[,2] %in% sv.names.short),,drop=F]
edges.short = unique(edges.short)

# Variables
partition.add = partition
edges.add = edges
if(!is.null(partition.add)){
  max.part = max(partition.add)  
} else {
  max.part = 0
}

sv.names.ingraph = unique(c(edges.add))

for(s.sv in sv.names.short){
  # pokaz(s.sv)
  sv.edges.connect = edges.short[rowSums(edges.short == s.sv) == 1,, drop=F]
  sv.edges.connect = sv.edges.connect[((sv.edges.connect[,1] %in% sv.names.ingraph) + 
                                         (sv.edges.connect[,2] %in% sv.names.ingraph)) == 1,,drop=F]
  sv.node.connect = c(sv.edges.connect)
  sv.node.connect = setdiff(sv.node.connect, s.sv)
  sv.node.connect = intersect(sv.node.connect, names(partition.add))
  
  
  if(length(sv.node.connect) == 0) {
    partition.add[s.sv] = max.part + 1
    max.part = max.part + 1
    
    # Add node name to the graph
    sv.names.ingraph = c(sv.names.ingraph, s.sv)
    next
  }
  sv.node.connect.part = partition.add[sv.node.connect]
  part.cnt = table(sv.node.connect.part)
  if(max(part.cnt)/sum(part.cnt) > part.cutoff){
    i.part = as.numeric(names(part.cnt)[which.max(part.cnt)])
    partition.add[s.sv] = i.part
    sv.part = names(partition.add)[partition.add == i.part]
    sv.edges.add = sv.edges.connect[(sv.edges.connect[,1] %in% sv.part) & 
                                      (sv.edges.connect[,2] %in% sv.part),, drop=F]
    edges.add = rbind(edges.add, sv.edges.add)
    
    # Add node name to the graph
    sv.names.ingraph = c(sv.names.ingraph, s.sv)
  }
}

edges = edges.add

# ***********************************************************************
# ---- Filtration ----
pokaz('Filtration small components...')

### Remove small components
g.comp <- getGraphComponents(edges)
memb.comp = g.comp$membership

min.cl.size = 4
idx.moderate = which(g.comp$csize >= min.cl.size)
sv.name.moderate = names(g.comp$membership)[g.comp$membership %in% idx.moderate]
edges = edges[(edges[,1] %in% sv.name.moderate),]

if(nrow(edges) == 0){
  pokazAttention('No graph was generated, not enough SVs')
  quit(save = "no", status = 0)
}

# ***********************************************************************
# ---- Final partitioning ----

pokaz('Final partitioning...')

igraph_g <- graph_from_edgelist(as.matrix(edges), directed = TRUE)
igraph_g <- as.undirected(igraph_g, mode = "collapse")

# Clustering
louvain_result <- cluster_louvain(igraph_g)
partition <- membership(louvain_result)
partition = sort(partition)

saveRDS(partition, paste0(path.sv, 'sv_partition_solved.rds'))
write.table(as.matrix(partition), paste0(path.sv, 'sv_partition_solved.txt'), 
            row.names = T, col.names = F, sep = '\t', quote = F)



# ***********************************************************************
# ---- Save ----
## Save edges
saveRDS(edges, paste0(path.sv, 'edges_solved.rds'))
write.table(edges, paste0(path.sv, 'edges_solved.txt'), quote = F, sep = '\t', row.names = F, col.names = F)

if(!flag.plot){
  quit(save = "no", status = 0)
}

# ***********************************************************************
# ---- Plot ----

# Create graph
g <- network(edges, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
g.names = network.vertex.names(g)

g %v% "partition" <- as.character(partition[g.names])

set.seed(239)
p.refined <- ggnet2(
  g,
  label = FALSE,
  edge.color = "grey70",
  node.size = 1,
  color = "partition",
  arrow.gap = 0.01,
  arrow.size = 2
) + theme(legend.position = "none") + scale_color_viridis_d(name = "Partition")

p.refined

savePNG(p.refined, path = path.figures, name = paste0('graph_05_solved'),
        width = 6, height = 6)

# ***********************************************************************
# ---- Colored by length ----

g.names.len = as.numeric(sapply(g.names, function(s) strsplit(s, '\\|')[[1]][2]))
g.names.len.bin <- as.character(cut(g.names.len, breaks = len.bins, right = FALSE, labels = len.labels))
g %v% "colors" = as.character(g.names.len.bin)

set.seed(239)
p <- ggnet2(g,
            edge.color = "grey70", 
            node.size = 1,
            color = 'colors',
            palette = color.len
) + 
  guides(size = F)  + theme(aspect.ratio = 1) + 
  guides(color = guide_legend(title = "Length, bp"))

p

savePNG(p, path = path.figures, name = paste0('graph_06_colored'),
        width = 7, height = 6)


# ***********************************************************************
# ---- With labels ----


g.label <- as.character(partition[g.names])

g.label = rep('', length(g.names))
names(g.label) = g.names
g.label[names(memb.comp)] = as.character(memb.comp)
g.label[duplicated(g.label)] = ''

g %v% "labels" = g.label


set.seed(239)
p <- ggnet2(g, label = 'labels', 
            edge.color = "grey70", 
            node.size = 1,
            color = 'colors',
            palette = color.len
) + 
  guides(size = F)  + theme(aspect.ratio = 1) + 
  guides(color = guide_legend(title = "Length, bp"))
p


savePNG(p, path = path.figures, name = paste0('graph_07_label'),
        width = 7, height = 6)





