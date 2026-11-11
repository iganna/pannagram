# Families of SVs, as in sv_03_graph_build2.R, but every edge of the graph is
# checked to be end-to-end: the nested sequence should lie at one of the ends of
# the covering one, not just somewhere inside of it: the covered span should be
# not further than --end.gap from one of its ends. The check is applied when
# the edges are built, before any other filtration of the graph, and uses the
# covered spans from simsearch -positions, no dotplots.
# The dominant edges are selected before the shortcuts are removed: a shortcut
# is an edge bypassed by a longer path, and the path should not go along the
# edges which are to be removed as non-dominant right after that.
# The support of an edge is asymmetric (--dominant.to): most of the members of
# the nested node should be attached, but not to most of the covering ones,
# because only a part of them has the ends at the right place.

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
  make_option("--path.figures", type = "character", default = "", help = "Path to folder with figures"),
  make_option("--path.sv", type = "character", default = NULL, help = "Path to sv dir"),
  make_option("--file.nestedness", type = "character", default = NULL, help = "File with nestedness"),
  make_option("--similarity", type = "integer", default = 85, help = "Similarity"),
  make_option("--coverage", type = "integer", default = 85, help = "Coverage"),
  make_option("--cores", type = "integer", default = 1, help = "Number of cores to use for parallel processing"),
  make_option("--flag.plot", type = "logical", default = TRUE, help = "Enable plotting (default: TRUE)"),
  make_option("--suff.out", type = "character", default = "_ends", help = "Suffix for the output files"),
  make_option("--dominant.to", type = "double", default = 0, help = "Minimal fraction of the members of the covering (parent) node which should take part in an edge; for the nested (child) node it is always 0.7"),
  make_option("--end.gap", type = "integer", default = 100, help = "Maximal distance, bp, from the covered span to the nearest end of the covering sequence"),
  make_option("--path.code", type = "character", default = "", help = "Directory with the sources, if the package is not rebuilt yet")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# The package may be not rebuilt yet: take the functions from the repository
if(nzchar(opt$path.code)){
  p = file.path(opt$path.code, 'analys/graph_func.R')
  if(file.exists(p)) source(p)
}

# ***********************************************************************
# Paths

path.sv <- opt$path.sv
if(!dir.exists(path.sv)) stop(paste0('No SV dir!', path.sv))

path.figures <- opt$path.figures
if(!dir.exists(path.figures)) stop(paste0('No SV figures dir', path.figures))

file.nestedness <- opt$file.nestedness
pokaz(file.nestedness)

# ***********************************************************************
# ---- Variables ----
sim.cutoff = opt$similarity
cov.cutoff = opt$coverage

suff.cutoffs = paste0(sim.cutoff, '_', cov.cutoff)
suff.out = paste0(suff.cutoffs, opt$suff.out)
pokaz('Cutoff suffix', suff.cutoffs)

dominant.effect = 0.7

# Filtration of the graph
min.len = 200
min.copy = 2
min.comp.size = 3

# End-to-end check of the edges
end.gap = opt$end.gap
dominant.to = opt$dominant.to

# Plot variables
seed.value = 239
plot.size = 7

# Binning
source(system.file("analys/sv_variables.R", package = "pannagram"))

# Variables for plotting
flag.plot = opt$flag.plot
flag.save.plot = F
i.plot = 1

# Echo
show.echo = T

# ***********************************************************************
# ---- Reading the data ----
if(show.echo) pokaz('Reading the data...')

if(!file.exists(file.nestedness)){
  pokazAttention('All SVs are different, can not biuld a graph.')
  quit(save = "no", status = 0)
}

nestedness = read.table(file.nestedness, header = T)
seqs = readFasta(file.path(path.sv, 'seq_sv_large.fasta'))

# End-to-end scores (0/1) from the covered spans, instead of scoreFlankCoverage
flank.scores = getEndScores(nestedness, end.gap = end.gap)

# ***********************************************************************
# ---- Get initial edges ----
if(show.echo) pokaz('Get initial edges...')

nestedness = filterNestedness(nestedness,
                              min.len = min.len,
                              show.echo=show.echo)

edges.init = getGraphFromNestednessEnds(nestedness, cov.cutoff = cov.cutoff,
                                        end.gap = end.gap)

if(show.echo){
  pokaz('Number of edges by nestedness:', nrow(getGraphFromNestedness(nestedness, cov.cutoff = cov.cutoff)))
  pokaz('Number of end-to-end edges:', nrow(edges.init))
}

if(flag.plot && (nrow(edges.init) > 0)){
# if(F){
  g <- network(edges.init, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
  g.names = network.vertex.names(g)
  
  set.seed(seed.value)
  p <- ggnet2(g, label = F, edge.color = "black",
              node.size = 0.5,
              color = '#468B97',
              arrow.gap = 0.01, arrow.size = 2,
  )
  
  name.output = sprintf("graph%s_%02d_init", opt$suff.out, i.plot)
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
nestedness.major = filterNestedness(nestedness,
                                    cov.cutoff = cov.cutoff,
                                    min.copy = min.copy)

edges.major = getGraphFromNestednessEnds(nestedness.major, cov.cutoff = cov.cutoff,
                                         end.gap = end.gap)

if(nrow(edges.major) == 0){
  pokaz('Not echough SVs to build the graph')
  quit(save = "no")
}

if(flag.plot){
  suppressMessages(suppressWarnings({
    
    g <- network(edges.major, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
    
    set.seed(seed.value)
    p <- ggnet2(g, label = F, edge.color = "black",
                node.size = 0.5,
                color = '#468B97',
                arrow.gap = 0.01, arrow.size = 2,
    )
    
    name.output = sprintf("graph%s_%02d_no_singletons", opt$suff.out, i.plot)
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

# if(T){
if(flag.plot){
  
  if(show.echo) pokaz('Create a compact graph...')
  
  graph.compact = getGraphCompact(edges.major)
  edges.compact = graph.compact$edges
  
  if(nrow(edges.compact) > 0){
    
    
    suppressMessages(suppressWarnings({
      
      g <- network(edges.compact, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
      g.names = network.vertex.names(g)
      
      # Add cluster size
      g %v% "size" <- graph.compact$nodes.size[g.names]
      
      set.seed(seed.value)
      p <- ggnet2(g, label = F, edge.color = "black",
                  node.size = 'size',
                  color = '#468B97',
                  arrow.gap = 0.01, arrow.size = 2,
      ) + theme(legend.position = "none")
      p
      
      name.output = sprintf("graph%s_%02d_compact", opt$suff.out, i.plot)
      i.plot = i.plot + 1
      savePNG(p, path = path.figures, name = name.output,
              width = plot.size, height = plot.size)
      
      rmSafe(g)
      rmSafe(p)
      
    }))
  }
}

# ***********************************************************************
# ---- Remain only dominant edges ----

edges.no.dominant <- filterEdgesDominant(edges.major, dominant.effect.to = dominant.to)


# ***********************************************************************
## ---- Update the compact structure and visualize again  ----

# if(T){
if(flag.plot){
  
  graph.compact = getGraphCompact(edges.no.dominant)
  edges.compact = graph.compact$edges
  
  if(nrow(edges.compact) > 0){
    
    suppressMessages(suppressWarnings({
      
      g <- network(edges.compact, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
      g.names = network.vertex.names(g)
      
      # Add cluster size
      g %v% "size" <- graph.compact$nodes.size[g.names]
      
      set.seed(seed.value)
      p <- ggnet2(g, label = F, edge.color = "black",
                  node.size = 'size',
                  color = '#468B97',
                  arrow.gap = 0.01, arrow.size = 2,
      ) + theme(legend.position = "none")
      p
      
      name.output = sprintf("graph%s_%02d_compact_dominant", opt$suff.out, i.plot)
      i.plot = i.plot + 1
      savePNG(p, path = path.figures, name = name.output,
              width = plot.size, height = plot.size)
      
      rmSafe(g)
      rmSafe(p)
      
    }))
  }
}

# ***********************************************************************
# ---- Remove Shortcut edges ----

edges.no.shortcut <- filterEdgesShortcut(edges.no.dominant)

# ***********************************************************************
## ---- Update the compact structure and visualize again  ----
# if(T){
if(flag.plot){
  graph.compact = getGraphCompact(edges.no.shortcut)
  edges.compact = graph.compact$edges
  
  if(nrow(edges.compact) > 0){
    
    suppressMessages(suppressWarnings({
      
      g <- network(edges.compact, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
      g.names = network.vertex.names(g)
      
      # Add cluster size
      g %v% "size" <- graph.compact$nodes.size[g.names]
      
      set.seed(seed.value)
      p <- ggnet2(g, label = F, edge.color = "black",
                  node.size = 'size',
                  color = '#468B97',
                  arrow.gap = 0.01, arrow.size = 2,
      ) + theme(legend.position = "none")
      p
      
      name.output = sprintf("graph%s_%02d_compact_no_shortcut", opt$suff.out, i.plot)
      i.plot = i.plot + 1
      savePNG(p, path = path.figures, name = name.output,
              width = plot.size, height = plot.size)
      
      rmSafe(g)
      rmSafe(p)
      
    }))
  }
}

# ***********************************************************************
# ---- Solve forks ----

edges.major.no.forks <- solveForkNodes(edges = edges.no.shortcut,
                                       seqs = seqs, flank.scores = flank.scores,
  cutoff.remain.edges.to = dominant.to)

# ***********************************************************************
## ---- Update the compact structure and visualize again ----
# if(T){
if(flag.plot){
  
  graph.compact = getGraphCompact(edges.major.no.forks)
  edges.compact = graph.compact$edges
  
  if(nrow(edges.compact) > 0){
    
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
      
      name.output = sprintf("graph%s_%02d_compact_no_repeats", opt$suff.out, i.plot)
      i.plot = i.plot + 1
      savePNG(p, path = path.figures, name = name.output,
              width = plot.size, height = plot.size)
      
      rmSafe(g)
      rmSafe(p)
      
    }))
  }
}


# ***********************************************************************
# ---- Solve Umbrella ----

edges.major.no.umbrella <- solveUmbrellaNodes(edges = edges.major.no.forks,
                                              seqs = seqs, show.echo=T,
                                              flank.scores = flank.scores,
  cutoff.remain.edges.to = dominant.to)


# ***********************************************************************
## ---- Update the compact structure and visualize again ----
# if(T){
if(flag.plot){
  graph.compact = getGraphCompact(edges.major.no.umbrella)
  edges.compact = graph.compact$edges
  
  if(nrow(edges.compact) > 0){
    
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
      
      name.output = sprintf("graph%s_%02d_compact_no_umbrella", opt$suff.out, i.plot)
      i.plot = i.plot + 1
      savePNG(p, path = path.figures, name = name.output,
              width = plot.size, height = plot.size)
      
      rmSafe(g)
      rmSafe(p)
      
    }))
  }
}


# ***********************************************************************
# ---- Put back all of the SVs ----

# file.intermediate = paste0(path.sv, "tmp_workspace.RData")
# save(list = ls("edges.major.no.umbrella", "edges.init"), file = file.intermediate)

if(show.echo)  pokaz('Put all SVs back with dominant effect...')

edges.solved = putEdgesBack(edges = edges.major.no.umbrella, 
                            edges.init = edges.init,
                            dominant.effect = dominant.effect, 
                            show.echo=T)

pokaz('Remain Dominants')
edges.solved <- filterEdgesDominant(edges.solved, dominant.effect.to = dominant.to)
pokaz('Filter shortcuts')
edges.solved <- filterEdgesShortcut(edges.solved)
pokaz('Filter Forks')
edges.solved <- solveForkNodes(edges = edges.solved, seqs = seqs, flank.scores = flank.scores,
  cutoff.remain.edges.to = dominant.to)
pokaz('Filter Umbrella')
edges.solved <- solveUmbrellaNodes(edges = edges.solved, seqs = seqs, flank.scores = flank.scores,
  cutoff.remain.edges.to = dominant.to)
pokaz('Filter components')
edges.solved <- filterEdges(edges.solved, min.comp.size = min.comp.size)

if(nrow(edges.solved) == 0){
  pokaz('Not echough SVs to build the graph')
  quit(save = "no")
}

components.info = getGraphComponents(edges.solved)
components.solved = components.info$membership


# ***********************************************************************
# ---- Attribute the lost Nodes is possible ----

components.lost = attributeNodes(edges.solved, edges.init)

# ***********************************************************************
# ---- Save ----

# Clustering safe
saveRDS(components.solved, 
        paste0(path.sv, 'sv_families_',suff.out,'.rds'))
write.table(as.matrix(components.solved), 
            paste0(path.sv, 'sv_families_',suff.out,'.txt'), 
            row.names = T, col.names = F, sep = '\t', quote = F)
write.table(as.matrix(components.lost), 
            paste0(path.sv, 'sv_families_',suff.out,'_candidates.txt'), 
            row.names = T, col.names = F, sep = '\t', quote = F)

## Save edges
saveRDS(edges.solved, 
        paste0(path.sv, 'edges_families_',suff.out,'.rds'))
write.table(edges.solved, 
            paste0(path.sv, 'edges_families_',suff.out,'.txt'), quote = F, sep = '\t', row.names = F, col.names = F)

if(!flag.plot){
  quit(save = "no", status = 0)
}

# ***********************************************************************
# ---- Plot ----

# if(T){
if(flag.plot){
  suppressMessages(suppressWarnings({
    
    g <- network(edges.solved, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
    g.names = network.vertex.names(g)
    
    g %v% "partition" <- as.character(components.solved[g.names])
    
    set.seed(seed.value)
    p <- ggnet2(
      g,
      label = FALSE,
      edge.color = "grey70",
      node.size = 1,
      arrow.size = 5,
      color = "partition",
    ) + theme(legend.position = "none") + scale_color_viridis_d(name = "Partition")
    
    
    name.output = sprintf("graph%s_%02d_families", opt$suff.out, i.plot)
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

# if(T){
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
    
    name.output = sprintf("graph%s_%02d_families_colored", opt$suff.out, i.plot)
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

# if(T){
if(flag.plot){
  suppressMessages(suppressWarnings({
    
    g.label = rep('', length(g.names))
    names(g.label) = g.names
    g.label[names(components.solved)] = as.character(components.solved)
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
    
    name.output = sprintf("graph%s_%02d_families_labeled", opt$suff.out, i.plot)
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

