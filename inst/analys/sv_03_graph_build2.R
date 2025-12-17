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
  make_option("--path.figures", type = "character", default = "", help = "Path to folder with figures"),
  make_option("--path.sv", type = "character", default = NULL, help = "Path to sv dir"),
  make_option("--file.nestedness", type = "character", default = NULL, help = "File with nestedness"),
  make_option("--coverage", type = "integer", default = 85, help = "Coverage"),
  make_option("--cores", type = "integer", default = 1, help = "Number of cores to use for parallel processing"),
  make_option("--flag.plot", type = "logical", default = TRUE, help = "Enable plotting (default: TRUE)")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# Paths

path.sv <- opt$path.sv
if(!dir.exists(path.sv)) stop(paste0('No SV dir!', path.sv))

path.figures <- opt$path.figures
if(!dir.exists(path.figures)) stop(paste0('No SV figures dir', path.figures))

file.nestedness <- opt$file.nestedness
if(!file.exists(file.nestedness)) stop(paste0('File with nestedness does not exist', file.nestedness))

# ***********************************************************************
# ---- Variables ----
cov.cutoff = opt$coverage
dominant.effect = 0.7

# Filtration of the graph
min.len = 200
min.copy = 2
min.comp.size = 3

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

# ***********************************************************************
# ---- Get initial edges ----
if(show.echo) pokaz('Get initial edges...')

# edges = res.sim[, c("name.query", "name.target")]
# edges = unique(edges)

nestedness = filterNestedness(nestedness,
                              min.len = min.len,
                              show.echo=show.echo)

edges.init = getGraphFromNestedness(nestedness, cov.cutoff = cov.cutoff)

# if(flag.plot){
if(F){
  g <- network(edges.init, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
  g.names = network.vertex.names(g)
  
  set.seed(seed.value)
  p <- ggnet2(g, label = F, edge.color = "black",
              node.size = 0.5,
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
nestedness.major = filterNestedness(nestedness,
                                    cov.cutoff = cov.cutoff,
                                    min.copy = min.copy)

edges.major = getGraphFromNestedness(nestedness.major, cov.cutoff = cov.cutoff)

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
    
    name.output = sprintf("graph_%02d_no_singletons", i.plot)
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

graph.compact = getGraphCompact(edges.major)
edges.compact = graph.compact$edges

if(nrow(edges.compact) > 0){
  
  if(T){
    # if(flag.plot){
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
      
      name.output = sprintf("graph_%02d_compact", i.plot)
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

edges.no.shortcut <- filterEdgesShortcut(edges.major)

# ***********************************************************************
# ---- Update the compact structure and visualize again  ----

graph.compact = getGraphCompact(edges.no.shortcut)
edges.compact = graph.compact$edges

if(nrow(edges.compact) > 0){
  if(T){
    # if(flag.plot){
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
      
      name.output = sprintf("graph_%02d_compact_no_shortcut", i.plot)
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
                                         seqs = seqs)

# ***********************************************************************
# ---- Update the compact structure and visualize again ----

graph.compact = getGraphCompact(edges.major.no.forks)
edges.compact = graph.compact$edges

if(nrow(edges.compact) > 0){
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
}


# ***********************************************************************
# ---- Solve Umbrella ----

edges.major.no.umbrella <- solveUmbrellaNodes(edges = edges.major.no.forks,
                                               seqs = seqs, show.echo=T)


# ***********************************************************************
# ---- Update the compact structure and visualize again ----

graph.compact = getGraphCompact(edges.major.no.umbrella)
edges.compact = graph.compact$edges

if(nrow(edges.compact) > 0){
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

edges.solved <- filterEdges(edges.solved, min.comp.size = min.comp.size)
edges.solved <- filterEdgesShortcut(edges.solved)
edges.solved <- solveForkNodes(edges = edges.solved, seqs = seqs)
edges.solved <- solveUmbrellaNodes(edges = edges.solved, seqs = seqs)

components.info = getGraphComponents(edges.solved)
components.solved = components.info$membership


# ***********************************************************************
# ---- Save ----

# Clustering safe
saveRDS(components.solved, paste0(path.sv, 'sv_families.rds'))
write.table(as.matrix(components.solved), paste0(path.sv, 'sv_families.txt'), 
            row.names = T, col.names = F, sep = '\t', quote = F)

## Save edges
saveRDS(edges.solved, paste0(path.sv, 'edges_families.rds'))
write.table(edges.solved, paste0(path.sv, 'edges_families.txt'), quote = F, sep = '\t', row.names = F, col.names = F)

if(!flag.plot){
  quit(save = "no", status = 0)
}

# ***********************************************************************
# ---- Plot ----

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
      color = "partition",
    ) + theme(legend.position = "none") + scale_color_viridis_d(name = "Partition")
    
    
    name.output = sprintf("graph_%02d_families", i.plot)
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
    
    name.output = sprintf("graph_%02d_families_colored", i.plot)
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
    
    name.output = sprintf("graph_%02d_families_labeled", i.plot)
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

