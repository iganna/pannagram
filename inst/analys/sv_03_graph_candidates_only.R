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
  make_option("--similarity", type = "integer", default = 85, help = "Similarity"),
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
sim.cutoff = opt$similarity
cov.cutoff = opt$coverage
suff.cutoffs = paste0(sim.cutoff, '_', cov.cutoff)

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
# ---- Get initial and solved edges ----
if(show.echo) pokaz('Get initial edges...')

nestedness = filterNestedness(nestedness,
                              min.len = min.len,
                              show.echo=show.echo)

edges.init = getGraphFromNestedness(nestedness, cov.cutoff = cov.cutoff)

edges.solved = readRDS(paste0(path.sv, 'edges_families_', suff.cutoffs, '.rds'))

components.lost = attributeNodes(edges.solved, edges.init)

# ***********************************************************************
# ---- Save ----

write.table(as.matrix(components.lost), 
            paste0(path.sv, 'sv_families_',suff.cutoffs,'_candidates.txt'), 
            row.names = T, col.names = F, sep = '\t', quote = F)
