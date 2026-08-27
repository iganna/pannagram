# Per-family slices of the nestedness table.
#
# The graph of families is built from the nestedness table, but the table itself
# is huge (millions of rows) and scanning it once per family is not an option.
# Here it is cut once into per-family pieces, so that everything downstream -
# the choice of the representatives, the orientation, the checks of the topology -
# reads only what belongs to its own family.
#
# The three inputs are the outputs of the previous steps, so the script can be run
# on a project that has already been processed, without rebuilding the graph.

suppressMessages({
  library(optparse)
  library(pannagram)
  library(data.table)
})

args = commandArgs(trailingOnly = TRUE)

option_list = list(
  make_option("--path.sv",    type = "character", default = NULL, help = "Path to the sv dir of the project"),
  make_option("--similarity", type = "integer",   default = 85,   help = "Similarity of the families"),
  make_option("--coverage",   type = "integer",   default = 85,   help = "Coverage of the families"),
  make_option("--path.out",   type = "character", default = NULL, help = "Where to write, default <path.sv>/graphs_<sim>_<cov>/")
)

opt = parse_args(OptionParser(option_list = option_list), args = args)

path.sv = opt$path.sv
if(is.null(path.sv) || !dir.exists(path.sv)) stop(paste0('No SV dir: ', path.sv))
suff = paste0(opt$similarity, '_', opt$coverage)

# The thresholds are in the name of the directory: a run with other -sim/-cov
# is a different graph, and it should not overwrite this one.
path.out = if(is.null(opt$path.out)) file.path(path.sv, paste0('graphs_', suff)) else opt$path.out
dir.create(path.out, showWarnings = FALSE, recursive = TRUE)
file.fam  = file.path(path.sv, paste0('sv_families_', suff, '.txt'))
file.edge = file.path(path.sv, paste0('edges_families_', suff, '.txt'))
file.nest = file.path(path.sv, paste0('nestedness_sv_large_', suff, '.txt'))
for(f in c(file.fam, file.edge, file.nest)) if(!file.exists(f)) stop(paste0('No file: ', f))

# ***********************************************************************
# ---- Reading ----

# The names of the SVs contain "|", so the separator has to be given explicitly:
# guessing it turns "SVgr_1_id_000572|238" into two columns.
fam = fread(file.fam, sep = '\t', header = FALSE, data.table = FALSE)
colnames(fam) = c('sv', 'family')
fam$family = as.character(fam$family)
pokaz('SVs in families:', nrow(fam), 'families:', length(unique(fam$family)))

edges = fread(file.edge, sep = '\t', header = FALSE, data.table = FALSE)[, 1:2]
colnames(edges) = c('v1', 'v2')
pokaz('Edges of the final graph:', nrow(edges))

nest = fread(file.nest, sep = '\t', header = TRUE, data.table = FALSE)
pokaz('Rows of the nestedness table:', nrow(nest))

# ***********************************************************************
# ---- Cutting ----

name2fam = setNames(fam$family, fam$sv)

f.q = name2fam[nest$name.query]
f.t = name2fam[nest$name.target]
keep = !is.na(f.q) & !is.na(f.t) & (f.q == f.t)
nest = nest[keep, , drop = FALSE]
nest$family = f.q[keep]
pokaz('Rows inside the families:', nrow(nest))

# Which of the pairs survived the cleaning of the graph. The edges are directed:
# `getGraphFromNestedness` puts query->target and target->query separately, each
# with its own condition. So both directions are marked, and everything
# downstream decides for itself what to call a mutual pair - both directions
# present, one of them, or something else.

by.fam = split(seq_len(nrow(nest)), nest$family)
pokaz('Writing', length(by.fam), 'files to', path.out)

info = data.frame(family = names(by.fam), n.members = 0L, n.pairs = 0L,
                  n.edge = 0L, n.both = 0L, n.sv.both = 0L, stringsAsFactors = FALSE)
n.mem = table(fam$family)

e.key = paste(edges$v1, edges$v2)
nest$in.graph = paste(nest$name.query, nest$name.target) %in% e.key
nest$in.graph.rev = paste(nest$name.target, nest$name.query) %in% e.key

for(i in seq_along(by.fam)){
  f.id = names(by.fam)[i]
  d = nest[by.fam[[i]], , drop = FALSE]
  write.table(d, file.path(path.out, paste0('fam_', f.id, '.txt')),
              sep = '\t', quote = FALSE, row.names = FALSE)

  both = d$in.graph & d$in.graph.rev
  info$n.members[i]  = n.mem[[f.id]]
  info$n.pairs[i]    = nrow(d)
  info$n.edge[i]     = sum(d$in.graph | d$in.graph.rev)
  info$n.both[i]     = sum(both)
  info$n.sv.both[i]  = length(unique(c(d$name.query[both], d$name.target[both])))
}

file.info.out = file.path(path.out, 'index.txt')
write.table(info, file.info.out, sep = '\t', quote = FALSE, row.names = FALSE)

pokaz('Index:', file.info.out)
pokaz('Families where no pair has both directions in the graph:', sum(info$n.both == 0))
pokaz('Share of members having a pair with both directions:',
      round(sum(info$n.sv.both) / sum(info$n.members), 3))
