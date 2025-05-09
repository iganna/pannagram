# Get SV positions, GFF files, dencity files and consensys sequences
# Find SVs and create GFF file

suppressMessages({ library(Biostrings)
  library(rhdf5)
  library('foreach')
  library(doParallel)
  library("optparse")
  library(pannagram)
  library(crayon)
  library(ggplot2)
  
  library(igraph)
  library(ggnet)
  library(network)
})


args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.cons"), type = "character", default = NULL, help = "path to directory with the consensus"),
  make_option(c("--cores"),     type = "integer",   default = 1, help = "number of cores to use for parallel processing")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

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

len.aa.min = 300

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

sv.partition = readRDS(paste0(path.sv, 'sv_partition_solved.rds'))
sv.seqs = readFasta(paste0(path.sv, 'seq_sv_big.fasta'))

sv.seqs = sv.seqs[names(sv.partition)]

pokaz('Numer of SVs is', length(sv.seqs))

# ***********************************************************************
# ---- Make a collapsed graph ----





















