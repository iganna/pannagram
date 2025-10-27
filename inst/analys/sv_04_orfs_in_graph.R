# Get SV positions, GFF files, dencity files and consensys sequences
# Find SVs and create GFF file

suppressMessages({
  library(Biostrings)
  library(rhdf5)
  library(foreach)
  library(doParallel)
  library(optparse)
  library(pannagram)
  library(crayon)
  library(ggplot2)
  
  library(igraph)
  library(GGally)
  library(network)
})


args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.features.msa"), type = "character", default = NULL, help = "path to msa dir (features)"),
  make_option("--path.sv", type = "character", default = NULL, help = "Path to sv dir"),
  make_option(c("--path.cons"), type = "character", default = NULL, help = "path to directory with the consensus"),
  make_option(c("--cores"),     type = "integer",   default = 1, help = "number of cores to use for parallel processing")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

path.features.msa <- opt$path.features.msa
if(!dir.exists(path.features.msa)) stop(paste0('Consensus folder does nto exist', path.features.msa))

path.sv = opt$path.sv
if(!dir.exists(path.sv)) stop(paste0('No SV dir!', path.sv))

# ***********************************************************************
# ---- Values ----

len.aa.min = 200

# ***********************************************************************
# ---- Reading the data ----

file.partition = paste0(path.sv, 'sv_partition_solved.rds')
file.seqs = paste0(path.sv, 'seq_sv_large.fasta')

if(!file.exists(file.partition)){
  pokazAttention('No partitioning was generated, ORFs will not be generated')
  quit(save = "no", status = 0)
}

sv.partition = readRDS(file.partition)
sv.seqs = readFasta(file.seqs)

sv.seqs = sv.seqs[names(sv.partition)]

sv.seqs = sv.seqs[nchar(sv.seqs) >= len.aa.min * 3]

pokaz('Min length of ORFs:', len.aa.min)
pokaz('Min length of SVs:', len.aa.min * 3)

if(length(sv.seqs) == 0){
  pokaz('No SVs to analyse')
  quit(save = "no")
} else {
  pokaz('Number of SVs to analyse is', length(sv.seqs))
}

# ***********************************************************************
# ---- Get ORFs for every SV ----

orfs = c()
for(i.s in 1:length(sv.seqs)){
  s = sv.seqs[i.s]
  ref = orfFinder(s)
  orf.tmp = ref$orf[nchar(ref$orf) > len.aa.min]
  if(length(orf.tmp) == 0) next
  names(orf.tmp) = paste(names(s), names(orf.tmp), sep='|')
  orfs = c(orfs, orf.tmp)
}

pokaz('Numer of ORFs is', length(orfs))
writeFasta(orfs, paste0(path.sv, 'sv_in_graph_orfs.fasta'))
