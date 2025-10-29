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
  make_option("--path.sv", type = "character", default = NULL, help = "Path to sv dir"),
  make_option(c("--cores"),     type = "integer",   default = 1, help = "number of cores to use for parallel processing"),
  make_option(c("--len.orf.min"), type = "integer", default = 65, help = "Minimum Length of ORFs")  # ~200/3
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

path.sv = opt$path.sv
if(!dir.exists(path.sv)) stop(paste0('No SV dir!', path.sv))

# ***********************************************************************
# ---- Values ----

len.orf.min = opt$len.orf.min
pokaz('Minimum length of ORF is', len.orf.min)

len.sv.min = len.orf.min * 3

# ***********************************************************************
# ---- Reading the data ----

file.seqs = paste0(path.sv, 'seq_sv_large.fasta')

sv.seqs = readFasta(file.seqs)

sv.seqs = sv.seqs[nchar(sv.seqs) > len.sv.min]

pokaz('Number of SVs to analyse is', length(sv.seqs))

# ***********************************************************************
# ---- Get ORFs for every SV ----

orfs = c()
for(i.s in 1:length(sv.seqs)){
  s = sv.seqs[i.s]
  ref = orfFinder(s)
  orf.tmp = ref$orf[nchar(ref$orf) >= len.orf.min]
  if(length(orf.tmp) == 0) next
  names(orf.tmp) = paste(names(s), names(orf.tmp), sep='|')
  orfs = c(orfs, orf.tmp)
}

pokaz('Number of ORFs is', length(orfs))
writeFasta(orfs, paste0(path.sv, 'sv_large_orfs.fasta'))
