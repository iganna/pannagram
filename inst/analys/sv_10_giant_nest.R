# Stitch BLAST hits of giant SV chunks back into the SVs and get the nestedness

suppressMessages({
  library(optparse)
  library(pannagram)
})
source(system.file("sim/sim_func.R", package = "pannagram"))

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option("--file.giant",  type = "character", default = NULL, help = "FASTA with giant SVs"),
  make_option("--file.blast",  type = "character", default = NULL, help = "BLAST of chunks on giant SVs"),
  make_option("--file.out",    type = "character", default = NULL, help = "Output nestedness table"),
  make_option("--similarity",  type = "numeric",   default = 85,   help = "Similarity threshold"),
  make_option("--coverage",    type = "numeric",   default = 85,   help = "Coverage threshold")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

file.giant = opt$file.giant
if(!file.exists(file.giant)) stop(paste0('No file with giant SVs!', file.giant))

file.blast = opt$file.blast
if(!file.exists(file.blast)) stop(paste0('No file with BLAST results!', file.blast))

file.out = opt$file.out
if(is.null(file.out)) stop('Output file is not specified')

sim.cutoff = opt$similarity
coverage = opt$coverage
pokaz('Similarity:', sim.cutoff, 'Coverage:', coverage)

# ***********************************************************************
# ---- Reading the data ----

sv.seqs = readFasta(file.giant)
len.sv = setNames(nchar(sv.seqs), names(sv.seqs))

v = readBlast(file.blast)
if(is.null(v) || nrow(v) == 0){
  pokazAttention('No similarity in giant SVs')
  quit(save = "no", status = 0)
}

# ---- Chunks back to SVs ----
# Chunk name: <SV name>|<position of the chunk start in the SV>
pos.chunk = as.numeric(sub('.*\\|', '', v$V1))
v$V1 = sub('\\|[0-9]+$', '', v$V1)
v$V2 = v$V2 + pos.chunk - 1
v$V3 = v$V3 + pos.chunk - 1
v$V9 = len.sv[v$V1]

v = v[v$V6 >= sim.cutoff,]
v = v[v$V1 != v$V8,]

if(nrow(v) == 0){
  pokazAttention('No similarity in giant SVs')
  quit(save = "no", status = 0)
}

# ---- Nestedness ----

res = findNestedness(v, use.strand = T)

pokaz('Number of pairs before the coverage cutoff', nrow(res))

res$len1 = len.sv[res$V1]
res$len8 = len.sv[res$V8]

res$p1 = res$C1 / res$len1
res$p8 = res$C8 / res$len8

res = res[,c('V1', 'V8', 'dir', 'len1', 'len8', 'p1', 'p8', 'B1', 'E1', 'B8', 'E8')]
colnames(res) <- c('name.query', 'name.target', 'strand', 'length.query', 'length.target', 'coverage.query', 'coverage.target',
                   'beg.query', 'end.query', 'beg.target', 'end.target')

# Incorporate coverage
res = res[(res$coverage.query > coverage / 100) | (res$coverage.target > coverage / 100),]
res = res[order(-res$coverage.query),]

pokaz('Number of pairs after the coverage cutoff', nrow(res))

write.table(res,
            file.out,
            sep       = "\t",
            quote     = FALSE,
            row.names = F)
