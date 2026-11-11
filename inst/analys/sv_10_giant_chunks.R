# Cut giant SVs into chunks for the self-BLAST

suppressMessages({
  library(optparse)
  library(pannagram)
})

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option("--file.giant",  type = "character", default = NULL, help = "FASTA with giant SVs"),
  make_option("--file.chunks", type = "character", default = NULL, help = "Output FASTA with chunks"),
  make_option("--chunk.len",   type = "integer",   default = 1000, help = "Length of chunks")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

file.giant = opt$file.giant
if(!file.exists(file.giant)) stop(paste0('No file with giant SVs!', file.giant))

file.chunks = opt$file.chunks
if(is.null(file.chunks)) stop('File with chunks is not specified')

chunk.len = opt$chunk.len
pokaz('Length of chunks', chunk.len)

# ***********************************************************************
# ---- Chunks ----

sv.seqs = readFasta(file.giant)
pokaz('Number of giant SVs', length(sv.seqs))

# Chunk name: <SV name>|<position of the chunk start in the SV>
seqs.chunks = c()
for(s.name in names(sv.seqs)){
  s.len = nchar(sv.seqs[s.name])
  pos.beg = seq(1, s.len, by = chunk.len)
  pos.end = c(pos.beg[-1] - 1, s.len)

  # Short tail goes to the previous chunk
  n.pos = length(pos.beg)
  if((n.pos > 1) && ((pos.end[n.pos] - pos.beg[n.pos] + 1) < chunk.len / 2)){
    pos.end[n.pos - 1] = s.len
    pos.beg = pos.beg[-n.pos]
    pos.end = pos.end[-n.pos]
  }

  s.chunks = substring(sv.seqs[s.name], pos.beg, pos.end)
  names(s.chunks) = paste(s.name, pos.beg, sep = '|')
  seqs.chunks = c(seqs.chunks, s.chunks)
}

pokaz('Number of chunks', length(seqs.chunks))
writeFasta(seqs.chunks, file.chunks)
