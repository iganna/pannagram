# This script gets gff file and the genome and return the set of candidate sequences for merging

source(system.file("utils/utils.R", package = "pannagram"))

pokazStage('Analyse counts..')

library(optparse)

option_list = list(
  make_option(c("--file.cnt"), type="character", default="", 
              help="Path to the GFF file", metavar="character"),
  make_option(c("--file.genome"), type="character", default="", 
              help="Path to the genome file", metavar="character"),
  make_option(c("--file.seqs"), type="character", default="", 
              help="Path to the sequences file", metavar="character"),
  make_option(c("--file.fix"), type="character", default="", 
              help="Path to the sequences file", metavar="character"),
  make_option(c("--copy.number"), type="integer", default=3, 
              help="Allowed minimal copu-number of genes", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# print(opt)

# ---- Parsing of parameters ----
opt = parse_args(opt_parser)

if (is.null(opt$file.cnt)) {
  stop("Error: --file.cnt is required.")
}
if (is.null(opt$file.genome)) {
  stop("Error: --file.genome is required.")
}
if (is.null(opt$file.seqs)) {
  stop("Error: --file.seqs is required.")
}
if (is.null(opt$file.fix)) {
  stop("Error: --file.fix is required.")
}

file.cnt = opt$file.cnt
file.genome = opt$file.genome
file.seqs = opt$file.seqs
file.fix = opt$file.fix
copy.number = opt$copy.number

# ---- Read the genome ----

if(!file.exists(file.genome)) stop('Genome file doesn’t exist')
genome = readFastaMy(file.genome)

genome.list = list()
for(i.chr in 1:length(genome)){
  genome.list[[i.chr]] = seq2nt(genome[i.chr])
}

pokaz('Chromosome lengths:', unname(nchar(genome)))

# ---- Sequence counts ----

res = read.table(file.cnt, 
                 row.names = 1, header = 1, stringsAsFactors = F)

res$name = rownames(res)
rownames(res) = NULL
res$id = as.numeric(sapply(res$name, function(s) strsplit(s, '\\|')[[1]][2]))
res$chr = as.numeric(gsub("Chr", '', sapply(res$name, function(s) strsplit(s, '\\|')[[1]][3])))

# Remove singletons
res = res[res$total >= copy.number,]
n.col.total = which(colnames(res) == 'total')
res = res[,-(1:(n.col.total - 1))]

# Sort according to the initial gff
res = res[order(res$id),]

# ---- Split merged and non-merged sequences ----
idx.merged = which(diff(res$id) == 1)
idx.merged = sort(unique(c(idx.merged, idx.merged + 1)))

# write.table(res[-idx.merged,], file.fix, append = T, quote = F, sep = '\t', col.names = F, row.names = F)
write.table(res, file.fix, append = T, quote = F, sep = '\t', col.names = F, row.names = F)  # Save all!

# ---- Merge further ----

if(length(idx.merged) == 0){
  pokaz('No further merging')
  quit(save = "no")
}
  
res = res[idx.merged,]

seqs.merge = c()
for(irow in 1:(nrow(res) - 1)){
  if((res$id[irow] + 1) != res$id[irow + 1]) next
  i.chr = as.numeric(gsub("Chr", "", strsplit(res$name[irow], '\\|')[[1]][3]))
  pos = sort(c(as.numeric(gsub("Chr", "", strsplit(res$name[irow], '\\|')[[1]][4:5])),
          as.numeric(gsub("Chr", "", strsplit(res$name[irow+1], '\\|')[[1]][4:5]))))
  pos1 = pos[1]
  pos2 = pos[4]
  
  s = nt2seq(genome.list[[i.chr]][pos1:pos2])
  s.name = paste0('te_merge|',
                  res$id[irow], '|',
                  'Chr', i.chr, '|', 
                  pos1,'|', 
                  pos2,'|', 
                  nchar(s))
  seqs.merge[s.name] = s
}

if(length(seqs.merge) > 0){
  pokaz('Total number of sequences:', length(seqs.merge))
  writeFastaMy(seqs.merge, file.seqs)  
} else {
  pokaz('No sequences were found for merging')  
}









