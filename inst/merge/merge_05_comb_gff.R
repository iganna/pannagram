# This script gets gff file and the genome and return the set of candidate sequences for merging

library(ggplot2)
library(pannagram)

library(optparse)

option_list = list(
  make_option("--path.out",        type="character", default="", help="Path to the output folder"),
  make_option("--file.gff",        type="character", default="", help="Initial gff file"),
  make_option("--file.gff.parent", type="character", default="", help="Gff file with parents")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# print(opt)

# ---- Parsing of parameters ----
opt = parse_args(opt_parser)

if (is.null(opt$path.out)) {
  stop("Error: --path.out is required.")
}
if (is.null(opt$file.gff)) {
  stop("Error: --file.gff is required.")
}

if (is.null(opt$file.gff.parent)) {
  stop("Error: --file.gff.parent is required.")
}


path.out = opt$path.out
file.gff = opt$file.gff
file.gff.parent = opt$file.gff.parent


pokaz('Initial GFF', basename(file.gff))
pokaz('Parents GFF', basename(file.gff.parent))

gff = read.table(file.gff, stringsAsFactors = F)

gff.merge = read.table(file.gff.parent, stringsAsFactors = F)
gff.merge = gff.merge[order(gff.merge$V5),]
gff.merge = gff.merge[order(gff.merge$V4),]
gff.merge = gff.merge[order(gff.merge$V1),]

pokaz(nrow(gff.merge))

n.register = nchar(as.character(nrow(gff.merge)))

id.merge = paste0('TE_merge_', sprintf(paste0("%0", n.register, "d"), 1:nrow(gff.merge)))

gff.merge$V9 = paste0("ID=", id.merge, ';', gff.merge$V9)

for(i.m in 1:nrow(gff.merge)){
  idx.merge = which((gff$V1 == gff.merge$V1[i.m]) & 
                      (gff$V4 >= gff.merge$V4[i.m]) & 
                      (gff$V5 <= gff.merge$V5[i.m]))
  if(length(idx.merge) < 2) stop(paste0("something is wrong with", i.m))
  
  gff$V9[idx.merge] = sub("Name=", paste0("Parent=", id.merge[i.m], ';Name='), gff$V9[idx.merge])
  
}

gff.out = rbind(gff[,1:9], gff.merge[,1:9])
gff.out = gff.out[order(-gff.out$V5),]
gff.out = gff.out[order(gff.out$V4),]
gff.out = gff.out[order(gff.out$V1),]


pokaz(basename(file.gff))

write.table(gff.out, 
            file = paste0(path.out, sub('.gff', '_merges.gff', basename(file.gff))),
            quote = F,
            row.names = F,
            col.names = F,
            sep = '\t')
  
  
  
  
  
  
  
  

