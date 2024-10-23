# This script gets gff file and the genome and return the set of candidate sequences for merging

library(optparse)

source(system.file("utils/utils.R", package = "pannagram"))

pokazStage('Extract candidates for merging..')


option_list = list(
  make_option(c("--file.gff"), type="character", default="", 
              help="Path to the GFF file", metavar="character"),
  make_option(c("--file.genome"), type="character", default="", 
              help="Path to the genome file", metavar="character"),
  make_option(c("--file.seqs"), type="character", default="", 
              help="Path to the sequences file", metavar="character"),
  make_option(c("--patterns"), type="character", default="", 
              help="Pattern to search for", metavar="character"),
  make_option(c("--len.max"), type="integer", default=50000, 
              help="Maximum length for filtering", metavar="integer"),
  make_option(c("--len.gap"), type="integer", default=1000, 
              help="Gap length for filtering", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# print(opt)

# ---- Parsing of parameters ----
opt = parse_args(opt_parser)

if (is.null(opt$file.gff)) {
  stop("Error: --file.gff is required.")
}
if (is.null(opt$file.genome)) {
  stop("Error: --file.genome is required.")
}
if (is.null(opt$file.seqs)) {
  stop("Error: --file.seqs is required.")
}

file.gff = opt$file.gff
file.genome = opt$file.genome
file.seqs = opt$file.seqs
patterns = opt$patterns
len.max = opt$len.max
len.gap = opt$len.gap

# ---- Prepare parameters ----

if(nchar(patterns) != ''){
  patterns = strsplit(patterns, ',')[[1]]
  pokaz('Keywords in types:', patterns)  
} else {
  pokaz('No keywords to filter types')  
}


# ---- Read the data ----
if(!file.exists(file.gff)) stop(paste('Gff file does not exist:', file.gff))

gff = read.table(file.gff, stringsAsFactors = F)
gff = gff[gff$V3 != 'centromeric_repeat',]

pokaz('Chromosomes:', unique(gff$V1))

gff$id = sapply(gff$V9, function(s){
  res = strsplit(s, ';')[[1]]
  res = paste0(res[1:3], collapse = '|')
  res = gsub('ID=', '', res)
  res = gsub('Name=', '', res)
  res = gsub('Classification=', '', res)
  return(res)
}  )

# ---- Remain "Parent" ----
idx.contains <- grepl("Parent", gff$V9)
if(sum(!idx.contains) > 0){
  gff = gff[!idx.contains,]  
}

# ---- Remove very long hits ----
idx.remain = (gff$V5 - gff$V4 + 1) <= len.max
if(sum(idx.remain) > 0){
  gff = gff[idx.remain,]
} else {
  stop(paste0('All hits are longer than ', len.max))
}

# ---- Remain types by the patterns ----
idx.remain = rep(F, nrow(gff))

if(nchar(patterns) != ''){  # If pattern is setup
  for(p in patterns){
    idx.remain = idx.remain | grepl(p, gff$V3)
  }
}

if(sum(idx.remain) > 0){
  gff = gff[idx.remain,]
  pokaz('Number of remained hits:', nrow(gff))
  print(table(gff$V3))
} else {
  stop('Patterns were not found in types of gff hits (3D column in the gff file)')
}

# ---- Sorting ----
rownames(gff) = NULL

# gff$idx.sort = 1:nrow(gff)

gff = gff[order(gff$V5),]
gff = gff[order(gff$V4),]
gff = gff[order(gff$V1),]

# is.unsorted(gff$idx.sort)
# which(diff(gff$idx.sort) != 1)
# gff$idx.sort = NULL

gff$idx = 1:nrow(gff)

# ---- Attributes of hits ----

gff$len = gff$V5 - gff$V4 + 1

gff$chr = as.numeric(gsub('Chr', '', gff$V1))

# ---- Read the genome ----

if(!file.exists(file.genome)) stop('Genome file doesn’t exist')
genome = readFastaMy(file.genome)

genome.list = list()
for(i.chr in 1:length(genome)){
  genome.list[[i.chr]] = seq2nt(genome[i.chr])
}

pokaz('Chromosome lengths:', unname(nchar(genome)))

# ---- Get all hits ----
# 
# seqs.all = c()
# for(i.chr in 1:5){
#   s.chr = genome.list[[i.chr]]
#   gff.chr = gff[gff$V1 == paste('Chr', i.chr, sep = ''),]
#   
#   seqs = c()
#   for(irow in 1:nrow(gff.chr)){
#     seqs[irow] = nt2seq(s.chr[gff.chr$V4[irow]:gff.chr$V5[irow]])
#   }
#   names(seqs) = paste(gff.chr$id, gff.chr$V1,sep = '|')
#   
#   names(seqs) = gff.chr$id
#   
#   seqs.all = c(seqs.all, seqs)
# }

# ---- Get sequences for merging ----

types = unique(gff$V3)
seqs.merge = c()
for(type in types){
  pokaz('Family', type)
  gff.type = gff[gff$V3 == type,]
  gff.type$dist = c(gff.type$V4[-1] - gff.type$V5[-nrow(gff.type)], 0)
  
  # Distances between Chromosomes set to Inf
  gff.type$dist[which(gff.type$V1[-1] != gff.type$V1[-nrow(gff.type)])] = Inf
  
  # Number of overlapping hits
  n.over = sum(gff.type$dist < 0)
  # pokaz('Number of overlaps:', n.over)
  
  for(i.chr in 1:5){
    pokaz('Chromosome', i.chr)
    s.chr = genome.list[[i.chr]]
    gff.chr = gff.type[gff.type$V1 == paste('Chr', i.chr, sep = ''),]
    
    idx.merge = which(gff.chr$dist <= len.gap)
    idx.merge = setdiff(idx.merge, nrow(gff.chr))
    pokaz('Number of merged hits:', length(idx.merge))
    
    if(length(idx.merge) == 0) next
    
    seqs = c()
    for(i.merge in idx.merge){
      pos1 = gff.chr$V4[i.merge]
      pos2 = gff.chr$V5[i.merge+1]
      if(pos2 < pos1) stop('Wrong positions')
      
      if(pos2 - pos1 > 50000) stop('length is too long')
      
      # if((pos2 == 14566046) && (pos1 == 14563556)){
      #   file.ws = "tmp_workspace.RData"
      # 
      #     all.local.objects <- ls()
      #     save(list = all.local.objects, file = file.ws)
      #   stop()
      # }
      
      seqs = c(seqs, nt2seq(s.chr[pos1:pos2]))
    }
    
    if(length(seqs) == 0) next
    
    # Names
    names(seqs) = paste0('te_merge|',
                         gff.chr$idx[idx.merge], '|',
                         gff.chr$V1[idx.merge],'|', 
                         gff.chr$V4[idx.merge],'|', 
                         gff.chr$V5[idx.merge+1],'|', 
                         nchar(seqs))
    seqs.merge = c(seqs.merge, seqs)
  }
}

if(length(seqs.merge) > 0){
  pokaz('Total number of sequences:', length(seqs.merge))
  writeFastaMy(seqs.merge, file.seqs)  
} else {
  pokaz('No sequences were found for merging')  
}






