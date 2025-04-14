# Wind a match between chromosomes and get the blocks of interest

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(pannagram)
  library(crayon)
})

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("--path.aln"),      type = "character", default = NULL, help = "Path to the output directory with alignments"),
  make_option(c("--path.resort"),   type = "character", default = NULL, help = "Path to the output directory with alignments"),
  
  make_option(c("--ref"),           type = "character", default = NULL, help = "Name of the reference genome"),
  
  make_option(c("--cores"),         type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option(c("--path.log"),      type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"),     type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Variables ----

max.len = 10^6
len.blast = 50

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- Values of parameters ----

# Number of cores
num.cores <- opt$cores

path.aln      <- ifelse(!is.null(opt$path.aln), opt$path.aln, stop('Folder with Alignments is not specified'))
path.resort   <- ifelse(!is.null(opt$path.resort), opt$path.resort, stop('Folder combinations to resort the genomes'))
base.acc      <- ifelse(!is.null(opt$ref), opt$ref, stop('Reference genome is not specified'))

pokaz('path resort', path.resort)

# Create folders for the alignment results
if(!dir.exists(path.resort)) dir.create(path.resort)

# Accessions
files.aln <- list.files(path.aln, pattern = ".*\\.rds$", full.names = F)
files.aln = sub("_maj.rds", "", files.aln)
accessions = c()
for(f in files.aln){
  s = strsplit(f, '_')[[1]]
  s = s[-(length(s)-(0:1))]
  s = paste0(s, collapse = '_')
  accessions = c(accessions, s)
  accessions = unique(accessions)
}

# ***********************************************************************
# ---- Correspondence ----

ref = 'GCA_002079055.1'
path.aln = "~/Library/CloudStorage/OneDrive-Personal/iglab/projects/pannagram_meta/yeast_wild/alignment/alignments_GCA_002079055.1/"
path.chr = "~/Library/CloudStorage/OneDrive-Personal/iglab/projects/pannagram_meta/yeast_wild/alignment/chromosomes/"
acc = "GCA_002079175.1"


ref = 'MN47'
acc = '1741'
path.aln = "~/Library/CloudStorage/OneDrive-Personal/iglab/projects/pannagram_meta/arabidopsis/alignment/"
path.chr = "~/Library/CloudStorage/OneDrive-Personal/iglab/projects/pannagram_meta/arabidopsis/chromosomes/"

file.ref.len = paste0(path.chr, ref, '_chr_len.txt', collapse = '')
ref.len = read.table(file.ref.len, header = 1)

min.overlap = 0.03
min.overlap.fragment = 0.01

for(acc in accessions){
  pokaz('Accession', acc)
  
  file.acc.len = paste0(path.chr, ref, '_chr_len.txt', collapse = '')
  acc.len = read.table(file.acc.len, header = 1)
  
  files.aln <- list.files(path.aln, pattern = paste0(acc,".*\\.rds$"), full.names = F)
  # pokaz(files.aln)
  files.sizes <- file.size(paste0(path.aln, files.aln))
  combinations = sub(paste0(acc, "_"), "", files.aln)
  combinations = sub("_maj.rds", "", combinations)
  
  result <- (do.call(rbind, strsplit(combinations, "_")))
  result <- data.frame(Col1 = as.numeric(result[, 1]), Col2 = as.numeric(result[, 2]))
  
  n.acc = max(result[,1])
  n.ref = max(result[,2])
  
  
  corresp.pure = c()
  corresp.combined = c()
  for(i.chr.ref in 1:n.ref){
    x.all = c()
    i.chr.ref.len = ref.len$len[i.chr.ref]
    for(i.chr.acc in 1:n.acc){
      file.aln = paste0(path.aln, acc, '_', i.chr.acc, '_', i.chr.ref, '_maj.rds')
      if(!file.exists(file.aln)) next
      # pokaz(file.aln)
      x = readRDS(file.aln)
      x$i.chr.acc = i.chr.acc
      x$dir = (x$V4 > x$V5) * 1
      x = glueZero(x)
      
      p = sum(abs(x$V5 - x$V4) + 1) / i.chr.ref.len
      if(p < min.overlap) next
      x.all = rbind(x.all, x)
      # mx.coverage[i.chr.ref, i.chr.acc] = sum(abs(x$V3 - x$V2) + 1)
    }
    i.chr.corresp = unique(x.all$i.chr.acc)
    if(length(i.chr.corresp) == 1){
      i.chr.acc = i.chr.corresp
      file.aln = paste0(path.aln, acc, '_', i.chr.acc, '_', i.chr.ref, '_maj.rds')
      x = readRDS(file.aln)
      x$dir = (x$V4 > x$V5) * 1
      pos.plus = sum((x$V3 - x$V2 + 1) * (1 - x$dir))
      pos.minus = sum((x$V3 - x$V2 + 1) * x$dir)
      if(pos.minus > pos.plus){
        i.chr.acc = -i.chr.acc
      }
      corresp.pure = rbind(corresp.pure, 
                           c(i.chr.ref, i.chr.acc))
    } else {
      pos = matrix(0, 
                   nrow = length(i.chr.corresp),
                   ncol = i.chr.ref.len)
      for(i in 1:length(i.chr.corresp)){
        i.chr.acc = i.chr.corresp[i]
        file.aln = paste0(path.aln, acc, '_', i.chr.acc, '_', i.chr.ref, '_maj.rds')
        x = readRDS(file.aln)
        x = x[order(x$V7),]
        for(irow in 1:nrow(x)){
          pos[i,x$V4[irow]:x$V5[irow]] = x$V6[irow]
        }
      }
      pos.idx = 1:ncol(pos)
      idx.remain = colSums(pos) > 0
      pos.idx = pos.idx[idx.remain]
      pos = pos[,idx.remain]
      
      min.len = min(round(i.chr.ref.len * min.overlap.fragment), 10000)

      df.all = c()
      for(i in 1:nrow(pos)){
        df = findOnes((pos[i,] != 0) * 1)
        df$len = df$end - df$beg + 1
        
        df$i.ref = i.chr.ref
        df$i.acc =  i.chr.corresp[i]
        
        df.all = rbind(df.all, df)
      }
      df.all = df.all[order(df.all$beg),]
      
      # Gradually remove the smallest
      while(min(df.all$len) < min.len){
        idx = which.min(df.all$len)[1]
        df.all = df.all[-idx,,drop=F]
        idx.merge = which(diff(df.all$i.acc) == 0)
        if(length(idx.merge) == 0) next
        idx.merge = idx.merge[1]
        df.all$end[idx.merge] = df.all$end[idx.merge + 1]
        df.all$len[idx.merge] = df.all$end[idx.merge] - df.all$beg[idx.merge] + 1
        df.all = df.all[-(idx.merge + 1),,drop=F]
        if(nrow(df.all) == 0) stop("no correspondence is left")
        # print(df.all)
      }
      df.all$beg = pos.idx[df.all$beg]
      df.all$end = pos.idx[df.all$end]
      df.all$len = df.all$end - df.all$beg + 1
      print(df.all)
      
      # Save
      corresp.combined = rbind(corresp.combined, df.all)
    }
  }
  colnames(corresp.pure) = c('i.ref', 'i.acc')
  corresp.pure = as.data.frame(corresp.pure)
  
  # Extend corresp.combined with the overlap
  for(i.cor in 1:nrow(corresp.pure)){
    i.chr.acc = abs(corresp.pure$i.acc[i.cor])
    i.chr.ref = corresp.pure$i.ref[i.cor]
    
    corresp.combined = rbind(corresp.combined,
                             c(1, ref.len$len[i.chr.ref], ref.len$len[i.chr.ref], i.chr.ref,  i.chr.acc))  
    
    
  }
  corresp.combined = corresp.combined[order(corresp.combined$i.acc),]
  corresp.combined
  
  # TODO: add corresponding coordinates of accessions genomes
  
  # # TODO: Set up the direction for every block
  # for(i.acc in i.chr.corresp){
  #   
  # }
  
  
  # Write chromosomes in a right order in a right merge
  
  
}




