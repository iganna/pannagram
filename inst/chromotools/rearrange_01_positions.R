# Wind a match between chromosomes and get the blocks of interest

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(pannagram)
  library(crayon)
})

source(system.file("chromotools/rearrange_func.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("--path.aln"),       type = "character", default = NULL, help = "Path to the output directory with alignments"),
  make_option(c("--path.chr"),       type = "character", default = NULL, help = "Path to the output directory with chromosomes"),
  make_option(c("--path.processed"), type = "character", default = NULL, help = "Path to the output directory with processed chromosomes"),
  
  make_option(c("--ref"),            type = "character", default = NULL, help = "Name of the reference genome"),
  
  make_option(c("--cores"),          type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option(c("--path.log"),       type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"),      type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)


# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- Values of parameters ----

# Number of cores
num.cores <- opt$cores

path.aln       <- ifelse(!is.null(opt$path.aln), opt$path.aln, stop('Folder with Alignments is not specified'))
path.chr       <- ifelse(!is.null(opt$path.chr), opt$path.chr, stop('Folder with Chromosomes is not specified'))
path.processed <- ifelse(!is.null(opt$path.processed), opt$path.processed, stop('Output folder with processed genomes is not specified'))
ref            <- ifelse(!is.null(opt$ref), opt$ref, stop('Reference genome is not specified'))

checkDir(path.aln)
checkDir(path.chr)
checkDir(path.processed)

pokaz('Output folder with processed genomes:', path.processed)

# Create folders for the alignment results
if(!dir.exists(path.processed)) dir.create(path.processed)

# Accessions
pokaz('Folder with alignments:', path.aln)
files.aln <- list.files(path.aln, pattern = ".*_maj\\.rds$", full.names = F)
files.aln = sub("_maj.rds", "", files.aln)
accessions = c()
for(f in files.aln){
  s = strsplit(f, '_')[[1]]
  s = s[-(length(s)-(0:1))]
  s = paste0(s, collapse = '_')
  accessions = c(accessions, s)
  accessions = unique(accessions)
}

pokaz('Accessions:', accessions)

# ***********************************************************************
# ---- Test ----

# ref = 'GCA_002079055.1'
# path.aln = "~/Library/CloudStorage/OneDrive-Personal/iglab/projects/pannagram_meta/yeast_wild/alignment/alignments_GCA_002079055.1/"
# path.chr = "~/Library/CloudStorage/OneDrive-Personal/iglab/projects/pannagram_meta/yeast_wild/alignment/chromosomes/"
# acc = "GCA_002079175.1"
# 
# 
# ref = 'MN47'
# acc = '1741'
# path.aln = "~/Library/CloudStorage/OneDrive-Personal/iglab/projects/pannagram_meta/arabidopsis/alignment/"
# path.chr = "~/Library/CloudStorage/OneDrive-Personal/iglab/projects/pannagram_meta/arabidopsis/chromosomes/"


# ***********************************************************************
# ---- Variables ----

file.ref.len = paste0(path.chr, ref, '_chr_len.txt', collapse = '')
if(!file.exists(file.ref.len)) stop('File', file.ref.len, 'does not exist.')
ref.len = read.table(file.ref.len, header = 1)
min.overlap = 0.01
min.overlap.fragment = 0.01


# ***********************************************************************
# ---- Main loop ----
for(acc in accessions){
  pokaz('Accession', acc)
  
  file.acc.len = paste0(path.chr, acc, '_chr_len.txt', collapse = '')
  acc.len = read.table(file.acc.len, header = 1)
  
  files.aln <- list.files(path.aln, pattern = paste0(acc,".*_maj\\.rds$"), full.names = F)
  # pokaz(files.aln)
  files.sizes <- file.size(paste0(path.aln, files.aln))
  combinations = sub(paste0(acc, "_"), "", files.aln)
  combinations = sub("_maj.rds", "", combinations)
  
  result <- (do.call(rbind, strsplit(combinations, "_")))
  result <- data.frame(Col1 = as.numeric(result[, 1]), Col2 = as.numeric(result[, 2]))
  
  n.acc = max(result[,1])
  n.ref = max(result[,2])
  # pokaz(n.acc, n.ref)
  
  pokaz("Correspondence to genome", ref, "(reference)...")
  corresp.pure = c()
  corresp.combined = c()
  for(i.chr.ref in 1:n.ref){
    x.all = c()
    i.chr.ref.len = ref.len$len[i.chr.ref]
    min.len = min(round(i.chr.ref.len * min.overlap.fragment), 10000)
    
    for(i.chr.acc in 1:n.acc){
      
      file.aln = paste0(path.aln, acc, '_', i.chr.acc, '_', i.chr.ref, '_maj.rds')
      # pokaz(file.aln)
      if(!file.exists(file.aln)) next
      
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
      corresp.pure = rbind(corresp.pure, 
                           c(i.chr.ref, i.chr.acc))
    } else {
      pos = matrix(0, 
                   nrow = length(i.chr.corresp),
                   ncol = i.chr.ref.len)
      for(i in 1:length(i.chr.corresp)){
        i.chr.acc = i.chr.corresp[i]
        file.aln = paste0(path.aln, acc, '_', i.chr.acc, '_', i.chr.ref, '_maj.rds')
        # pokaz(file.aln)
        x = readRDS(file.aln)
        x = x[order(x$V7),]
        for(irow in 1:nrow(x)){
          pos[i,x$V4[irow]:x$V5[irow]] = x$V6[irow]
        }
      }
      
      # save(list = ls(), file = "tmp_workspace_1.RData")
      df.all = findBestChromosome (pos, 
                                   i.chr.ref.len, 
                                   i.chr.ref, 
                                   i.chr.corresp, 
                                   min.len) 
      
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
                             data.frame(beg = 1, 
                                        end = ref.len$len[i.chr.ref],
                                        len = ref.len$len[i.chr.ref],
                                        i.ref = i.chr.ref,
                                        i.acc = i.chr.acc
                                        ))
    
    
  }
  # print(corresp.combined)
  corresp.combined = corresp.combined[order(corresp.combined$i.acc),]
  # corresp.combined$id = 1:nrow(corresp.combined)
  
  i.chr.split = unique(corresp.combined$i.acc[duplicated(corresp.combined$i.acc)])
  corresp.remain = corresp.combined[!(corresp.combined$i.acc %in% i.chr.split), ]
  corresp.acc2ref = c()
  
  # Split accession chromosomes
  pokaz("Split accession",acc,"chromosomes...")
  if(length(i.chr.split) > 0){
    for(i.chr.acc in i.chr.split){
      
      i.chr.acc.len = acc.len$len[i.chr.acc]
      i.chr.corresp = corresp.combined$i.ref[corresp.combined$i.acc == i.chr.acc]
      
      pos = matrix(0, 
                   nrow = length(i.chr.corresp),
                   ncol = i.chr.acc.len)
      for(i in 1:length(i.chr.corresp)){
        i.chr.ref = i.chr.corresp[i]
        file.aln = paste0(path.aln, acc, '_', i.chr.acc, '_', i.chr.ref, '_maj.rds')
        x = readRDS(file.aln)
        x = x[order(x$V7),]
        for(irow in 1:nrow(x)){
          pos[i,x$V2[irow]:x$V3[irow]] = x$V6[irow]
        }
      }
      # save(list = ls(), file = "tmp_workspace_2.RData")
      df.all = findBestChromosome (pos, 
                                   i.chr.acc.len, 
                                   i.chr.acc, 
                                   i.chr.corresp, 
                                   min.len) 
      idx <- match(c('i.ref', 'i.acc'), colnames(df.all))
      colnames(df.all)[idx] <- c('i.acc', 'i.ref')
      
      corresp.acc2ref = rbind(corresp.acc2ref, df.all)
    }    
  }
  
  # Add the rest
  if(nrow(corresp.remain) != 0){
    for(irow in 1:nrow(corresp.remain)){
      i.chr.acc = corresp.remain$i.acc[irow]
      i.chr.ref = corresp.remain$i.ref[irow]
      i.chr.acc.len = acc.len$len[i.chr.acc]
      
      corresp.acc2ref = rbind(corresp.acc2ref,
                              data.frame(beg = 1, 
                                         end = i.chr.acc.len,
                                         len = i.chr.acc.len,
                                         i.acc = i.chr.acc,
                                         i.ref = i.chr.ref))
    }  # for
  }  # if
  
  
  pokaz("Get final correspondence...")
  corresp.acc2ref = corresp.acc2ref[order(corresp.acc2ref$i.ref),]
  corresp.acc2ref$pos = 0
  
  for(irow in 1:nrow(corresp.acc2ref)){
    i.chr.acc = corresp.acc2ref$i.acc[irow]
    i.chr.ref = corresp.acc2ref$i.ref[irow]
    # pokaz(irow, i.chr.acc, i.chr.ref)
    
    file.aln = paste0(path.aln, acc, '_', i.chr.acc, '_', i.chr.ref, '_maj.rds')
    x = readRDS(file.aln)
    x = x[(x$V2 >= corresp.acc2ref$beg[irow]) & (x$V3 <= corresp.acc2ref$end[irow]),]
    
    pos.mean = round((mean(x$V4) + mean(x$V5)) / 2)
    x$dir = (x$V4 > x$V5) * 1
    pos.plus = sum((x$V3 - x$V2 + 1) * (1 - x$dir))
    pos.minus = sum((x$V3 - x$V2 + 1) * x$dir)
    if(pos.minus > pos.plus){
      pos.mean = -pos.mean
    } 
    corresp.acc2ref$pos[irow] = pos.mean

  }
  
  # Save
  print(corresp.acc2ref)
  saveRDS(corresp.acc2ref, 
          paste0(path.processed, 'corresp_',acc,'_to_',ref, '.rds'))
  write.table(corresp.acc2ref, 
              paste0(path.processed, 'corresp_',acc,'_to_',ref, '.txt'), col.names = T, row.names = F, quote = F, sep = '\t')
  
}



# 
# 
# genome.ref = c()
# for(i.chr.acc in 1:n.acc){
#   genome.ref[i.chr.acc] = readFasta(paste0(path.chr, )) 
# }
# # # TODO: Set up the direction for every block
# # for(i.acc in i.chr.corresp){
# #   
# # }
# 
# 
# # Write chromosomes in a right order in a right merge


