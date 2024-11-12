# Create annotation groups


# ***********************************************************************
# ---- Libraries and dependencies ----
library(crayon)
library(rhdf5)
source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func.R", package = "pannagram"))
source(system.file("analys/analys_func.R", package = "pannagram"))

# ***********************************************************************
# ---- Setup ----

# 
# # Set the working path for annotations
# path.msa = '/Users/annaigolkina/Library/CloudStorage/OneDrive-Personal/pushkin/genomes_GR3013/data/msa/'
# path.annot = '/Users/annaigolkina/Library/CloudStorage/OneDrive-Personal/pushkin/genomes_GR3013/data/gff/'
# path.res = '/Users/annaigolkina/Library/CloudStorage/OneDrive-Personal/pushkin/genomes_GR3013/data/gff_common/'
# if (!dir.exists(path.res)) {
#   dir.create(path.res, recursive = TRUE)
# } 

library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("--ref"),   type = "character", default = "",   help = "Reference prefix"),
  make_option(c("--aln.type"),   type = "character", default = NULL, help = "Alignment type"),
  make_option(c("--path.msa"),   type = "character", default = NULL, help = "Path to MSA files"),
  make_option(c("--path.annot"), type = "character", default = NULL, help = "Path to annotation files"),
  make_option(c("--path.res"),   type = "character", default = NULL, help = "Path to result files"),
  make_option(c("--s.chr"),      type = "character", default = '_Chr', help = "Chromosome name format")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))


# ***********************************************************************

# ---- HDF5 ----
source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

# ***********************************************************************


path.msa <- opt$path.msa
path.annot <- opt$path.annot
path.res <- opt$path.res
s.chr <- opt$s.chr
ref <- opt$ref
aln.type <- opt$aln.type


# Check that required parameters are provided
if (is.null(path.msa) || path.msa == "") stop("Error: 'path.msa' parameter must be specified.")
if (is.null(path.annot) || path.annot == "") stop("Error: 'path.annot' parameter must be specified.")
if (is.null(path.res) || path.res == "") stop("Error: 'path.res' parameter must be specified.")

# Set the defaul values of other parameters
if (ref == ''){
  ref.suff = ref
} else {
  ref.suff = paste0('_ref_',ref)
}

if(is.null(aln.type)){
  aln.type = aln.type.msa
}


# Create results directory if it does not exist
if (!dir.exists(opt$path.res)) {
  dir.create(opt$path.res, recursive = TRUE)
}


# ---- Variables ----

s.pannagram = 'Pannagram'
s.strand = c('+', '-')

# ---- Accessions ----

files.gff <- list.files(path = path.annot, pattern = "\\.gff$", full.names = FALSE)

accessions <- sub("\\.gff$", "", files.gff)

accessions = setdiff(accessions, '0')
accessions = setdiff(accessions, '22001')
accessions = accessions[1:27]

pokaz('Amount of accessions:', length(accessions))
pokaz('  ', accessions)


# ***********************************************************************
# Accessions 220011
# acc.new = "22001_mod"
# acc.prev = "220011"
# for(i.chr in 1:5){
#   file.msa = paste0(path.msa, aln.type, i.chr, '_', i.chr, '.h5')
#   print(h5ls(file.msa))
#   # v = h5read(file.msa, paste0(gr.accs.e, acc.new))
#   # h5write(v, file.msa, paste0(gr.accs.e, acc.prev))
# }

# ***********************************************************************
# ---- Merge all GFF files into a common structure ----
# Split gene and the rest (CDSs, exons....) -> save separately

# Initialize an empty list to store all GFF data

# ***********************************************************************
# ---- Convert of initial of GFF files ----

gff.main.pan = c()
gff.main.own = c()
for(acc in accessions){
  pokaz('Accession', acc)
  file.pan.gff = paste0(path.res, acc,'_pangen_raw.gff')
  file.acc.gff = paste0(path.annot, acc,'.gff')
  
  gff.acc = read.table(file.acc.gff, stringsAsFactors = F)
  gff.acc$acc = acc
  gff.acc$idx.init = 1:nrow(gff.acc)
  
  if(file.exists(file.pan.gff)){
    gff.acc.pan = read.table(file.pan.gff, stringsAsFactors = F)
    # gff.acc.pan$acc = acc
  } else {
    pokaz('Conversion of accession', acc)
    gff.acc.pan = gff2gff(path.cons = path.msa,
                          gff.acc,
                          acc1 = acc,
                          acc2 = s.pannagram,
                          n.chr = 5,
                          aln.type = aln.type,
                          s.chr = s.chr,
                          exact.match = T,
                          remain = T)
    
    # Save GFF
    # write.table(gff.acc.pan[,1:9], file.raw.gff, row.names = F, col.names = F, quote = F, sep = '\t')
    write.table(gff.acc.pan, file.pan.gff, row.names = F, col.names = F, quote = F, sep = '\t')
  }
  gff.main.pan = rbind(gff.main.pan, gff.acc.pan)
  gff.main.own = rbind(gff.main.own, gff.acc)
}

colnames(gff.main.pan) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "acc", "idx.init")
colnames(gff.main.own) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "acc", "idx.init")

gff.main.pan$acc = as.character(gff.main.pan$acc)
gff.main.own$acc = as.character(gff.main.own$acc)

gff.main.pan$chr = as.numeric(sub(paste0(s.pannagram, s.chr), '', gff.main.pan$V1))

idx.chr = grep(s.chr, gff.main.own$V1)
gff.main.own = gff.main.own[idx.chr,]

gff.main.own$chr = NA
for(acc in accessions){
  pokaz('Accession', acc)
  pattern.acc = paste0(acc, s.chr)
  
  idx.acc = gff.main.own$acc == acc
  gff.main.own$chr[idx.acc] = as.numeric(sub(pattern.acc, '', gff.main.own$V1[idx.acc]))
}
if(sum(is.na(gff.main.own$chr)) > 0) stop('Something is wrong with chromosomes')

# save(list = ls(), file = "tmp_workspace_all.RData")

# ***********************************************************************
# Get length of pangenome coordinates
len.pan = c()
for(i.chr in 1:5){
  file.msa = paste0(path.msa, aln.type, i.chr, '_', i.chr, '.h5')
  tmp = h5ls(file.msa)
  tmp = tmp[tmp$otype == "H5I_DATASET",]
  len.tmp = unique(tmp$dim)
  if(length(len.tmp) != 1){
    stop('Pangenome coordinated do not match')
  } else {
    len.tmp = as.numeric(len.tmp)
    len.pan[i.chr] = len.tmp
  }
}

# save(list = ls(), file = "tmp_workspace_easy.RData")
# stop()


# ***********************************************************************
# Groups
gff.gene.pan = gff.main.pan[gff.main.pan$V3 == 'gene',]
gff.gene.pan$group = 0
gff.gene.pan$idx = 1:nrow(gff.gene.pan)


gr.shift = 0
for(i.chr in 1:5){
  for(s.s in s.strand){
    pokaz('Chr', i.chr, s.s)
    idx.tmp = (gff.gene.pan$V7 == s.s) & (gff.gene.pan$chr == i.chr)
    gff.tmp = gff.gene.pan[idx.tmp,]
    
    gff.tmp = gff.tmp[order(gff.tmp$V4),]
    
    gff.tmp$group = 0
    gff.tmp$group[1] = 1
    gr.end = gff.tmp$V5[1]
    for(irow in 2:nrow(gff.tmp)){
      if(gff.tmp$V4[irow] < gr.end){
        gff.tmp$group[irow] = gff.tmp$group[irow-1]
        gr.end = max(gr.end, gff.tmp$V5[irow])
      } else {
        gff.tmp$group[irow] = gff.tmp$group[irow-1] + 1
        gr.end = gff.tmp$V5[irow]
      }
    }
    gff.tmp$group = gff.tmp$group + gr.shift
    gr.shift = max(gff.tmp$group)
    
    gff.gene.pan$group[gff.tmp$idx] = gff.tmp$group
    
    if(min(gff.tmp$group) == 0) stop('Groups with 0 index')
    pokaz('Number of groups', length(unique(gff.tmp$group)))
    
    # Checkup
    pos = rep(0,max(gff.tmp$V5) + 1)
    for(irow in 1:nrow(gff.tmp)){
      pos[gff.tmp$V4[irow]:gff.tmp$V5[irow]] = 1
    }
    tmp = findOnes(pos)
    if(nrow(tmp) != length(unique(gff.tmp$group))) stop('Groups are wrongly defined')
    
  }
}

if(sum(gff.gene.pan$group == 0) > 0) stop('Zero-groups gound')


# save(list = ls(), file = "tmp_workspace_anno2.RData")


# ***********************************************************************
# ---- Confusing groups ----

# Confusing groups
gff.gene.pan$group = paste0('gr_', gff.gene.pan$group)
gr.acc.cnt = tapply(gff.gene.pan$acc, gff.gene.pan$group, 
                    function(x){
                      x.tbl <- table(x)
                      n.dup <- sum(x.tbl > 1)
                      return(n.dup)
                    })

gr.confusing = names(gr.acc.cnt)[gr.acc.cnt != 0]

an.blocks = c()
for(s.gr in gr.confusing){
  gff.gr = gff.gene.pan[gff.gene.pan$group == s.gr, , drop=F]
  colnames(gff.gr)[1:9] = c('V1', 'V2', 'type', 'beg', 'end', 'V6', 'strand', 'V8', 'info')
  
  pos.beg = min(gff.gr$beg)
  pos.end = max(gff.gr$end)
  pos.shift = pos.beg - 1
  gff.gr$beg = gff.gr$beg - pos.shift
  gff.gr$end = gff.gr$end - pos.shift
  gr.len = pos.end - pos.beg + 1
  
  # orfplot(gff.gr, y = gff.gr$acc)
  
  # Alignment matrix
  acc.break <- c()  # Consider breaks only for those accessions, that are presented twice
  acc.all = unique(gff.gr$acc)
  mx.cover = matrix(0, nrow = length(acc.all), ncol = gr.len, dimnames = list(acc.all, NULL))
  for(acc.tmp in acc.all){
    gff.gr.acc = gff.gr[gff.gr$acc == acc.tmp,]
    for(irow in 1:nrow(gff.gr.acc)){
      mx.cover[acc.tmp, (gff.gr.acc$beg[irow]:gff.gr.acc$end[irow])  ] = 1
      
      if(irow > 1){
        acc.break = rbind(acc.break, 
                          data.frame(beg = gff.gr.acc$end[irow-1] + 1, 
                                     end = gff.gr.acc$beg[irow] - 1))
      }
    }
  }
  p.cover = colSums(mx.cover != 0) / length(acc.all)
  split.cutoff = 0.5
  
  
  # Previous method
  mx.increase = which(colSums((mx.cover[,-1] > mx.cover[,-ncol(mx.cover)]) * 1) != 0)
  mx.decrease = which(colSums((mx.cover[,-1] < mx.cover[,-ncol(mx.cover)]) * 1) != 0)
  
  if(length(mx.increase) == 0) stop('1')
  if(length(mx.decrease) == 0) stop('2')
  
  mx = rbind(cbind(mx.increase, 1),
             cbind(mx.decrease, -1))
  mx = mx[order(mx[,1]),]
  diff.beg = mx[which((mx[-1,2] == 1) & (mx[-nrow(mx),2] == -1)) + 1, 1]
  
  if(length(diff.beg) == 0)  stop('3')
  
  
  an.blocks.split = data.frame(beg = c(1, diff.beg), 
                               end = c(diff.beg - 1, ncol(mx.cover)))
  
  # Fit the "end" for every block
  for(i.bl in 1:(nrow(an.blocks.split) - 1)){
    bl.range = an.blocks.split$end[i.bl]:(an.blocks.split$beg[i.bl]+1)
    for(i.end in bl.range){
      if(i.end %in% mx.decrease) break
      # stop()
      if(sum(mx.cover[,i.end] != mx.cover[,i.end-1]) == 0) {
        an.blocks.split$end[i.bl] = i.end - 1
      } else {
        break
      }
    }
  }
  
  while(T){
    idx = which(an.blocks.split$beg == an.blocks.split$end)
    if(length(idx) == 0) break
    idx = idx[1]
    an.blocks.split$beg[idx + 1] = an.blocks.split$beg[idx]
    an.blocks.split = an.blocks.split[-idx,]
  }
  
  # Test every block and between blocks
  # Get number of accessions by each block
  n.acc.block = c()
  for(i.bl in 1:nrow(an.blocks.split)){
    n.acc.block[i.bl] = sum(rowSums(mx.cover[,an.blocks.split$beg[i.bl]:an.blocks.split$end[i.bl]]) > 0)
  }
  
  # Get number of accessions by each gap between blocks
  # Split or merge
  idx.merge = c()
  for(i.bl in 1:(nrow(an.blocks.split)-1)){
    n.gap = sum(rowSums(mx.cover[,(an.blocks.split$end[i.bl]+1):
                                   (an.blocks.split$beg[i.bl+1]-1),drop=F]) > 0)
    # print(n.gap)
    if (2 * n.gap >= max(n.acc.block[i.bl], n.acc.block[i.bl + 1])){
      idx.merge = c(idx.merge, i.bl)
    }
  }
  
  rownames(an.blocks.split) = NULL
  if(length(idx.merge) > 0){
    for(i.m in (idx.merge+1)){
      an.blocks.split$beg[i.m] = an.blocks.split$beg[i.m - 1]
    }
    an.blocks.split = an.blocks.split[-idx.merge,]
  }
  
  
  # ***********************************************************************
  
  # # Visualise
  # p = orfplot(gff.gr, y = gff.gr$acc) +
  #   annotate("rect", xmin = an.blocks.split$beg, xmax = an.blocks.split$end,
  #            ymin = -Inf, ymax = Inf, alpha = 0.2, fill = 'blue')
  # 
  # path.figures = '/Volumes/Samsung_T5/vienn/test/a27/intermediate/consensus/figures/'
  # png(paste(path.figures,'gr_', s.gr, '.png', sep = ''),
  #     width = 6, height = 6, units = "in", res = 300)
  # print(p)     # Plot 1 --> in the first page of PDF
  # dev.off()
  # 
  # ***********************************************************************

  an.blocks.split$beg = an.blocks.split$beg + pos.shift
  an.blocks.split$end = an.blocks.split$end + pos.shift
  an.blocks.split$chr = gff.gr$chr[1]
  an.blocks.split$group = gff.gr$group[1]
  an.blocks.split$strand = gff.gr$strand[1]
  
  an.blocks = rbind(an.blocks, an.blocks.split)
  
}

an.blocks.init = an.blocks


# ***********************************************************************
# ---- Non-confusing groups ----

idx.good = !(gff.gene.pan$group %in% gr.confusing)
pos.good.beg = tapply(gff.gene.pan$V4[idx.good], gff.gene.pan$group[idx.good], min)
pos.good.end = tapply(gff.gene.pan$V5[idx.good], gff.gene.pan$group[idx.good], max)
strand.good = tapply(gff.gene.pan$V7[idx.good], gff.gene.pan$group[idx.good], unique)
chr.good = tapply(gff.gene.pan$chr[idx.good], gff.gene.pan$group[idx.good], unique)
gr.good = group=names(pos.good.beg)

an.blocks.good = data.frame(beg = pos.good.beg[gr.good], 
                            end = pos.good.end[gr.good],
                            chr = chr.good[gr.good],
                            group = gr.good,
                            strand = strand.good[gr.good])

an.blocks.all = rbind(an.blocks, 
                      an.blocks.good)

# ***********************************************************************
# Check that groups are not overlapped

for(ichr in 1:5){
  for(s.s in s.strand){
    pokaz('Chr', i.chr, s.s)
    tmp = an.blocks.all[(an.blocks.all$chr == i.chr) &
                          (an.blocks.all$strand == s.s),]
    tot.len = sum(tmp$end - tmp$beg + 1)
    pos.blocks = fillBegEnd(len.pan[i.chr], tmp)
    
    if(sum(pos.blocks != 0) != tot.len) stop('Overlap')
  }
}


# ***********************************************************************
# ---- Exons ----

gff.exons.pan = gff.main.pan[gff.main.pan$V3 == 'exon',]
gff.exons.own = gff.main.own[gff.main.own$V3 == 'exon',]



# ***********************************************************************
# ---- Assign exons to annotation groups and form gene genes ----
gff.new.pan = c()
gff.new.own = c()

for(i.chr in 1:5){
  for(s.s in s.strand){
    for(acc in accessions){
      
      pokaz(i.chr, s.s, acc)
      # Reading
      gff.exons = gff.exons.pan[(gff.exons.pan$chr == i.chr) & 
                                (gff.exons.pan$acc == acc) &
                                (gff.exons.pan$V7 == s.s),]
      
      an.blocks.tmp = an.blocks.all[(an.blocks.all$chr == i.chr) &
                            (an.blocks.all$strand == s.s),]
      pos.blocks = fillBegEnd(len.pan[i.chr], an.blocks.tmp)

      gff.exons$an.beg = pos.blocks[gff.exons$V4]
      gff.exons$an.end = pos.blocks[gff.exons$V5]
      
      # To extend annotation
      # gff.exons$an.end[gff.exons$an.end == 0] = gff.exons$an.beg[gff.exons$an.end == 0]
      # gff.exons$an.beg[gff.exons$an.beg == 0] = gff.exons$an.end[gff.exons$an.beg == 0]
      
      gff.exons = gff.exons[gff.exons$an.beg != 0,]
      gff.exons = gff.exons[gff.exons$an.end != 0,]
      gff.exons = gff.exons[gff.exons$an.beg == gff.exons$an.end,]
      
      gff.exons = gff.exons[order(gff.exons$V4),]
      # checkTranslocations(gff.exons$an.beg)
      
      # Create a new GFF annotation
      
      # Genious counts of IDs
      gff.exons$exon.id <- ave(gff.exons$an.beg, gff.exons$an.beg, FUN = seq_along)
      
      # ---- Form pan annotation ----
      
      # Gff for exons - the same for both annotations
      gff.exons$V9 = paste('ID=',
                             'AT',i.chr,'Gr',which(s.strand == s.s),sprintf("%07.0f", gff.exons$an.beg),
                             '.', acc,
                             '.exon', sprintf("%02.0f",  gff.exons$exon.id),
                             ';Parent=' ,
                             'AT',i.chr,'Gr',which(s.strand == s.s),sprintf("%07.0f", gff.exons$an.beg),
                             '.', acc,
                             sep = '')
      gff.exons$V2 = s.pannagram
      
      gff.mrna.an.gr = tapply(gff.exons$an.beg, gff.exons$an.beg, unique)
      gff.mrna = data.frame(V1 = gff.exons$V1[1],
                            V2 = s.pannagram,
                            V3 = 'mRNA',
                            V4 = tapply(gff.exons$V4, gff.exons$an.beg, min),
                            V5 = tapply(gff.exons$V5, gff.exons$an.beg, max),
                            V6 = '.',
                            V7 = s.s,
                            V8 = '.',
                            V9 = paste('ID=',
                                       'AT',i.chr,'Gr',which(s.strand == s.s),
                                       sprintf("%07.0f", gff.mrna.an.gr),
                                       '.', acc,
                                       sep = ''))
      
      gff.tmp = rbind(gff.mrna, gff.exons[,1:9])
      gff.tmp = gff.tmp[order(gff.tmp$V4),]
      
      gff.new.pan = rbind(gff.new.pan, gff.tmp)
    
      # ---- Form the own genomes annotation
     
      
      indices <- match(gff.exons$idx.init, gff.exons.own$idx.init)
      gff.acc = gff.exons.own[indices,]
      gff.acc$exon.id = gff.exons$exon.id
      gff.acc$an.beg = gff.exons$an.beg
      
      gff.exons = gff.acc
      
      gff.exons$V9 = paste('ID=',
                           'AT',i.chr,'Gr',which(s.strand == s.s),sprintf("%07.0f", gff.exons$an.beg),
                           '.', acc,
                           '.exon', sprintf("%02.0f",  gff.exons$exon.id),
                           ';Parent=' ,
                           'AT',i.chr,'Gr',which(s.strand == s.s),sprintf("%07.0f", gff.exons$an.beg),
                           '.', acc,
                           sep = '')
      gff.exons$V2 = s.pannagram
      
      gff.mrna.an.gr = tapply(gff.exons$an.beg, gff.exons$an.beg, unique)
      gff.mrna = data.frame(V1 = gff.exons$V1[1],
                            V2 = s.pannagram,
                            V3 = 'mRNA',
                            V4 = tapply(gff.exons$V4, gff.exons$an.beg, min),
                            V5 = tapply(gff.exons$V5, gff.exons$an.beg, max),
                            V6 = '.',
                            V7 = s.s,
                            V8 = '.',
                            V9 = paste('ID=',
                                       'AT',i.chr,'Gr',which(s.strand == s.s),
                                       sprintf("%07.0f", gff.mrna.an.gr),
                                       '.', acc,
                                       sep = ''))
      
      gff.tmp = rbind(gff.mrna, gff.exons[,1:9])
      gff.tmp = gff.tmp[order(gff.tmp$V4),]
      
      gff.new.own = rbind(gff.new.own, gff.tmp)
      
    }
    
    
  }
}

write.table(gff.new.own, '/Volumes/Samsung_T5/vienn/test/a27/intermediate/consensus/annotation_new/gff_own.gff', 
            row.names = F,
            col.names = F,
            sep = '\t',
            quote = F)


write.table(gff.new.pan, '/Volumes/Samsung_T5/vienn/test/a27/intermediate/consensus/annotation_new/gff_pan.gff', 
            row.names = F,
            col.names = F,
            sep = '\t',
            quote = F)


for(acc in accessions){
  x = read.table(paste0('/Users/annaigolkina/Library/CloudStorage/OneDrive-Personal/vienn/pacbio/1001Gplus/01_data/04_annotation/02_pannagram/genes_v05/genes_v05_',acc,'.gff'))
  y = gff.new.own[grep(acc, gff.new.own$V1),]
  
  y = y[y$V3 == 'mRNA',]
  x = x[x$V3 == 'mRNA',]  
  pokaz()
}



tmp = y[,4:5]
pos.y = rep(0, 40000000)
for(irow in 1:nrow(tmp)){
  pos.y[tmp$V4[irow]:tmp$V5[irow]] = 1
}


tmp = x[,4:5]
# pos.x = rep(0, 40000000)
idx = c()
for(irow in 1:nrow(tmp)){
  # pos.x[tmp$V4[irow]:tmp$V5[irow]] = 1
  
  m = mean(pos.y[tmp$V4[irow]:tmp$V5[irow]])
  if(m != 1){
    idx = rbind(idx, c(irow, m))  
    # stop()
  }
  
}


table(pos.x, pos.y)




