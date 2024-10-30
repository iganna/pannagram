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
#   # print(h5ls(file.msa))
#   # v = h5read(file.msa, paste0(gr.accs.e, acc.new))
#   # h5write(v, file.msa, paste0(gr.accs.e, acc.prev))
# }

# ***********************************************************************
# ---- Merge all GFF files into a common structure ----
# Split gene and the rest (CDSs, exons....) -> save separately

# Initialize an empty list to store all GFF data

file.gff.main = paste0(path.res, 'gff_main.rds')
if(!file.exists(file.gff.main)){
  gff.main = c()
  
  # Loop through each accession in the list of accessions to work with
  for(acc in accessions){
    
    pokaz('Accession', acc)
    
    # Read the GFF file 
    gff = read.table(paste0(path.annot, acc,'.gff'), stringsAsFactors = F)
    colnames(gff) = c('V1', 'V2', 'type', 'beg', 'end', 'V6', 'strand', 'V8', 'info')
    n.gff = nrow(gff)
    gff = gff[gff$type != 'chromosome',]
    
    gff = gff[order(gff$beg),]
    gff = gff[order(gff$V1),]
    
    gff$idx = 1:nrow(gff)
    
    # ---- Check chromosome names format ----
    # If less than 70% - stop, maybe format is wrong.
    # And show examples which do not match with the pattern
    
    gff = extractChrByFormat(gff, s.chr)
    
    # Accession ID  
    gff$acc = acc
    
    # Separate annotations into non-overlapping categories
    
    gff.p = gff[gff$strand == '+',]
    gff.m = gff[gff$strand == '-',]
    
    res.p = findIncludeAndOverlap(gff.p)
    res.m = findIncludeAndOverlap(gff.m)
    
    idx.parental = rbind(res.p$idx.include,
                         res.m$idx.include)
    
    # Save
    gff.main = rbind(gff.main, res.p$gff.cut)
    gff.main = rbind(gff.main, res.m$gff.cut)
    
  }
  # Save the combined annotations to an RDS file for future use
  saveRDS(gff.main, file.gff.main)
} else {
  pokaz('Reading pre-calculations...')
  gff.main = readRDS(file.gff.main)
}

# ***********************************************************************
# ---- Convert of initial of GFF files ----
idx.remain = c()
for(acc in accessions){
  file.raw.gff = paste0(path.res, acc,'_pangen_raw.gff')
  file.raw.txt = paste0(path.res, acc,'_pangen_raw.txt')
  if(file.exists(file.raw.gff) & file.exists(file.raw.txt)) next
  pokaz('Conversion of accession', acc)
  gff.acc = read.table(paste0(path.annot, acc,'.gff'), stringsAsFactors = F)
  gff.acc = gff.main[gff.main$acc == acc,]
  gff.acc$idx.main = gff.acc$idx
  
  gff.acc.pan = gff2gff(path.cons = path.msa,
                        gff.acc,
                        acc1 = acc,
                        acc2 = s.pannagram,
                        n.chr = 5,
                        aln.type = aln.type,
                        s.chr = s.chr,
                        remain = T)
  
  idx.remain = c(idx.remain, gff.acc.pan$idx.main)
  
  write.table(gff.acc.pan, file.raw.txt,
              row.names = F, col.names = T, quote = F, sep = '\t')
  
  write.table(gff.acc.pan[,1:9], file.raw.gff,
              row.names = F, col.names = F, quote = F, sep = '\t')
  
}

file.gff.main.idx.remain = paste0(path.res, 'gff_main_idx_remain.rds')
if(length(idx.remain) != 0){
  saveRDS(idx.remain, file.gff.main.idx.remain)  
} else {
  idx.remain = readRDS(file.gff.main.idx.remain)
}

gff.main = gff.main[gff.main$idx %in% idx.remain,]

# ***********************************************************************

for(i.chr in 1:5){
  
  pokaz(paste('Chromosome', i.chr))
  file.chr.gff = paste0(path.res, 'chr_', i.chr,'.gff')
  if(file.exists(file.chr.gff)) next

  # Read basic info from the alignment
  s.comb = paste(i.chr, i.chr, sep = '_')
  file.msa = paste0(path.msa, aln.type, s.comb, ref.suff,'.h5')
  
  groups = h5ls(file.msa)
  # accessions = groups$name[groups$group == gr.accs.b]
  msa.lens = as.numeric(groups$dim[groups$group == gr.accs.b])
  len.pan = unique(msa.lens)
  if(length(len.pan) != 1) stop('something is wrong with the alignment LENGTHs')
  
  # ---- Gene Minus and gene Plus ----

  file.gene.blocks = paste0(path.res, 'gene_blocks_chr_', i.chr, '.RData')
    
  if(!file.exists(file.gene.blocks)){
    
    gene.mins = rep(0, len.pan)
    gene.plus = rep(0, len.pan)
    gene.blocks = c()
    
    for(acc in accessions){
      pokaz('Accession', acc)
      # Idx in the pangenome coordinate
      v = h5read(file.msa, paste0(gr.accs.e, acc))
      
      v = cbind(v, 1:len.pan)
      v = v[!is.na(v[,1]),]
      v = v[v[,1] != 0 ,]
      
      v[v[,1] < 0,] = v[v[,1] < 0,] * (-1)  # Fix inversion
      
      # Idx in the accession coordinate
      v.acc = rep(0, max(v[,1]))
      v.acc[v[,1]] = v[,2]
      len.acc = max(v[,1])
      
      # Gff acc
      gff.acc = gff.main[(gff.main$acc == acc) & (gff.main$chr == i.chr),]
      if(nrow(gff.acc) == 0) next
      gff.plus = gff.acc[gff.acc$strand == '+',]
      gff.mins = gff.acc[gff.acc$strand == '-',]
      
      # Genes, hoping that they are not overlapped
      g.plus = fillBegEnd(len.acc, gff.plus)
      g.mins = fillBegEnd(len.acc, gff.mins)
      g.plus[v.acc == 0] = 0
      g.mins[v.acc == 0] = 0
      
      # Change plus and minus depending on the strand
      idx.change = v.acc < 0
      pokaz('Changing the strand (bp):', sum(idx.change))
      if(sum(idx.change) > 0){
        g.tmp = g.plus[idx.change]
        g.plus[idx.change] = g.mins[idx.change]
        g.mins[idx.change] = g.tmp  
      }
      
      # Check that no split of genes
      g.id.split = setdiff(intersect(g.plus, g.mins), 0)
      if(length(g.id.split) > 0){
        pokazAttention('Number of split genes', length(g.id.split))
        g.plus[g.plus %in% g.id.split] = 0
        g.mins[g.mins %in% g.id.split] = 0
      }
      pokazAttention('Total number of lost genes', length(setdiff(gff.acc$idx, c(g.plus, g.mins))))
      
      # Blocks of genes
      g.plus.blocks = getGeneBlocks(g.plus, len.pan, v.acc)
      g.mins.blocks = getGeneBlocks(g.mins, len.pan, v.acc)
      
      # Fill pangenome coordinates with IDX of genes
      gene.plus.add = fillBegEnd(len.pan, g.plus.blocks)
      gene.mins.add = fillBegEnd(len.pan, g.mins.blocks)
      
      # Save it for Annogroups
      gene.plus = gene.plus + (gene.plus.add != 0) * 1
      gene.mins = gene.mins + (gene.mins.add != 0) * 1
      
      # Save begin-end for split annogroups
      g.plus.blocks$strand = '+'
      g.mins.blocks$strand = '-'
      blocks.df = rbind(g.plus.blocks, g.mins.blocks)
      blocks.df$acc = acc
      gene.blocks = rbind(gene.blocks, blocks.df)
    }  
    
    gene.list = list('+' = gene.plus, '-' = gene.mins)
    
    save(gene.list, gene.blocks, file = file.gene.blocks)
    rm(gene.plus)
    rm(gene.mins)
  } else {
    load(file.gene.blocks)
  }
  
  
  # save(list = ls(), file = "tmp_workspace_here.RData")
  # stop()

  # ----   Groups  ----
  
  pokaz('Groups')
  
  an.blocks.all = list()
  n.an = 0
  s.strand = names(gene.list)
  for(s.s in s.strand){
    # Find blocks
    an.blocks = findOnes((gene.list[[s.s]] > 0) * 1)
    an.blocks$strand = s.s
    an.blocks$idx = 1:nrow(an.blocks)
    
    # First big annogroups
    annogr = fillBegEnd(len.pan, an.blocks)
    gene.blocks.s = gene.blocks[gene.blocks$strand == s.s,]
    
    if(sum(annogr[gene.blocks.s$beg] !=  annogr[gene.blocks.s$end]) > 0) stop('Wrong annogroups assignment (#2)')
    
    # Assignment of big annogroups
    gene.blocks.s$annogr = annogr[gene.blocks.s$beg]
    if(min(gene.blocks.s$annogr) == 0) stop('Wrong annogroups assignment (#1)')
    
    # Find those which should be split
    cnt.acc = tapply(gene.blocks.s$acc, 
                     gene.blocks.s$annogr, 
                     function(x) max(table(x)))
    pokaz('Counts of genes in annogroups', table(cnt.acc))
    
    # Problemaric groups
    an.problem = as.numeric(names(cnt.acc)[cnt.acc != 1])
    for(i.gr in an.problem){
      
      # Define vector to determine splits
      gr.len = an.blocks$end[i.gr] - an.blocks$beg[i.gr] + 1
     
      # Get genes within this block
      genes.gr = gene.blocks.s[gene.blocks.s$annogr == i.gr,]
      genes.gr$beg = genes.gr$beg - an.blocks$beg[i.gr] + 1
      genes.gr$end = genes.gr$end - an.blocks$beg[i.gr] + 1
      
      genes.gr$len = genes.gr$end - genes.gr$beg
      # orfplot(genes.gr)
      
      # Define accessions having at least two genes
      acc.probl = unique(genes.gr$acc[duplicated(genes.gr$acc)])
      pos.split = rep(0, gr.len)
      for(acc.tmp in acc.probl){
        genes.gr.acc = genes.gr[genes.gr$acc == acc.tmp,]
        for(irow in 1:(nrow(genes.gr.acc) - 1)){
          # pokaz(acc.tmp, genes.gr.acc$end[irow],genes.gr.acc$beg[irow+1])
          if((genes.gr.acc$end[irow]+1) != genes.gr.acc$beg[irow+1]) {
            pos.split[(genes.gr.acc$end[irow]+1):(genes.gr.acc$beg[irow+1]-1)] = 
              pos.split[(genes.gr.acc$end[irow]+1):(genes.gr.acc$beg[irow+1]-1)] + 1
          }
        }
      }

      # Naive coverage. Was replaced.
      acc.all = unique(genes.gr$acc)
      mx.cover = matrix(0, nrow = length(acc.all), ncol = gr.len, dimnames = list(acc.all, NULL))
      for(acc.tmp in acc.all){
        genes.gr.acc = genes.gr[genes.gr$acc == acc.tmp,]
        for(irow in 1:nrow(genes.gr.acc)){
          mx.cover[acc.tmp, genes.gr.acc$beg[irow]:genes.gr.acc$end[irow]] = 1
        }
      }
      
    
      pos.cover = gene.list[[s.s]][an.blocks$beg[i.gr]:an.blocks$end[i.gr]]
      if(max((pos.cover + pos.split)) > length(accessions)) {
        pokazAttention('Problem with the split')
      }
      
      
      # ***********************************************************************
      
      # Forth
      mx.increase = which(colSums((mx.cover[,-1] > mx.cover[,-ncol(mx.cover)]) * 1) != 0)
      mx.decrease = which(colSums((mx.cover[,-1] < mx.cover[,-ncol(mx.cover)]) * 1) != 0)
      
      if(length(mx.increase) == 0) next
      if(length(mx.decrease) == 0) next
      
      mx = rbind(cbind(mx.increase, 1),
                 cbind(mx.decrease, -1))
      mx = mx[order(mx[,1]),]
      diff.beg = mx[which((mx[-1,2] == 1) & (mx[-nrow(mx),2] == -1)) + 1, 1]
      
      if(length(diff.beg) == 0) next
      
      
      an.blocks.split = data.frame(beg = c(1, diff.beg), 
                                   end = c(diff.beg - 1, length(pos.cover)))
      
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
                                                 (an.blocks.split$beg[i.bl+1]-1)]) > 0)
        # print(n.gap)
        if (2 * n.gap >= max(n.acc.block[i.bl], n.acc.block[i.bl + 1])){
          idx.merge = c(idx.merge, i.bl)
        }
      }
      
      if(length(idx.merge) > 0){
        an.blocks.split$end[idx.merge] = an.blocks.split$end[idx.merge + 1]
        an.blocks.split = an.blocks.split[(-1) * (idx.merge + 1),]
      }

      
      # ***********************************************************************
      # Plot
      # orfplot(genes.gr) +
      #   # geom_vline(xintercept = 6448, color = "green") +
      #   geom_vline(xintercept = an.blocks.split$beg, color = "red") +
      #   geom_vline(xintercept = an.blocks.split$end, color = "blue") 
        
      
      # ***********************************************************************
      if(nrow(an.blocks.split) > 1){
        
        an.blocks.split = an.blocks.split + an.blocks$beg[i.gr] +1
        an.blocks$end[i.gr] = an.blocks.split$end[1]
        an.blocks.split = an.blocks.split[-1,]
        an.blocks.split$strand = s.s
        
        an.blocks.split$idx = 0
    
        an.blocks = rbind(an.blocks, an.blocks.split)
        
      }
      
    }
    
    an.blocks = an.blocks[order(an.blocks$beg),]
    an.blocks$idx.prev = an.blocks$idx
    an.blocks$idx = (1:nrow(an.blocks)) + n.an
    n.an = nrow(an.blocks)
    
    an.blocks.all[[s.s]] = an.blocks
    
  } # s.s
  
  
  
  # ---- The rest of the annotation into groups ----
  
  
  # Assign exons to annotation groups and form gene genes
  
  for(acc in accessions){
    # Reading
    gff.exons = read.table(paste0(path.res, acc,'_pangen_raw.txt'), 
                           stringsAsFactors = F,
                           header = 1)
    gff.exons = gff.exons[gff.exons$V3 == 'CDS',]
    
    for(s.s in s.strand){
      
      pos.blocks = fillBegEnd(len.pan, an.blocks.all[[s.s]])
      
      gff.exons.s = gff.exons[(gff.exons$V7 == s.s) & (gff.exons$chr == i.chr),]
      
      gff.exons.s$an.beg = pos.blocks[gff.exons.s$V4]
      gff.exons.s$an.end = pos.blocks[gff.exons.s$V5]
      
      # To extend annotation
      # gff.exons.s$an.end[gff.exons.s$an.end == 0] = gff.exons.s$an.beg[gff.exons.s$an.end == 0]
      # gff.exons.s$an.beg[gff.exons.s$an.beg == 0] = gff.exons.s$an.end[gff.exons.s$an.beg == 0]
      
      gff.exons.s = gff.exons.s[gff.exons.s$an.beg != 0,]
      gff.exons.s = gff.exons.s[gff.exons.s$an.end != 0,]
      gff.exons.s = gff.exons.s[gff.exons.s$an.beg == gff.exons.s$an.end,]
      
      gff.exons.s = gff.exons.s[order(gff.exons.s$V4),]
      checkTranslocations(gff.exons.s$an.beg)
      
      # Create a new GFF annotation
      
      # Genious counts of IDs
      gff.exons.s$exon.id <- ave(gff.exons.s$an.beg, gff.exons.s$an.beg, FUN = seq_along)
      
      # Gff for exons - the same for both annotations
      gff.exons.s$V9 = paste('ID=',
                             'AT',i.chr,'Gr',which(s.strand == s.s),sprintf("%07.0f", gff.exons.s$an.beg),
                             '.', acc,
                             '.exon', sprintf("%02.0f",  gff.exons.s$exon.id),
                             ';Parent=' ,
                             'AT',i.chr,'Gr',which(s.strand == s.s),sprintf("%07.0f", gff.exons.s$group),
                             '.', acc,
                             sep = '')
      gff.exons.s$V2 = s.pannagram
      
      
      gff.mrna.an.gr = tapply(gff.exons.s$an.beg, gff.exons.s$an.beg, unique)
      gff.mrna = data.frame(V1 = gff.exons.s$V1[1],
                            V2 = s.pannagram,
                            V3 = 'gene',
                            V4 = tapply(gff.exons.s$V4, gff.exons.s$an.beg, min),
                            V5 = tapply(gff.exons.s$V5, gff.exons.s$an.beg, max),
                            V6 = '.',
                            V7 = s.s,
                            V8 = '.',
                            V9 = paste('ID=',
                                       'AT',i.chr,'Gr',which(s.strand == s.s),
                                       sprintf("%07.0f", gff.mrna.an.gr),
                                       '.', acc,
                                       sep = ''))
      
      
      gff.own = rbind(gff.mrna[,1:9], gff.exons.s[,1:9])
      gff.own = gff.own[order(gff.own$V4),]
      
      gff.mrna$V3 = 'mrna'
      gff.mrna$V9 = paste('ID=',
                          'AT',i.chr,'Gr',which(s.strand == s.s),
                          sprintf("%07.0f", gff.mrna.an.gr),
                          '.', acc,
                          ';Parent=' ,
                          'AT',i.chr,'Gr',which(s.strand == s.s),
                          sprintf("%07.0f", gff.mrna.an.gr),
                          sep = '')
      
      gff.pan = rbind(gff.mrna[,1:9], gff.exons.s[,1:9])
      gff.pan = gff.pan[order(gff.pan$V4),]
      
    }
    
    
  } # acc
  

}


# ---- Manual testing ----
# # TEST: Check presence of annogroups in previous annotation
# x = paste(an.blocks$beg, an.blocks$end, sep = '_')
# x = x[cnt.acc == 1]
# y = paste(pal$V4, pal$V5, sep = '_')
# 
# z = setdiff(x, y)


# ---- OLD CODE (2024) ----



#           ATTEMPTS FOR SOLVE CONFUSING GROUPS
# ***********************************************************************
# # Third 
# 
# n.acc.tmp = length(unique(genes.gr$acc))
# # pos.cover[pos.split == 0] = n.acc.tmp
# 
# # plot(pos.cover[-1] / pos.cover[-length(pos.cover)])
# diff.cover = pos.cover[-1] / pos.cover[-length(pos.cover)]
# diff.beg = sort(unique(c(which(diff.cover <= 0.5), 
#                          which(diff.cover >= 2) + 1)))
# 
# diff.beg = which(diff.cover >= 2) + 1
# 
# if(length(diff.beg) == 0) next
# 
# an.blocks.split = data.frame(beg = c(1, diff.beg), 
#                              end = c(diff.beg - 1, length(pos.cover)))
# 
# # Fis the "end" for every block
# for(i.bl in 1:(nrow(an.blocks.split) - 1)){
#   bl.range = an.blocks.split$end[i.bl]:(an.blocks.split$beg[i.bl]+1)
#   for(i.end in bl.range){
#     # stop()
#     if(sum(mx.cover[,i.end] != mx.cover[,i.end-1]) == 0) {
#       an.blocks.split$end[i.bl] = i.end - 1
#     } else {
#       break
#     }
#   }
# }
# 
# while(T){
#   idx = which(an.blocks.split$beg == an.blocks.split$end)
#   if(length(idx) == 0) break
#   idx = idx[1]
#   an.blocks.split$beg[idx + 1] = an.blocks.split$beg[idx]
#   an.blocks.split = an.blocks.split[-idx,]
# }


# ***********************************************************************
# # Second attempt
# n.acc.tmp = length(unique(genes.gr$acc))
# pos.cover[pos.split == 0] = n.acc.tmp
# an.blocks.split = findOnes((pos.cover/n.acc.tmp >= 0.5) * 1)


# ***********************************************************************
# First attempt
# # ratio.split.cover = pos.split / (pos.cover + pos.split)
# ratio.split.cover = pos.split / max((pos.cover + pos.split))
# 
# an.blocks.split = findOnes((ratio.split.cover <= 0.5) * 1)
# an.blocks.split = an.blocks.split + an.blocks$beg[i.gr] - 1





