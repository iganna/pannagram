# This script gets gff file and the genome and return the set of candidate sequences for merging

library(ggplot2)
library(pannagram)

library(optparse)

option_list = list(
  make_option("--path.out",      type="character", default="", help="Path to the output folder"),
  make_option("--file.genome",   type="character", default="", help="File with the reference genome"),
  make_option("--file.gff",      type="character", default="", help="Initial gff file"),
  make_option("--plot",          type="logical", default=TRUE, help="Enable plotting")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# print(opt)

# ---- Parsing of parameters ----
opt = parse_args(opt_parser)

if (is.null(opt$path.out)) {
  stop("Error: --path.out is required.")
}
if (is.null(opt$file.genome)) {
  stop("Error: --file.genome is required.")
}

if (is.null(opt$file.gff)) {
  stop("Error: --file.gff is required.")
}


path.out = opt$path.out
file.genome = opt$file.genome
file.gff = opt$file.gff
to.plot = opt$plot

if(!to.plot){
  pokaz("Visualisation is not needed")
}

# ---- Paths . Files ----

path.figures = paste0(path.out, 'figures/')
path.figures.m = paste0(path.figures, 'merges/')
if (!dir.exists(path.figures)) {
  dir.create(path.figures)
}
if (!dir.exists(path.figures.m)) {
  dir.create(path.figures.m)
}

path.res = paste0(path.out, 'simseqrch_seqs_fix/')

file.cnt = paste0(path.res, 'simsearch.GCA_028009825_80_95.cnt')

# path for Mafft
path.work = paste0(path.out, 'tmp/')
if (!dir.exists(path.work)) {
  dir.create(path.work)
}

path.gff.out = paste0(path.out, 'gff_out/')
if (!dir.exists(path.gff.out)) {
  dir.create(path.gff.out)
}

# ---- Read ----

tbl.cnt = read.table(file.cnt, stringsAsFactors = F)
tbl.cnt = tbl.cnt[tbl.cnt$total >= 4,]

pokaz("Number#1:", dim(tbl.cnt))

# Genome
genome = readFasta(file.genome)
genome.list = lapply(genome, seq2nt)
names(genome.list) = sapply(names(genome), function(s) strsplit(s, ' ')[[1]][1])


# ---- Read gff ----

gff = read.table(file.gff, stringsAsFactors = F)

gff.init = gff

categories2removed = c('centromeric_repeat',
                       'repeat_region',
                       'satellite_DNA',
                       'telomeric_repeat',
                       'cytosolic_25S_rRNA',
                       'cytosolic_5_8S_rRNA',
                       'cytosolic_18S_rRNA',
                       'cytosolic_5S_rRNA',
                       'rDNA_intergenic_spacer_element',
                       'repeat_fragment',
                       'rRNA_3_external_transcribed_sequence',
                       'rRNA_5_external_transcribed_sequence',
                       'rRNA_internal_transcribed_spacer1',
                       'rRNA_internal_transcribed_spacer2',
                       'tRNA_SINE_retrotransposon',
                       'target_site_duplication')

gff = gff[!(gff$V3 %in% categories2removed),]
table(gff$V3)

pokaz('Chromosomes:', unique(gff$V1))


gff$id = sapply(gff$V9, function(s){
  res = strsplit(s, ';')[[1]]
  res = paste0(res[1:3], collapse = '|')
  res = gsub('ID=', '', res)
  res = gsub('Name=', '', res)
  res = gsub('Classification=', '', res)
  return(res)
}  )


idx.contains <- grepl("Parent", gff$V9)
gff = gff[!idx.contains,]

idx.contains <- grepl("rDNA", gff$V9)
if(sum(idx.contains) != 0){
  stop('Something is wrong with rDNA')
}


# Remove very long hits
gff$len = (gff$V5 - gff$V4 + 1)
len.max = 50000
# len.min = 100
gff[gff$len > len.max,]

gff = gff[gff$len <= len.max,]
# gff = gff[gff$len >= len.min,]

# Result
pokaz('Number of hits:', nrow(gff))

# Sorting
rownames(gff) = NULL
gff$idx.sort = 1:nrow(gff)
gff = gff[order(gff$V5),]
gff = gff[order(gff$V4),]
gff = gff[order(gff$V1),]
# is.unsorted(gff$idx.sort)

gff$idx.sort = NULL
gff$idx = 1:nrow(gff)

# Chromosome ID
gff$chr = as.numeric(gsub('Chr', '', gff$V1))

# ---- Stat ----

m.names = rownames(tbl.cnt)

m.df <- data.frame(
  name = m.names, 
  chr = as.numeric(sapply(m.names, function(s) sub('Chr', '', strsplit(s, '\\|')[[1]][3]))), 
  beg = as.numeric(sapply(m.names, function(s) strsplit(s, '\\|')[[1]][4])), 
  end = as.numeric(sapply(m.names, function(s) strsplit(s, '\\|')[[1]][5]))
)

# ---- Check that blast hits for merged fragments don't overlap ----

file.merged.gff = paste0(path.res, 'simsearch.GCA_028009825_80_95.gff')

pokaz('File with results:', file.merged.gff)

gff.merge = read.table(file.merged.gff, stringsAsFactors = F)
idx.suspicious = c()

for(i.m in 1:nrow(m.df)){
  gff.tmp = gff.merge[grepl(m.df$name[i.m], gff.merge$V9, fixed = TRUE),]
  for(s.chr in unique(gff.tmp$V1)){
    gff.tmp.chr = gff.tmp[gff.tmp$V1 == s.chr,]
    if(nrow(gff.tmp.chr) == 1) next
    
    pos = c()
    for(jrow in 1:nrow(gff.tmp.chr)){
      pos = c(pos, gff.tmp.chr$V4[jrow]:gff.tmp.chr$V5[jrow])
    }
    
    # Corrected if statement
    if(sum(duplicated(pos)) != 0){
      idx.suspicious = c(idx.suspicious, i.m)
    }
  }
}

if(length(idx.suspicious) > 0){
  m.df = m.df[-idx.suspicious,]  
}

# ---- Number of merged ----

idx.merge = list()
for(i.m in 1:nrow(m.df)){
  idx.merge[[i.m]] = which((gff$chr == m.df$chr[i.m]) & (gff$V4 >= m.df$beg[i.m]) & (gff$V5 <= m.df$end[i.m]))
}

m.df$n = unlist(lapply(idx.merge, length))

table(m.df$n)

df <- as.data.frame(table(m.df$n))
df$Var1 = factor(df$Var1, levels = sort(unique(df$Var1)))

p = ggplot(df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = '#72BF78') +
  geom_text(aes(label = Freq), vjust = -0.5) +
  labs(x = "Number of merged fragments", y = "Amount") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  expand_limits(y = max(df$Freq) * 1.1) 


pdf(paste(path.figures, 'n_merged.pdf', sep = ''), width = 5, height = 4)
print(p)     # Plot 1 --> in the first page of PDF
dev.off()

# ---- Visualisation ----

colors <- colorRampPalette(c('#117554','#6EC207', "#FFEB00", "#4379F2"))
# colors <- colorRampPalette(c("#F87474", "#3AB0FF"))

check.again.out = c()
m.df$check = 0

for(i.m in which(m.df$n > 1)){
  pokaz(i.m)
  # for(i.m in 1:nrow(m.df)){
  gff.tmp = gff.merge[grepl(m.df$name[i.m], gff.merge$V9, fixed = TRUE),]

  idx.tmp = idx.merge[[i.m]]
  pos = sort(c(gff$V4[idx.tmp], gff$V5[idx.tmp]))
  
  seqs = nt2seq(genome.list[[m.df$chr[i.m]]][m.df$beg[i.m]:m.df$end[i.m]])
  seqs.names = c('init')
  for(irow in 1:nrow(gff.tmp)){
    s = genome.list[[gff.tmp$V1[irow]]][gff.tmp$V4[irow]:gff.tmp$V5[irow]]
    if(gff.tmp$V7[irow] == '-'){
      s = revCompl(s)
    }
    s = nt2seq(s)
    seqs = c(seqs, s)
    seqs.names = c(seqs.names, paste0('chr', 
                                      which(names(genome.list) == gff.tmp$V1[irow]), '|',
                                      gff.tmp$V4[irow], '|',
                                      gff.tmp$V5[irow], '|',
                                      gff.tmp$V5[irow] - gff.tmp$V4[irow] + 1))
  }
  names(seqs) = seqs.names
  
  file.clean = paste0(path.work, 'seqs_clean_', i.m, '.fasta')
  writeFasta(seqs, file.clean)
  
  aln.fasta = paste0(path.work, 'aln_', i.m, '.fasta')
  
  if(!file.exists(aln.fasta)){
    system(paste('mafft  --quiet --op 3  --ep 0.1 --treeout ', file.clean, '>', aln.fasta,  sep = ' '))
  }
  
  aln.mx = aln2mx(readFasta(aln.fasta))
  
  # save(list = ls(), file = "tmp_workspace2.RData")
  
  p = msadiff(aln.mx) 
  p.nt = msaplot(aln.mx) 
  
  # Positions in the alignment
  pos.aln = rep(0, ncol(aln.mx))
  pos.aln[aln.mx[1,] != '-'] = 1:nchar(seqs[1])
  
  # ---- Check positions ----
  pref = ""
  for(jrow in c(1,length(idx.tmp))){
    
    pre1 = gff$V4[idx.tmp[jrow]] - min(pos) + 1
    pre2 = gff$V5[idx.tmp[jrow]] - min(pos) + 1
    
    p1 = which(pos.aln == pre1)
    p2 = which(pos.aln == pre2)
    
    aln.mx.annot = aln.mx[,p1:p2]
    aln.mx.annot = aln.mx.annot[,aln.mx.annot[1,] != '-', drop = F]
    n.gap = colSums(aln.mx.annot == '-')
    p.gap = colSums(aln.mx.annot == '-') / nrow(aln.mx) 
    tot.gap = sum( p.gap > 0.5) / length(p.gap)
    if(tot.gap > 0.8) {
      pref = 'check_'
      check.again.out = c(check.again.out, i.m)
      if(jrow == 1){
        m.df$check[i.m] = m.df$check[i.m] + 1
      } else {
        m.df$check[i.m] = m.df$check[i.m] + 2
      }
    }
  }
  
  if(!to.plot){
    next
  }
  
  pokaz('Visualisation')
  
  # ---- Visualisation ----
  
  # Vertical annotation
  gradient_colors <- colors(length(idx.tmp))
  
  # Dotplot
  s = seqs[1]
  p.dot = dotplot.s(s, s, 15, 12)
  
  for(jrow in 1:length(idx.tmp)){
    
    pre1 = gff$V4[idx.tmp[jrow]] - min(pos) + 1
    pre2 = gff$V5[idx.tmp[jrow]] - min(pos) + 1
    
    p1 = which(pos.aln == pre1)
    p2 = which(pos.aln == pre2) 
    
    p = p +
      annotate("rect", xmin = p1, xmax = p2,
               ymin = -Inf, ymax = Inf, alpha = 0.2, fill = gradient_colors[jrow]) +
      annotate("rect", xmin = p1, xmax = p2, 
               ymin = 0, ymax = 0.5, alpha = 0.8,
               fill = gradient_colors[jrow])
    
    p.nt = p.nt +
      annotate("rect", xmin = p1, xmax = p2, 
               ymin = 0, ymax = 0.5,  alpha = 0.8,
               fill = gradient_colors[jrow])
    
    p.dot = p.dot + 
      annotate("rect", xmin = pre1, xmax = pre2 - 15,
               ymin = -Inf, ymax = Inf, alpha = 0.2, fill = gradient_colors[jrow])
    
  }
  
  p = p + ggtitle(paste0(c(paste(i.m, m.df$chr[i.m], m.df$beg[i.m], m.df$end[i.m] ),
                           paste(gff$id[idx.tmp], gff$chr[idx.tmp],gff$V4[idx.tmp], gff$V5[idx.tmp], sep = ' ')),
                         collapse = '\n')) + 
    theme(plot.title = element_text(size = 10)) 
  
  png(paste(path.figures.m, pref,  'aln_',i.m,'_msadiff.png', sep = ''), 
      width = 8, height = 6, units = "in", res = 300)
  print(p.nt)     # Plot 1 --> in the first page of PDF
  dev.off()
  
  png(paste(path.figures.m, pref, 'aln_',i.m,'_msaplot.png', sep = ''), 
      width = 8, height = 6, units = "in", res = 300)
  print(p)     # Plot 1 --> in the first page of PDF
  dev.off()
  
  png(paste(path.figures.m, pref, 'aln_',i.m,'_dot.png', sep = ''), 
      width = 8, height = 6, units = "in", res = 300)
  print(p.dot)     # Plot 1 --> in the first page of PDF
  dev.off()
  
}

# ---- Save ----

saveRDS(m.df, paste0(path.gff.out, 'm_df.rds'))

# ---- Form the GFF ----

m.df$id = 1:nrow(m.df)

idx.merge = list()
for(i.m in 1:nrow(m.df)){
  idx.merge[[i.m]] = which((gff$chr == m.df$chr[i.m]) & (gff$V4 >= m.df$beg[i.m]) & (gff$V5 <= m.df$end[i.m]))
}

idx.remove = which((m.df$check != 0) & (m.df$n == 2))
idx.remain = setdiff(which(m.df$check != 0), idx.remove)
for(i.m in idx.remain){
  idx.tmp = idx.merge[[i.m]]
  if(m.df$check[i.m] == 1){
    idx.tmp = idx.tmp[-1]
  } else {
    idx.tmp = idx.tmp[-length(idx.tmp)]
  }
  idx.merge[[i.m]] = idx.tmp
  gff.tmp = gff[idx.tmp,]
  
  beg.new = min(gff.tmp$V4)
  end.new = max(gff.tmp$V5)
  
  if(m.df$check[i.m] == 1){
    if(beg.new == m.df$beg[i.m]) stop("Wrong begin 1")
    if(end.new != m.df$end[i.m]) stop("Wrong end 1")
  } else {
    if(beg.new != m.df$beg[i.m]) stop("Wrong begin 2")
    if(end.new == m.df$end[i.m]) stop("Wrong end 2")
  } 
  m.df$beg[i.m] = beg.new
  m.df$end[i.m] = end.new
}

if(length(idx.remove) > 0){
  m.df = m.df[-idx.remove,]  
}
m.df = m.df[m.df$check != 3,]
rownames(m.df) = NULL


m.df$pref = ''
m.df$type = ''
m.df$dir = ''
for(i.m in 1:nrow(m.df)){
  idx.tmp = idx.merge[[i.m]]
  gff.tmp = gff[idx.tmp,]
  
  info = unname(sapply(gff.tmp$V9, function(s) sub("Name=", "", strsplit(s, ';')[[1]][2])))
  info = unique(info)
  info = sort(info)
  # if(length(info) != 1) stop()
  m.df$pref[i.m] = paste0(info, collapse = '|')
  
  # Type
  info = unique(gff.tmp$V3)
  info = sort(info)
  # if(length(info) != 1) stop()
  m.df$type[i.m] = paste0(info, collapse = '|')
  
  
  # Direction
  info = unique(gff.tmp$V7)
  info = sort(info)
  # if(length(info) != 1) stop()
  m.df$dir[i.m] = paste0(info, collapse = '|')
}

m.df$type[m.df$type == "Gypsy_LTR_retrotransposon|LTR_retrotransposon"] = "LTR_retrotransposon"
pokaz(unique(m.df$dir))
m.df$dir[!(m.df$dir %in% c('-', '+'))] = '.'
pokaz('Unique types:', unique(m.df$type))

saveRDS(m.df, paste0(path.gff.out, 'm_df.rds'))

genome.names = names(genome.list)

m.gff = data.frame(
                   # V1 = genome.names[m.df$chr],
                   V1 = paste0('Chr', m.df$chr),
                   V2 = 'merging', 
                   V3 = m.df$type,
                   V4 = m.df$beg,
                   V5 = m.df$end,
                   V6 = m.df$n,
                   V7 = m.df$dir,
                   V8 = '.',
                   V9 = paste0('Name=', m.df$pref))

write.table(m.gff, 
            file = paste0(path.gff.out, 'gff_merged.gff'),
            quote = F,
            row.names = F,
            col.names = F,
            sep = '\t')
