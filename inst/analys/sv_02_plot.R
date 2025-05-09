# Get SV positions, GFF files, dencity files and consensys sequences
# Find SVs and create GFF file

suppressMessages({ library(Biostrings)
  library(rhdf5)
  library('foreach')
  library(doParallel)
  library("optparse")
  library(pannagram)
  library(crayon)
  library(ggplot2)
})


args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.cons"), type = "character", default = NULL, help = "path to directory with the consensus"),
  make_option(c("--cores"),     type = "integer",   default = 1, help = "number of cores to use for parallel processing")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# Paths

if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop(paste0('Consensus folder does nto exist', path.cons))

path.sv = paste0(path.cons, 'sv/')
if (!dir.exists(path.sv)) dir.create(path.sv)
if(!dir.exists(path.sv)) stop(paste0('SV folder does nto exist', path.cons))


path.figures = paste0(path.cons, 'plot_svs/')
if (!dir.exists(path.figures)) dir.create(path.figures)
if(!dir.exists(path.figures)) stop(paste0('Folder for SV figures does nto exist', path.figures))

# ***********************************************************************
# ---- Values ----
len.min = 15 


# Binning
len.bins <- c(0, 100, 200, 400, 800, 1000, 3000, 5000, 7000, 12000, Inf)
len.labels <- c("0-100", "100-200", "200-400", "400-800", "800-1k", "1k-3k", "3k-5k", "5k-7k", "7k-12k", "12k+")

color.len <- c(
  "0-100" = "#1f77b4",
  "100-200" = "#50B498",
  "200-400" = "#2ca02c",
  "400-800" = "#bcbd22",
  "800-1k" = "#ff7f0e",
  "1k-3k" = "#d62728",
  "3k-5k" = "#9467bd",
  "5k-7k" = "#e377c2",
  "7k-12k" = "#8c564b",
  "12k+" = "#7f7f7f"
)

# ***********************************************************************
# ---- Simple-Complex pie chart ----

file.sv.pos = paste0(path.sv, 'sv_pangen_pos.rds')
if(!file.exists(file.sv.pos)){
  stop('SVs were not generated.')
}
sv.all = readRDS(file.sv.pos)
sv.all$chr = as.numeric(sv.all$chr)

# Print stat
res.len = c()

thresholds = c(0,15,50,100, 1000)
for(thresh in thresholds){
  cnt = c(table(sv.all$single[sv.all$len >thresh]))
  tmp = c(sum(sv.all$len >thresh), cnt, as.numeric(sprintf("%.2f",cnt[2]/cnt[1])))
  res.len = rbind(res.len, tmp)
}
colnames(res.len) = c('all', 'complex', 'simple', 'ratio')
rownames(res.len) = paste('len >', thresholds, 'bp', sep = '')
rownames(res.len)[1] = 'all SVs'

print(res.len)


df = reshape2::melt(t(apply(res.len[,2:3], 1, function(x) x/sum(x))))
df$group = factor(rep(rownames(res.len), 2), levels = rev(rownames(res.len)))


# create nested pie chart using ggplot
p = ggplot(df, aes(x = factor(group), y = value, fill = factor(Var2))) +
  geom_col() +
  scale_x_discrete(limits = rev(unique(df$group))) +
  coord_polar("y")  + labs(fill = "") +
  labs(fill = "SVs:") +
  scale_fill_manual(values = c('#F48484', '#B5D5C5'),
                    labels = c("complex", "simple")) +
  theme_minimal() + 
  geom_text(data = df[1:length(thresholds),], aes(x = factor(group), y = 0, 
                                                  # label = gsub(">=", "≥ ",group)), 
                                                  label = group), 
            size = 3, hjust = -0.05, color = '#454545') +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # legend.position = "inside",
    legend.position.inside = c(0.04, 0.93),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  ) +
  geom_text(aes(label = round(value, 2)), position = position_stack(vjust = 0.5)) +
  ylab('')+ xlab('') + theme(plot.margin = unit(c(-1, -1, -1, -1), "cm"))

# p

pdf(paste0(path.figures, 'sv_pie_chart.pdf'), width = 4, height = 4)
print(p)     # Plot 1 --> in the first page of PDF
dev.off()

write.table(df, paste0(path.figures, 'sv_pie_chart.txt'), row.names = F, sep = '\t')
saveRDS(df, paste0(path.figures, 'sv_pie_chart.rds'))
saveRDS(res.len, paste0(path.figures, 'sv_pie_chart_num.rds'))


# ***********************************************************************
# ---- Chromosomal distribution ----

g <- ggplot(sv.all[sv.all$len > len.min,], aes(x=beg, fill = as.factor(single))) + 
  geom_histogram(bins = 50, color='grey20') + 
  theme_minimal() + 
  facet_grid(rows = vars(chr)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA),
        strip.text.y = element_text(angle = 0)) +
  labs(x = "Pangenome coordinate", y = "Count", fill = "SVs:") + 
  scale_fill_manual(values = c('#F48484', '#B5D5C5'), labels = c("complex", "simple"))

# g

pdf(paste(path.figures, 'sv_chr_minlen',len.min,'_pangen.pdf', sep = ''), 
    width = 6, height = 3/5 * max(sv.all$chr))
print(g)     # Plot 1 --> in the first page of PDF
dev.off()


# ***********************************************************************
# ---- Length distribution ----

sv.se = sv.all[sv.all$single == 1,]
sv.se$len.gr =  cut(sv.se$len, breaks = len.bins, right = FALSE, labels = len.labels)

cnt = as.matrix(table(sv.se$freq.max[sv.se$len>=len.min], sv.se$len.gr[sv.se$len>=len.min]))

f.max = max(sv.se$freq.max)

if(f.max != 1){
  cnt = rowSums(cnt)
  df = data.frame(Var1 = 1:length(cnt), value = cnt)
  g <- ggplot(df, aes(Var1, value)) +
    # annotate(geom = "rect",xmin = -Inf, xmax = 3, ymin = -Inf, ymax = Inf,
    #          fill = 'grey60', alpha = 0.5) +
    # annotate(geom = "rect",xmin = 25, xmax = Inf, ymin = -Inf, ymax = Inf,
    #          fill = 'grey60', alpha = 0.5) +
    geom_line(size = 2) + theme_minimal() + 
    theme(legend.position='none',
          strip.text.y = element_text(angle = 0)) + 
    # geom_segment(aes(x = 1, y = 15000, xend = f.max, yend = 15000), 
    #              arrow = arrow(type = "closed", length = unit(0.2, "cm"), ends = "both"), 
    #              color = "grey20") +
    # annotate("text", x = 2.5, y = 14000, label = "Singleton\ninsertions", hjust = "left", vjust = 1, color = "grey20", 
    #          size = 2.5) +
    # annotate("text", x = 25.5, y = 14000, label = "Singleton\ndeletions", hjust = "right", vjust = 1, color = "grey20",
    #          size = 2.5)+
    # viridis::scale_color_viridis() +
    xlab('Frequency of presence') + 
    ylab('Absolute number') +
    theme(
      panel.background = element_rect(fill = "white", color = 'white'),
      plot.background = element_rect(fill = "white", color = 'white')
    ) 
  # g  
  
  pdf(paste(path.figures, 'sv_freq_hist.pdf', sep = ''), width = 2.6, height = 1.7)
  print(g)     # Plot 1 --> in the first page of PDF
  dev.off()
}


# ***********************************************************************
# ---- Length distribution bins ----

tbl = table(sv.se[(sv.se$len > len.min), c('len.gr', 'freq.max')])
tbl = tbl[rowSums(tbl) != 0,,drop=F]
tbl = apply(tbl, 2, function(x) x / sum(x))
df = reshape2::melt(tbl)
p <- ggplot(data=df, aes(x=freq.max, y = value, fill=len.gr)) +
  geom_bar(stat="identity") + theme_minimal()  + xlab('Frequency of presence') + ylab('Proportion') + 
  scale_fill_manual(values=color.len,  name ='sSVs len (bp)') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))  +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.direction="vertical",
        legend.box="horizontal",
        legend.box.just = c("top"), 
        legend.background = element_rect(fill=alpha('white', 0.75)),
        legend.margin = margin(2, 2, 2, 2))
p
# 

pdf(paste0(path.figures, 'sv_freq_hist_length_minlen', len.min ,'_norm.pdf'), 
    width = max(2, 7/5*max(sv.se$freq.max)), height = 3)
print(p)     # Plot 1 --> in the first page of PDF
dev.off()



# ---- Unormalised
tbl = table(sv.se[(sv.se$len > len.min), c('len.gr', 'freq.max')])
tbl = tbl[rowSums(tbl) != 0,,drop=F]
df.abs = reshape2::melt(tbl)
p <- ggplot(data=df.abs, aes(x=freq.max, y = value, fill=len.gr)) +
  geom_bar(stat="identity") + theme_minimal()  + xlab('Frequency of presence') + ylab('Proportion') + 
  scale_fill_manual(values=color.len,  name ='sSVs len (bp)') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))  +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.direction="vertical",
        legend.box="horizontal",
        legend.box.just = c("top"), 
        legend.background = element_rect(fill=alpha('white', 0.75)),
        legend.margin = margin(2, 2, 2, 2))
p
# 

pdf(paste0(path.figures, 'sv_freq_hist_length_minlen', len.min ,'_abs.pdf'), 
    width = max(2, 7/5*max(sv.se$freq.max)), height = 3)
print(p)     # Plot 1 --> in the first page of PDF
dev.off()

df$value.abs = df.abs$value
write.table(df, paste0(path.figures, 'sv_freq_hist_length_minlen', len.min ,'.txt'),sep = '\t', row.names = F)


