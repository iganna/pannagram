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
sim.cutoff = 0.9


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
# ---- Reading the data ----

file.sv.pos = paste0(path.sv, 'sv_pangen_pos.rds')
if(!file.exists(file.sv.pos)){
  stop('SVs were not generated.')
}
sv.all = readRDS(file.sv.pos)
sv.all$chr = as.numeric(sv.all$chr)

sv.se = sv.all[sv.all$single == 1,]
sv.se$len.gr =  cut(sv.se$len, breaks = len.bins, right = FALSE, labels = len.labels)

f.max = max(sv.se$freq.max)

res.cover.file = 'seq_sv_big_on_sv_cover.rds'
res.nest = readRDS(paste(path.sv, res.cover.file, sep = ''))

# ***********************************************************************
# ----Construct the graph ----

file.g.content = paste0(path.sv, 'g_content_sim',round(sim.cutoff * 100),'.rds')
# remove duplicates in orientations
if(!file.exists(file.g.content)){
  idx = which(duplicated(res.nest[,c('V1', 'V8')]))
  if(length(idx) > 0){
    
    res.nest$comb = paste0(res.nest$V1, '|', res.nest$V8)
    res.nest$id = 1:nrow(res.nest)
    res.check = res.nest[res.nest$comb %in% res.nest$comb[idx],]
    res.check = res.check[order(res.check$comb),]
    cnt = table(res.check$comb)
    if(sum(cnt != 2) > 0) stop('Somethig is wrong with combinations')
    
    res.check$score = res.check$C1 + res.check$C8
    score.higher = res.check$score[seq(1,nrow(res.check), 2)] > res.check$score[seq(2,nrow(res.check), 2)]
    
    idx.remove = c(res.check$id[seq(1,nrow(res.check), 2)][!score.higher],
                   res.check$id[seq(2,nrow(res.check), 2)][score.higher])
    
    if(length(idx.remove) > 0){
      res.nest = res.nest[-idx.remove,,drop=F]
    }
  }
  
  g.content = getGraphFromBlast(res.nest = res.nest, sim.cutoff = sim.cutoff, collapse = T)
  
  saveRDS(g.content, file.g.content)
} else {
  g.content = readRDS(file.g.content)
}

# ***********************************************************************
# ----Construct the graph ----

