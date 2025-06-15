# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
  library(ggplot2)
})

library(pannagram)

# Define blocks in the alignemnt

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--file.pi",    type = "character", default = NULL, help = "Path to consensus directory"),
  make_option("--path.figures", type = "character", default = "",   help = "Path to folder with figures"),
  make_option("--path.snp", type = "character", default = NULL, help = "Path to snp dir")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);



# ***********************************************************************
# ---- Values of parameters ----

path.snp <- opt$path.snp
if (!dir.exists(path.snp)) stop('Folder `path.snp` does not exist')

# Path with the consensus output
file.pi <- opt$file.pi
if(!file.exists(file.pi)) stop('File with pi doesn’t exist')
file.pi = basename(file.pi)


# Path with the figures folder
path.figures <- opt$path.figures
if(!dir.exists(path.figures)) stop('Consensus folder doesn’t exist')


file.dist = paste0(path.snp, file.pi, '.dist.dist' )
file.id = paste0(path.snp, sub('.vcf', '_output.ID.FORMAT', file.pi))


lines <- readLines(file.dist)

info = strsplit(lines, '\t')
n = length(info[[length(info)]]) + 1

dist.mx = matrix(0, nrow = n, ncol = n)
for(irow in 1:length(info)){
  info.tmp = info[[irow]]
  dist.mx[irow + 1, 1:length(info.tmp)] = as.numeric(info.tmp)
}

dist.mx = dist.mx + t(dist.mx)

labels = read.table(file.id)
labels = labels[1,-c(1,2)]
colnames(dist.mx) <- labels
rownames(dist.mx) <- labels

h = hclust(as.dist(dist.mx))
plot(h)

dist.mx = dist.mx[h$order,h$order]

p = heatplot(dist.mx)

pokaz('Plot...')
pokaz(paste0(path.figures, file.pi, '_dist.pdf'))
pdf(paste0(path.figures, file.pi, '_dist.pdf'), width = 7, height = 7)
print(p)
dev.off()

saveRDS(dist.mx, paste0(path.snp, file.pi, '_dist.rds'))



