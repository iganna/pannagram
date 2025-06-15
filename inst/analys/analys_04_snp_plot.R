# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
  library(ggplot2)
  library(pannagram)
})

# Define blocks in the alignemnt

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--file.pi",    type = "character", default = NULL, help = "Path to consensus directory"),
  make_option("--path.figures", type = "character", default = "",   help = "Path to folder with figures"),
  make_option("--len.wnd", type = "integer", default = 200000,   help = "Window length for plotting")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

len.wnd <- opt$len.wnd

# Path with the consensus output
file.pi <- opt$file.pi
if(!file.exists(file.pi)) stop('File with pi doesn’t exist')

# Path with the figures folder
path.figures <- opt$path.figures
if(!dir.exists(path.figures)) stop('Consensus folder doesn’t exist')

pi.acc = read.table(file.pi, header = 1)
pi.acc = piVCF(pi.acc, len.wnd)

p.acc = ggplot(pi.acc, aes(x = pos, y = pi)) +
  geom_point(color = "#176B87", size = 1) +
  geom_smooth(method = "loess", 
              color = "#A02334",
              fill = "#A02334",
              se = T,
              span = 0.2) +  
  ylab(expression(pi)) + 
  xlab(NULL) +
  theme_minimal() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

saveRDS(p.acc, paste0(path.figures, basename(file.pi), '_smooth.rds'))

savePDF(p.acc, path = path.figures, name = paste0(basename(file.pi), '_smooth'), width = 6, height = 1)


