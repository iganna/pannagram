# Load the necessary library
library(optparse)
source('sim/sim_func.R')
source('utils/utils.R')

# Define options
option_list = list(
  make_option(c("--out"), type = "character", default = NULL,
              help = "Path to the output coverage file", metavar = "FILE"),
  make_option(c("--sim"), type = "numeric", default = 90,
              help = "Similarity threshold", metavar = "NUMBER")
)

# Create the option parser
opt_parser = OptionParser(option_list = option_list)

# Parse the arguments
opt = parse_args(opt_parser)

# Check for the presence of all required arguments
output.file <- ifelse(!is.null(opt$out), opt$out, 
                      stop("Output file not specified", call. = FALSE))
sim.cutoff <- ifelse(!is.null(opt$sim), opt$sim, 
                     stop("Similarity threshold not specified", call. = FALSE))
sim.cutoff = as.numeric(sim.cutoff)


# ---- Testing ----

# blast.file = 'tmp.txt.blast.tmp'
# source('pannagram/utils.R')
# fasta.file = 'new_genes/new_genes.fasta'
# sim.cutoff = 0.85


# ---- Main ----

output.dir = dirname(output.file)

files <- list.files(path = output.dir, pattern = paste0(".*",sim.cutoff,"\\.cnt$"), full.names = T)
if(length(files) == 0) stop('No files with results')
# print(files)

total.cnt.list = list()
total.cnt.names = c()

for (i.file in 1:length(files)) {
  file = files[i.file]
  data <- read.table(file, header = TRUE, stringsAsFactors = F)
  total.cnt.list[[i.file]] = data[,ncol(data), drop=F]
  total.cnt.names = unique(c(total.cnt.names, rownames(data)))
}

mx.cnt = matrix(0, nrow = length(total.cnt.names), ncol = length(files), 
                dimnames = list(total.cnt.names, NULL))

for (i.file in 1:length(files)) {
  mx.cnt[rownames(total.cnt.list[[i.file]]),i.file] = total.cnt.list[[i.file]][,1]
}


# Colnames

acc.names <- sapply(basename(files), function(s){
  s = strsplit(s, '\\.')[[1]]
  return(s[length(s) - 1])
})

colnames(mx.cnt) = acc.names
# pokaz(acc.names)


# Save
write.table(mx.cnt, paste0(output.file, '.total_', sim.cutoff, '.cnt'), 
            quote = F, row.names = T, col.names = T, sep = '\t')






