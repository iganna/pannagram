
### imports
source("visualisation/visualisation.R")
source("pangen/synteny_func_plot.R")
source("utils/utils.R")

suppressMessages({
  require("optparse")
  library("parallel")
})

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
    make_option(c("--ref"),
        type = "character",
        default = NULL,
        help = "Name af reference sequence",
        metavar = "character"
    ),
    make_option(c("--path_ref"),
        type = "character",
        default = NULL,
        help = "Path to reference sequence",
        metavar = "character"
    ),
    make_option(c("--path_out"),
        type = "character",
        default = NULL,
        help = "Path to the output directory",
        metavar = "character"
    ),
    make_option(c("--path_chr"),
                type = "character",
                default = NULL,
                help = "Path to chromosomes",
                metavar = "character"
    ),
    make_option(c("--algn_path"),
        type = "character",
        default = NULL,
        help = "Path to the alignment directory",
        metavar = "character"),
    make_option(c("--path.log"), 
        type = "character", 
        default = NULL,
        help = "Path for log files",
        metavar = "character"),
    make_option(c("--log.level"), 
        type = "character", 
        default = NULL,
        help = "Level of log to be shown on the screen", 
        metavar = "character"),
    make_option(c("--cores"),
                type = "integer",
                default = 1,
                help = "Number of cores to use",
                metavar = "character"
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

# Logging
source('utils/chunk_logging.R') # a common code for all R logging

if (
    is.null(opt$ref) ||
    is.null(opt$path_ref) ||
    is.null(opt$path_chr) ||
    is.null(opt$path_out) ||
    is.null(opt$algn_path)
) stop("All mandatory arguments must be provided.")

ref <- opt$ref
path.ref <- opt$path_ref
path.chr <- opt$path_chr
path.out <- opt$path_out
path.aln<- opt$algn_path

# Extracting only ids from path
pattern <- ".*_[0-9]+_[0-9]+_maj\\.rds$"
files.aln <- list.files(path = path.aln, pattern = pattern, full.names = F)
acc.ids <- unique(sapply(files.aln, function(s) sub("^(.*?)_\\d+_\\d+_maj\\.rds$", "\\1", s)))
pokaz('Genomes analysed:', acc.ids, file=file.log.main, echo=echo.main)

# # Find the file with the reference genome
# ext <- c('fasta', 'fna', 'fa', 'fas')
# ref.name <- findGenomeFile(genome.pref = ref, 
#                            path.genome = path.ref,
#                            ext = ext)
# if (is.null(ref.name)) stop('No reference genome files found in the specified folder')

# Output folder with plots
pdf.path <- normalizePath(file.path(path.out, paste0("plots_", ref)), mustWork = FALSE)
dir.create(pdf.path, showWarnings = FALSE, recursive = TRUE)

print(path.aln)

for (acc in acc.ids){
  pokaz('Accession', acc, 
        file=file.log.main, echo=echo.main)
  
  # Lengths of chromosomes for the accession and reference
  chr.len = c()
  for(id in c(ref, acc)){
    file.chr.len = paste0(path.chr, id, '_chr_len.txt')  
    # pokaz(file.chr.len)
    if(file.exists(file.chr.len)){
      chr.len.id = read.table(file.chr.len, header=T)
      chr.len = rbind(chr.len, chr.len.id)  
    } else {
      for(i.chr in 1:1e10){
        file.chr =  paste0(path.chr, id, '_chr',i.chr,'.fasta')
        pokaz('Chromosomel file', file.chr)
        if(!file.exists(file.chr)) break
        pokaz('yes')
        seq.chr = readFastaMy(file.chr)
        chr.len.tmp = data.frame(acc = id,
                                 chr = i.chr,
                                 len = nchar(seq.chr))
        chr.len = rbind(chr.len, chr.len.tmp)
      }
      write.table(chr.len, file.acc.len, sep = '\t', col.names = T, row.names = F, quote = F)
    }
  }
  if (is.null(acc)) stop('No target genome files found in the specified folder')
  
  # Get ggplot with the synteny
  p <- plotSynAllChr(path.aln,
                            acc=acc,
                            ref=ref,
                            chr.len=chr.len)
  
  # Save
  pdf.name <- paste0(ref, "-", id)
  savePDF(p, path = pdf.path, name = pdf.name)
}

pokaz('Done.',
      file=file.log.main, echo=echo.main)

