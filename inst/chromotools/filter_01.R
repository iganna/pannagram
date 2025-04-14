# Alignment-1. Remaining syntenic (major) matches

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(pannagram)
  library(crayon)
})

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("--path.genomes"),   type = "character", default = NULL, help = "Path to the genome folder"),
  make_option(c("--path.filtered"),  type = "character", default = NULL, help = "Path to the output directory with alignments"),
  make_option(c("--words.remove"),   type = "character", default = NULL, help = "Substring in names to remove"),
  make_option(c("--words.remain"),   type = "character", default = NULL, help = "Substring in names to remain"),
  
  make_option(c("--cores"),         type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option(c("--path.log"),      type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"),     type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Variables ----

max.len = 10^6
len.blast = 50

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- Values of parameters ----

# Number of cores
num.cores <- opt$cores

path.genomes  <- ifelse(!is.null(opt$path.genomes), opt$path.genomes, stop('Folder with Genomes is not specified'))
path.filtered <- ifelse(!is.null(opt$path.filtered), opt$path.filtered, stop('Folder with Genomes is not specified'))

words.remove <- opt$words.remove
words.remain <- opt$words.remain

if(!is.null(words.remain)){
  words.remain = strsplit(words.remain, ',')[[1]]
  pokaz("Words to remain:", words.remain)
}

if(!is.null(words.remove)){
  words.remove = strsplit(words.remove, ',')[[1]]
  pokaz("Words to remove:", words.remove)
}

# Create folders for the alignment results
if(!dir.exists(path.genomes)) dir.create(path.genomes)

query.types <- c('fasta', 'fna', 'fa', 'fas')

# ***********************************************************************
# ---- Correspondence ----

files <- list.files(path.genomes, pattern = paste0("\\.(", paste(query.types, collapse = "|"), ")$"), full.names = F, ignore.case = TRUE)
# print(files)

df.names.removed = c()
for(f in files){
  
  pokaz("Analysing genome:", f)
  genome = readFasta(paste0(path.genomes, f))
  names.chr = names(genome)
  names.chr = gsub(" ", '_', names.chr)
  names(genome) = names.chr
  
  # Filtering
  if(!is.null(words.remain)){
    for(w in words.remain){
      names.chr = names.chr[grepl(w, names.chr)]
    }
  }
  if(!is.null(words.remove)){
    for(w in words.remove){
      names.chr = names.chr[!grepl(w, names.chr)]
    }
  }
  
  # Save
  names.remove = setdiff(names(genome), names.chr)
  if(length(names.remove) > 0){
    df = data.frame(genome = f, names = names.remove)
    df.names.removed = rbind(df.names.removed, df)  
  }
  
  writeFasta(genome[names.chr], paste0(path.filtered, f))
}

if(!is.null(df.names.removed)){
  write.table(df.names.removed, paste0(path.filtered, 'removed_chr_names.txt'), row.names = F, col.names = T, quote = F)  
}




