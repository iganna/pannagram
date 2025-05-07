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
  make_option(c("--words.keep"),     type = "character", default = NULL, help = "Substring in names to keep for sure"),
  
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
words.keep <- opt$words.keep

if(!is.null(words.remain)){
  words.remain = gsub(" ", "", words.remain)
  words.remain = strsplit(words.remain, ',')[[1]]
  pokaz("Words to remain:", words.remain)
}

if(!is.null(words.remove)){
  words.remove = gsub(" ", "", words.remove)
  words.remove = strsplit(words.remove, ',')[[1]]
  pokaz("Words to remove:", words.remove)
}

if(!is.null(words.keep)){
  words.keep = gsub(" ", "", words.keep)
  words.keep = strsplit(words.keep, ',')[[1]]
  pokaz("Words to keep for sure:", words.keep)
}

# Create folders for the alignment results
if(!dir.exists(path.genomes)) dir.create(path.genomes)

query.types <- c('fasta', 'fna', 'fa', 'fas')

# ***********************************************************************
# ---- Correspondence ----

files <- list.files(path.genomes, pattern = paste0("\\.(", paste(query.types, collapse = "|"), ")$"), full.names = F, ignore.case = TRUE)
# print(files)

df.names.removed = c()
df.names.remained = c()
chr.cnt = c()
for(f in files){
  
  pokaz("Analysing genome:", f)
  genome = readFasta(paste0(path.genomes, f))
  names.chr = names(genome)
  names.chr = gsub(" ", '_', names.chr)
  names(genome) = names.chr
  names.chr.init = names.chr
  
  # Filtering
  if(!is.null(words.remain)){
    names.chr.remain = c()
    for(w in words.remain){
      names.chr.remain = c(names.chr.remain, names.chr[grepl(w, names.chr)])
    }
    names.chr.remain = unique(names.chr.remain)
    names.chr.remain = sort(names.chr.remain)
    names.chr = names.chr.remain
  }
  pokaz(names.chr)
  
  if(!is.null(words.remove)){
    for(w in words.remove){
      names.chr = names.chr[!grepl(w, names.chr)]
    }
  }
  pokaz(names.chr)
  pokaz('--------')
  
  if(! is.null(words.keep)){
    names.chr.keep = c()
    for(w in words.keep){
      names.chr.keep = c(names.chr.keep, names.chr.init[grepl(w, names.chr.init)])
      pokaz(w, ':', names.chr.keep)
    }
    names.chr = c(names.chr, names.chr.keep)
    names.chr = unique(names.chr)
    names.chr = sort(names.chr)
  }
  pokaz(names.chr)
  
  # Save
  names.remove = setdiff(names(genome), names.chr)
  if(length(names.remove) > 0){
    df.names.removed = rbind(df.names.removed, 
                             data.frame(genome = f, names = names.remove))  
  }
  
  if(length(names.chr) > 0){
    df.names.remained = rbind(df.names.remained, 
                             data.frame(genome = f, names = names.chr))  
    writeFasta(genome[names.chr], paste0(path.filtered, f))
  } else {
    pokazAttention("Genome", f, "does not have chromosomes remained.")
  }
  
  chr.cnt = rbind(chr.cnt, 
                  data.frame(genome = f, remained = length(names.chr), removed = length(names.remove)))
  
}

if(!is.null(df.names.removed)){
  write.table(df.names.removed, paste0(path.filtered, 'chr_names_removed.txt'), row.names = F, col.names = T, quote = F)  
}
if(!is.null(df.names.remained)){
  write.table(df.names.remained, paste0(path.filtered, 'chr_names_remained.txt'), row.names = F, col.names = T, quote = F)  
}

write.table(chr.cnt, paste0(path.filtered, 'chr_counts.txt'), row.names = F, col.names = T, quote = F)  

print(chr.cnt)


