suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
  library(muscle)
  library(pannagram)
})

# source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.features.msa", type = "character", default = NULL, help = "Path to msa directory (features)"),
  make_option("--path.inter.msa",    type = "character", default = NULL, help = "Path to msa directory (internal)"),
  make_option("--path.chromosomes",  type = "character", default = NULL, help = "Path to directory with chromosomes"),
  make_option("--path.mafft.in",     type = "character", default = NULL, help = "Path to directory where to combine fasta files for mafft runs"),
  
  make_option("--max.len.gap",       type = "integer",   default = NULL, help = "Max length of the gap"),
  
  make_option("--cores",             type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",          type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",         type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser, args = args)

#TODO: SHOULD BE PARAMATERS
# len.large = 40000
len.short = 50
n.flank = 30
len.large.mafft = 15000

s.flank.beg = nt2seq(rep('A', n.flank))
s.flank.end = nt2seq(rep('T', n.flank))


# ***********************************************************************
# ---- Logging ----
source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----
source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

aln.type.in = aln.type.clean
aln.type.in = paste0(aln.type.in, '_')

# ***********************************************************************
# ---- Values of parameters ----

# Max len gap
if (is.null(opt$max.len.gap)) {
  stop("Error: max.len.gap is NULL")
} else {
  len.large <- opt$max.len.gap
}

# Number of cores for parallel processing
num.cores = opt$cores
if(is.null(num.cores)) stop('Wrong number of cores: NULL')
pokaz('Number of cores', num.cores, file=file.log.main, echo=echo.main)
if(num.cores > 1){
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
}

# Path with the MSA output (features)
path.features.msa <- opt$path.features.msa
path.inter.msa <- opt$path.inter.msa

if (is.null(path.features.msa) || is.null(path.inter.msa)) {
  stop("Error: both --path.features.msa and --path.inter.msa must be provided")
}

if (!dir.exists(path.features.msa)) stop('Features MSA directory doesn???t exist')
if (!dir.exists(path.inter.msa)) stop('Internal MSA directory doesn???t exist')

# Path to chromosomes
if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes

# Path to mafft input
if (!is.null(opt$path.mafft.in)) path.mafft.in <- opt$path.mafft.in

# ***************************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type.in, ".*")
files <- list.files(path = path.features.msa, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.in, "", files)
pref.combinations <- sub(".h5", "", pref.combinations)

if(length(pref.combinations) == 0) {
  stop('No files with the ref-based alignments are found')
}

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)


# ***************************************************************************
# ---- Functions ----
getSeqsByIdx <- function(idx, v.beg, v.end, acc, genome) {
  
  df <- data.frame(
    idx = idx,
    p1 = v.beg[idx, acc],
    p2 = v.end[idx, acc]
  )
  
  if (any(df$p1 * df$p2 < 0)) {
    stop("Error: beg and end positions have different sign")
  }
  
  df$strand <- c("+", "-", "")[
    (df$p1 == 0 | df$p2 == 0) * 2 +
      (df$p1 < 0 | df$p2 < 0) + 1
  ]
  
  df$p1.own <- 0
  df$p2.own <- 0
  
  ## Positive strand
  idx.pos <- df$strand == "+"
  df$p1[idx.pos] <- df$p1[idx.pos] + 1
  df$p2[idx.pos] <- df$p2[idx.pos] - 1
  df$p1.own[idx.pos] <- df$p1[idx.pos]
  df$p2.own[idx.pos] <- df$p2[idx.pos]
  
  ## Negative strand
  idx.neg <- df$strand == "-"
  tmp <- df$p1[idx.neg]
  df$p1[idx.neg] <- -df$p2[idx.neg] + 1
  df$p2[idx.neg] <- -tmp - 1
  df$p1.own[idx.neg] <- -df$p2[idx.neg]
  df$p2.own[idx.neg] <- -df$p1[idx.neg]
  
  df$len <- df$p2 - df$p1 + 1
  df$len[df$strand == ""] <- 0
  
  check.table <- table(df$strand, df$len != 0)
  if (sum(check.table[, 1] * check.table[, 2]) != 0) {
    stop("Check table did not pass")
  }
  
  seqs <- mapply(
    function(a, b) substr(genome, a, b),
    df$p1,
    df$p2,
    SIMPLIFY = TRUE
  )
  seqs = unlist(seqs)
  names(seqs) <- NULL
  
  if (sum(seqs[df$strand == ""] != "") > 0) {
    stop("Check empty strings did not pass")
  }
  
  ## Reverse complement
  seqs[idx.neg] <- sapply(seqs[idx.neg], revCompl)
  
  list(df = df, seqs = seqs)
}


# ***********************************************************************
# ---- MAIN program body ----

for(s.comb in pref.combinations){
  
  # Log files
  file.log.loop = paste0(path.log, 'test_loop_', s.comb, '.log')
  if(!file.exists(file.log.loop)) invisible(file.create(file.log.loop))
  
  # Check log Done
  if(checkDone(file.log.loop)) next
  
  pokaz('* Combination', s.comb, file=file.log.loop, echo=echo.loop)
  q.chr = strsplit(s.comb, '_')[[1]][1]
  
  # Load variables from the previous step
  file.ws = paste0(path.inter.msa, 'breaks_ws_', s.comb, '.RData')
  load(file.ws)
  
  # Define breaks type
  breaks$type = ''
  breaks$type[(breaks$single == 1) & ((breaks$idx.end - breaks$idx.beg - 1) == 0)] = 'single'
  breaks$type[(breaks$single != 1) & (breaks$len.acc <= len.short)] = 'short'
  breaks$type[(breaks$single != 1) & (breaks$len.acc > len.short) & (breaks$len.mean <= len.large.mafft)] = 'long'
  breaks$type[(breaks$single != 1) & (breaks$len.mean > len.large.mafft)] = 'extra'
  
  # TODO: weak place, there sould not be any breaks$type == ''
  breaks$type[breaks$type == ''] = '-'
  if(sum(breaks$type == '') != 0) stop('Not all of the breaks are classifyed')
  
  # Define indexes for short and singletons
  idx.singl = which(breaks$type == 'single')
  idx.short = which(breaks$type == 'short')
  
  # Define indexes for long sequences
  idx.large = which(breaks$type == 'long')
  idx.extra = which(breaks$type == 'extra')
  
  pokaz('Number of singl/short/large/extra',
        length(idx.singl), 
        length(idx.short), 
        length(idx.large), 
        length(idx.extra), file=file.log.loop, echo=echo.loop)
  
  # ----
  
  # IMPORTANT: THERE ARE SOME BREAKS WITH LOOK LIKE SINGLETONS < BUT THEY ARE NOT
  # if(sum(length(idx.singl) +
  #        length(idx.short) +
  #        length(idx.large) +
  #        length(idx.extra)) != nrow(breaks)) {
  #   save(list = ls(), file = "tmp_wrong_Chrckpoint7.RData")
  #   stop('Chrckpoint7')
  # } 
  
  # Save breaks
  # file.breaks.merged = paste0(path.inter.msa, 'breaks_merged_', s.comb,'.rds')
  # saveRDS(breaks, file.breaks.merged)
  
  ## ---- Save singletons ----
  saveRDS(list(pos.beg = v.beg[idx.singl,],
               pos.end = v.end[idx.singl,],
               ref.pos = data.frame(beg = breaks$idx.beg[idx.singl],
                                    end = breaks$idx.end[idx.singl]) ), 
          paste0(path.inter.msa, 'singletons_',s.comb,'.rds'), compress = F)
  
  ## ---- Save SHORT LONG EXTRA ----
  idx.list <- list(
    short = idx.short,
    large = idx.large,
    extra = idx.extra
  )
  
  for(acc in accessions){
    pokaz(acc)
    pokaz(acc, file=file.log.loop, echo=echo.loop)
    file.chromosome = paste(path.chromosomes, 
                            acc, 
                            '_chr', q.chr, '.fasta', sep = '')
    if(!file.exists(file.chromosome)){
      pokaz('File', file.chromosome, 'does not exist')
      stop()
    }
    genome = readFasta(file.chromosome)
    
    for (type in names(idx.list)) {
      pokaz(type)
      
      idx.save <- idx.list[[type]]
      
      if(length(idx.save) == 0) next
      
      res <- getSeqsByIdx(idx.save, v.beg, v.end, acc, genome)
      df.save   <- res$df
      seqs.save <- res$seqs
      if (any(df.save$len != nchar(seqs.save))) stop("No match in sequence length")
      
      if (type != "short") {
        seqs.save <- paste0(s.flank.beg, seqs.save, s.flank.end)
        seqs.save[df.save$len == 0] <- ""
      }
      
      file.seqs <- paste0(path.inter.msa, acc, "_", type, "_", s.comb, ".txt")
      file.df   <- paste0(path.inter.msa, acc, "_", type, "_", s.comb, "_df.rds")
      
      
      save(list = c("res", "acc", "seqs.save", "file.seqs"), file = "tmp_workspace.RData")
      
      writeLines(c(acc, seqs.save), file.seqs)
      saveRDS(df.save, file.df)
      
      rm(res, df.save, seqs.save)
    }
    
    rmSafe(genome)
  }
  
  pokaz('Done.', file=file.log.loop, echo=echo.loop)
}

if(num.cores > 1){
  stopCluster(myCluster)
}

warnings()

stop('TEST Script has been finished')