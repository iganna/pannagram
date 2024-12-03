# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
  # library(muscle) #BiocManager::install("muscle")
  library(pannagram)
  library(igraph)
  # library(Biostrings)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func.R", package = "pannagram"))
source(system.file("pangen/comb_func_extra2.R", package = "pannagram"))
source(system.file("pangen/synteny_func.R", package = "pannagram"))
# source("synteny_funcs.R")

# pokazStage('Step 10. Prepare sequences for MAFFT')

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.cons"),        type = "character", default = NULL, help = "path to consensus directory"),
  make_option(c("--path.chromosomes"), type = "character", default = NULL, help = "path to directory with chromosomes"),
  make_option(c("--path.extra"),       type = "character", default = NULL, help = "path to directory, where to combine fasta files for mafft runs"),
  make_option(c("-c", "--cores"),      type = "integer",   default = 1,    help = "number of cores to use for parallel processing"),
  make_option(c("--path.log"),         type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"),        type = "character", default = NULL, help = "Level of log to be shown on the screen"),
  make_option(c("--len.cutoff"),       type = "integer",   default = Inf, help = "Max break considered")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

# ***********************************************************************

aln.type.in = aln.type.add

# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores = opt$cores
if(is.null(num.cores)) stop('Whong number of cores: NULL')

# Max length of break to consider
len.cutoff = opt$len.cutoff

# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop('Consensus folder does not exist')

# Path to chromosomes
if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes

# Path to mafft input
if (!is.null(opt$path.extra)) path.extra <- opt$path.extra

path.work = paste0(path.extra, 'tmp/')
if (!dir.exists(path.work)) {
  dir.create(path.work)
}

# **************************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type.in, ".*")
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.in, "", files)
pref.combinations <- sub(".h5", "", pref.combinations)

if(length(pref.combinations) == 0) {
  stop('No files with the ref-based alignments are found')
}

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

echo = T
for(s.comb in pref.combinations){
  
  if(echo) pokaz('* Combination', s.comb)
  q.chr = strsplit(s.comb, '_')[[1]][1]
  
  file.comb = paste0(path.cons, aln.type.in, s.comb,'.h5')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)
  pokaz("Number of accessions", length(accessions))
  
  save(list = ls(), file =paste0(path.extra, "tmp_breaks.RData"))
  
  # Get breaks
  breaks.init = c()
  for(acc in accessions){
    pokaz(acc)
    s.acc = paste0(gr.accs.e, acc)
    v = h5read(file.comb, s.acc)
    v[is.na(v)] = 0
    breaks.acc = findBreaks(v)
    idx.bad.orientation = which(breaks.acc$val.beg > breaks.acc$val.end)
    if(length(idx.bad.orientation) > 0){
      pokazAttention("Wrong orientation. Number:", length(idx.bad.orientation))
      breaks.acc = breaks.acc[-idx.bad.orientation,,drop=F]
    }
    breaks.init = rbind(breaks.init, breaks.acc)
  }
  
  if(nrow(breaks.init) == 0){
    pokaz('Nothing to add')
    next
  } 
  
  idx.extra = which(breaks.init$len.comb > len.cutoff)
  if(length(idx.extra) > 0){
    breaks.init.extra = breaks.init[idx.extra,]
    breaks.init = breaks.init[-idx.extra,] 
  }
  
  if(nrow(breaks.init) == 0){
    pokaz('Nothing to add')
    next
  }
  
  # Sort and IDs
  breaks.init = breaks.init[order(-breaks.init$idx.end),]
  breaks.init = breaks.init[order(breaks.init$idx.beg),]
  breaks.init$id = 1:nrow(breaks.init)
  
  # Breaks from the previous data
  
  breaks <- mergeOverlaps(breaks.init)
  n.digits <- nchar(as.character(nrow(breaks)))
  format.digits <- paste0("%0", n.digits, "d")
  breaks$gr = 1:nrow(breaks)
  
  # Assign initial breaks to groups
  breaks.init$gr = 0
  for (i.b in 1:nrow(breaks)) {
    pos.b <- breaks$idx.beg[i.b]
    pos.e <- breaks$idx.end[i.b]
    breaks.init[(breaks.init$idx.beg >= pos.b) & (breaks.init$idx.end <= pos.e),]$gr = i.b
  }
  
  if(sum(breaks.init$gr == 0) > 0) stop('Some groups were wrongly defined')
  
  
  # ---- Get the consensus sequences and save then into files ----
  pokaz('Get consensus sequences')
  breaks$id.s = sapply(1:nrow(breaks), function(i.b) paste0('break_',s.comb, '_', sprintf(format.digits, i.b)))
  breaks.init$seq = ''
  
  for(acc in accessions){
    pokaz("Accession", acc)
    # Read the chromosome
    file.chromosome = paste(path.chromosomes, 
                            acc, 
                            '_chr', q.chr, '.fasta', sep = '')
    genome = readFasta(file.chromosome)
    genome = seq2nt(genome)
    
    # Read the alignment
    s.acc = paste0(gr.accs.e, acc)
    v = h5read(file.comb, s.acc)
    v[is.na(v)] = 0
    
    pokaz('Save common breaks..')
    for(i.b in 1:nrow(breaks)){
      
      if((breaks$idx.end[i.b] - breaks$idx.beg[i.b]) == 1) next
      
      pos.b <- breaks$idx.beg[i.b] + 1
      pos.e <- breaks$idx.end[i.b] - 1
      
      v.b = v[pos.b:pos.e]
      s.b = rep('-', length(v.b))
      s.b[v.b != 0] = genome[abs(v.b)]
      idx.rc = which(v.b < 0)
      if(length(idx.rc) > 0){
        s.b[idx.rc] = justCompl(s.b[idx.rc])
      }
      
      s.b = nt2seq(s.b)
      s.b.name = paste0('acc_', acc, '_', breaks$id.s[i.b])
      s.b = setNames(s.b, s.b.name)
      
      file.br.fasta = paste0(path.extra, breaks$id.s[i.b], '_group.fasta')
      if(!file.exists(file.br.fasta)){
        writeFasta(s.b, file.br.fasta) 
      } else {
        writeFasta(s.b, file.br.fasta, append=T)
      }
      
      # TODO
      file.br.idx = paste0(path.extra, breaks$id.s[i.b], '_group.txt')
      line <- paste(c(s.b.name, v.b), collapse = "\t")
      
      con <- file(file.br.idx, open = if (file.exists(file.br.idx)) "a" else "w") # Open the file in append mode if needed
      writeLines(line, con = con, sep = "\n", useBytes = TRUE)           # Write the lines
      close(con)                                  # Close the connection
      
    }
    
    pokaz('Save additional sequences..')
    for(i.b in which(breaks.init$acc == acc)){
      
      # if((breaks.init$idx.end[i.b] - breaks.init$idx.beg[i.b]) == 1) next
      
      pos.b <- breaks.init$val.beg[i.b] + 1
      pos.e <- breaks.init$val.end[i.b] - 1
      
      poses = pos.b:pos.e
      s.b = genome[abs(poses)]
      idx.rc = which(poses < 0)
      if(length(idx.rc) > 0){
        s.b[idx.rc] = justCompl(s.b[idx.rc])
      }
      s.b = nt2seq(s.b)
      s.b.name = paste0('acc_', acc, '_', breaks$id.s[i.b], '|', pos.b, '|', pos.e)
      s.b = setNames(s.b, s.b.name)
      
      breaks.init$seq[i.b] = s.b
      
      # Check length
      if(nchar(breaks.init$seq[i.b]) != breaks.init$len.acc[i.b]){
        save(list = ls(), file = "tmp_wrong.RData")
        stop('Wrong sequence is taken')
      }
      
      file.br.fasta = paste0(path.extra, breaks$id.s[breaks.init$gr[i.b]], '_add.fasta')
      if(!file.exists(file.br.fasta)){
        writeFasta(s.b, file.br.fasta) 
      } else {
        writeFasta(s.b, file.br.fasta, append=T)
      }
      
    }
  }
  
  pokaz('Save..')
  
  if(sum(breaks.init$seq == '') > 0) stop('Some sequences are empty')
  
  file.breaks.info = paste0(path.extra, "breaks_info_",s.comb,".RData")
  # pokaz(file.breaks.info)
  save(list = c("breaks.init", "breaks"), file =file.breaks.info)
  
  gc()
  
  H5close()
  gc()
}


warnings()

