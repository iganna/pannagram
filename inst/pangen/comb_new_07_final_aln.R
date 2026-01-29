# Combine all alignments together into the final one

suppressMessages({
  library(rhdf5)
  library(foreach)
  library(doParallel)
  library(optparse)
})

source(system.file("utils/utils.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.mafft.out",   type = "character", default = NULL, help = "Path to directory where mafft results are"),
  make_option("--path.features.msa",type = "character", default = NULL, help = "Path to msa directory (features)"),
  make_option("--path.inter.msa",   type = "character", default = NULL, help = "Path to msa directory (internal)"),
  make_option("--path.inter.synteny",type = "character", default = NULL, help = "Path to alignments the between synteny blocks"),
  make_option("--accessions",        type = "character", default = NULL, help = "File containing accessions to analyze"),
  
  make_option("--cores",            type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",         type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",        type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser, args = args)

n.flank = 30

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

aln.type.in = paste0(aln.type.clean, '_')
aln.type.out = paste0(aln.type.msa, '_')

# ***********************************************************************
# ---- Accessions ----

file.acc <- ifelse(!is.null(opt$accessions), opt$accessions, stop("File with accessions are not specified"))
accessions.specified <- as.character(read.table(file.acc, stringsAsFactors = FALSE)[, 1])

# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores <- opt$cores

if (!is.null(opt$path.mafft.out)) path.mafft.out <- opt$path.mafft.out

# Path with the MSA output (features)
path.features.msa <- opt$path.features.msa
path.inter.msa <- opt$path.inter.msa
path.inter.synteny <- opt$path.inter.synteny

if (is.null(path.features.msa) || is.null(path.inter.msa)) {
  stop("Error: both --path.features.msa and --path.inter.msa must be provided")
}

if (!dir.exists(path.features.msa)) stop('Features MSA directory does not exist')
if (!dir.exists(path.inter.msa)) stop('Internal MSA directory does not exist')

# ***********************************************************************

# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type.in, ".*")
files <- list.files(path = path.features.msa, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.in, "", files)
pref.combinations <- sub(".h5", "", pref.combinations)

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)


# ***********************************************************************
# ---- MAIN program body ----

stat.comb <- data.frame(comb = character(),
                        coverage = numeric(),
                        stringsAsFactors = FALSE)

pref.combinations = '5_5'
for(s.comb in pref.combinations){
  
  # Log files
  file.log.loop = paste0(path.log, 'loop_', s.comb, '.log')
  if(!file.exists(file.log.loop)) invisible(file.create(file.log.loop))
  
  # Check log Done
  if(checkDone(file.log.loop)) next
  
  pokaz('* Combination', s.comb, file=file.log.main, echo=echo.main)
  
  # ---- PRE-Resultant File ----
  
  file.res.pre = paste0(path.inter.msa, aln.type.out, s.comb,'.h5')
  # if (file.exists(file.res.pre)) file.remove(file.res.pre)
  if (!file.exists(file.res.pre)){
    h5createFile(file.res.pre)
    suppressMessages({
      h5createGroup(file.res.pre, gr.accs.e)
    })  
  } else {
    pokaz('Preliminary file exists')
  }
  
  
  # Paths
  path.short.aln = paste0(path.inter.synteny, 'short_',s.comb,'/')
  path.large.aln = paste0(path.inter.synteny, 'large_',s.comb,'/')

  # Get accessions
  file.comb = paste0(path.features.msa, aln.type.in, s.comb,'.h5')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  accessions = intersect(accessions, accessions.specified)
  n.acc = length(accessions)
  if (n.acc == 0) stop(paste("No accessions for combination", s.comb))
  
  len.aln.synteny = length(h5read(file.comb, paste0(gr.accs.e, accessions[1])))
  
  # ---- Read breaks ----
  # The most important is to get estimates of len.new!
  # Then read everything by accession and insert corresponding positions
  df.breaks = readRDS(paste0(path.inter.msa, "breaks_annotated_", s.comb,".rds"))
  df.breaks = df.breaks[,c("idx.beg", "idx.end", "type")]
  df.breaks$len.new = 0
  
  # ---- Read singletons results ----
  data.single = readRDS(paste0(path.inter.msa, "singletons_", s.comb, ".rds"))
  df.single = data.single$ref.pos
  
  idx.beg <- which(data.single$pos.beg != 0, arr.ind = TRUE)
  idx.end <- which(data.single$pos.end != 0, arr.ind = TRUE)
  if(sum(idx.beg != idx.end) != 0) stop('Wrong data for singletons')
  idx.beg = idx.beg[order(idx.beg[,1]),]
  if(nrow(idx.beg) != idx.beg[nrow(idx.beg),1]) stop('Indexes in singleton are wrong')
  
  df.single <- cbind(df.single, data.frame(
    acc.beg   = data.single$pos.beg[idx.beg],
    acc.end   = data.single$pos.end[idx.beg],
    acc = colnames(data.single$pos.beg)[idx.beg[, "col"]]
  ))
  df.single$extra.pos = df.single$acc.end - df.single$acc.beg - 1  # -1 is important!!!
  
  df.breaks$len.new[df.breaks$type == 'single'] = df.single$extra.pos
  
  # ---- Read len.new from short results ----
  
  data.short = readLines(paste0(path.short.aln, accessions[1], "_short_", s.comb, ".txt.aln.txt"))
  
  # REMOVE LATER:
  if((length(data.short) - 1) == sum(df.breaks$type == 'short')){
    data.short = data.short[-1]
  }
  
  if(length(data.short) != sum(df.breaks$type == 'short')) stop('Short alignments do not match')
  df.breaks$len.new[df.breaks$type == 'short'] = nchar(data.short)
  
  rm(data.short)
    
  # ---- Read len.new the Mafft results ----
  
  data.large = readLines(paste0(path.large.aln, accessions[1], "_large_", s.comb, ".txt.txt"))
  
  # REMOVE LATER:
  if((length(data.large) - 1) == sum(df.breaks$type == 'long')){
    data.large = data.large[-1]
  }
  
  if(length(data.large) != sum(df.breaks$type == 'long')) stop('Large alignments do not match')
  df.breaks$len.new[df.breaks$type == 'long'] = nchar(data.large)
  
  rm(data.large)
  
  save(list = ls(), file = "tmp_workspace.RData")
  
  # # TODO
  # # Checkup
  # tmp = table(df.breaks$type, df.breaks$len.new != 0)
  # if (any(tmp[,1] * tmp[,2]) != 0) stop('Wrong types and len.new')
  df.breaks = df.breaks[df.breaks$type != 'extra',]
  df.breaks = df.breaks[df.breaks$type != '-',]
  
  # ---- New combined coordinates ----
  
  # Remove those, who do not change the len.new: extra alignments
  # df.breaks = df.breaks[df.breaks$len.new > 0,]  # Put back this line, when the previous TODO is managed
  df.breaks$len = df.breaks$idx.end - df.breaks$idx.beg - 1
  df.breaks$extra = df.breaks$len.new - df.breaks$len
  
  # TODO
  df.breaks$fail = (df.breaks$len.new == 0) | (df.breaks$extra < 0)
  df.breaks$extra[df.breaks$fail] = 0
  
  if(any(df.breaks$extra < 0)){
    stop('Why negative extra?')
  } 
  
  # Mapping old coordinates to new
  idx.extra = rep(0, len.aln.synteny)
  idx.extra[df.breaks$idx.beg + 1] = df.breaks$extra
  idx.extra = cumsum(idx.extra)
  idx.map = (1:len.aln.synteny) + idx.extra
  
  if(any(idx.map <= 0)) stop('Wrong mapping: zeros or negative')
  
  df.breaks$new.beg = idx.map[df.breaks$idx.beg] + 1
  df.breaks$new.end = idx.map[df.breaks$idx.end] - 1
  df.breaks$len.new[df.breaks$fail] = df.breaks$new.end[df.breaks$fail] - df.breaks$new.beg[df.breaks$fail] + 1
  
  if(any(df.breaks$new.beg > df.breaks$new.end)) stop('Beging is highre than end')
  
  tmp = df.breaks$new.end - df.breaks$new.beg + 1
  if(any(df.breaks$len.new != tmp)) stop('Wrong mapping: length mismatch')
  
  # Length of new alignment
  len.aln.new = max(idx.map)
  
  if(len.aln.new != sum(df.breaks$extra) + len.aln.synteny) stop('Wrong mapping: pangenome length')
  
  # Idx which should be zero after mapping
  idx.zero = rep(0, len.aln.new)
  idx.zero[df.breaks$new.beg] = 1
  idx.zero[df.breaks$new.end] = idx.zero[df.breaks$new.end] - 1
  idx.zero = cumsum(idx.zero)
  idx.zero[df.breaks$new.beg] = 1
  idx.zero[df.breaks$new.end] = 1
  
  if(sum(idx.zero) != sum(df.breaks$len.new)) stop('Idx.zero are wrongly defined')
  
  # ---- Get results by accessions ----
  idx.all.acc.zeros = rep(0, len.aln.new)
  for(acc in accessions){
    
    # Log files
    file.log.loop.acc = paste0(path.log, 'loop_', s.comb, '_', acc, '.log')
    if(!file.exists(file.log.loop.acc)) invisible(file.create(file.log.loop.acc))
    
    # Check log Done
    if(checkDone(file.log.loop.acc)){
      pokaz('Accession', acc ,'was analysed before')
      v.new = h5read(file.res.pre, paste0(gr.accs.e, acc))
      idx.all.acc.zeros = idx.all.acc.zeros + (v.new == 0)
      next
    }
    
    pokaz('Accession', acc)
    
    # Read accession alignment
    v = h5read(file.comb, paste0(gr.accs.e, acc))
    
    # Map new alignment
    v.new = rep(0, len.aln.new)
    v.new[idx.map] = v
    v.new[idx.zero == 1] = 0
    
    # Fill up singletons
    df.br.tmp = df.breaks[df.breaks$type == 'single',]
    for(irow in which(df.single$acc == acc)){
      # stop()
      v.new[df.br.tmp$new.beg[irow]:df.br.tmp$new.end[irow]] = (df.single$acc.beg[irow]+1):(df.single$acc.end[irow]-1)
      
      # tmp = c(df.br.tmp$new.beg[irow] - 1,
      #         df.br.tmp$new.beg[irow],
      #         df.br.tmp$new.end[irow],
      #         df.br.tmp$new.end[irow] + 1)
      # pokaz(v.new[tmp])
    }
    
    # Check duplicates
    if(sum(duplicated(abs(v.new[v.new != 0]))) > 0) stop('Duplicated after singletons')
    
    # Fill up short alignments
    for(s.type in c('short', 'large')){
      pokaz('Insert type', s.type)
      
      if(s.type == 'large'){
        df.br.tmp = df.breaks[df.breaks$type == 'long',]  
      } else {
        df.br.tmp = df.breaks[df.breaks$type == s.type,]
      }
      
      file.df.acc = paste0(path.inter.synteny, acc, "_",s.type,"_", s.comb, "_df.rds")
      
      if(s.type == 'short'){
        s.ext = ".txt.aln.txt"
      } else {
        s.ext = ".txt.txt"
      }
      
      file.aln.acc = paste0(path.inter.synteny, s.type, "_", s.comb, '/' , 
                            acc, "_",s.type,"_", s.comb, s.ext)
      
      df.acc = readRDS(file.df.acc)
      aln.acc = readLines(file.aln.acc)
      
      # REMOVE LATER:
      if((length(aln.acc) - 1) == nrow(df.br.tmp)){
        aln.acc = aln.acc[-1]
      }
      if (length(unique(c(nrow(df.br.tmp),
                          length(aln.acc),
                          nrow(df.acc)))) != 1) {
        stop('Alignment wrong lengths')
      }
      
      # if(s.type == 'large'){
      #   save(list = ls(), file = "tmp_workspace_large.RData")
      # }
      
      for(i in 1:length(aln.acc)){
        # pokaz(i)
        
        if(df.br.tmp$fail[i]) next  # Kostyl
        if(df.br.tmp$extra[i] < 0) next
        
        # if(s.type == 'large') stop()
        if(df.acc$p1[i] == 0) next
        
        p.own = df.acc$p1.own[i]:df.acc$p2.own[i]
        
        p.insert = rep(0, df.br.tmp$len.new[i])
        idx.tmp.aln = c(gregexpr("[^-]", aln.acc[i])[[1]])
        
        if(length(idx.tmp.aln) != length(p.own)) {
  
          if((length(idx.tmp.aln) - 60) != length(p.own)){
            pokaz(i, length(idx.tmp.aln), length(p.own))
            next  # Kostyl
          } 
          
          idx.tmp.aln <- idx.tmp.aln[(n.flank + 1):(length(idx.tmp.aln) - n.flank)]
        }
        
        if(length(idx.tmp.aln) != length(p.own)) stop('Wrong length aligned')
        p.insert[idx.tmp.aln] = p.own
        
        v.new[(df.br.tmp$new.beg[i]):(df.br.tmp$new.end[i])] = p.insert
        
        # tmp = (df.br.tmp$new.beg[i] - 1) : (df.br.tmp$new.end[i] + 1)
        # pokaz(v.new[tmp])
      }
      
      # Check duplicates # Kostyl
      if(sum(duplicated(abs(v.new[v.new != 0]))) > 0){
        pokaz('Duplicated in', s.type, sum(duplicated(abs(v.new[v.new != 0]))))
        dup.values = abs(v.new[duplicated(abs(v.new))])
        v.new[abs(v.new) %in% dup.values] = 0
        pokaz('Duplicated after', sum(duplicated(abs(v.new[v.new != 0]))))
      } 
    }
    
    # Save positions which are zeros
    idx.all.acc.zeros = idx.all.acc.zeros + (v.new == 0)
    
    # Save
    suppressMessages({
      h5write(v.new, file.res.pre, paste0(gr.accs.e, acc)) })
      
    pokaz('Done.', file=file.log.loop.acc, echo=echo.loop)
  }
  
  idx.remain = (idx.all.acc.zeros != length(accessions))
  
  # Remove zeros
  pokaz('Resultant file: remove seros')
  for(acc in accessions){
    
    # ---- Resultant File ----
    file.res = paste0(path.features.msa, aln.type.out, s.comb,'.h5')
    # if (file.exists(file.res)) file.remove(file.res)
    if (!file.exists(file.res)) {
      h5createFile(file.res)
      suppressMessages({
        h5createGroup(file.res, gr.accs.e)
      })  
    }
    
    pokaz('Accession', acc)
    
    # Read accession alignment
    v = h5read(file.res.pre, paste0(gr.accs.e, acc))
    v.remain = v[idx.remain]
    
    # Save
    suppressMessages({
      h5write(v.remain, file.res, paste0(gr.accs.e, acc)) })
    
  }

  H5close()
  gc()
  
  pokaz('Done.', file=file.log.loop, echo=echo.loop)
  
}

warnings()

pokaz('Done.', file=file.log.main, echo=echo.main)

