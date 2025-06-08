# Combine all alignments together into the final one

suppressMessages({
  library(rhdf5)
  library(foreach)
  library(doParallel)
  library(optparse)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func_mafft_refine.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.mafft.out",   type = "character", default = NULL, help = "Path to directory where mafft results are"),
  make_option("--path.features.msa",type = "character", default = NULL, help = "Path to msa directory (features)"),
  make_option("--path.inter.msa",   type = "character", default = NULL, help = "Path to msa directory (internal)"),
  
  make_option("--cores",            type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",         type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",        type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser, args = args)

n.flank = 30
max.block.elemnt = 3 * 10^ 6

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

aln.type.in = aln.type.clean
aln.type.out = aln.type.msa

# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores <- opt$cores

if (!is.null(opt$path.mafft.out)) path.mafft.out <- opt$path.mafft.out

# Path with the MSA output (features)
path.features.msa <- opt$path.features.msa
path.inter.msa <- opt$path.inter.msa

if (is.null(path.features.msa) || is.null(path.inter.msa)) {
  stop("Error: both --path.features.msa and --path.inter.msa must be provided")
}

if (!dir.exists(path.features.msa)) stop('Features MSA directory doesn???t exist')
if (!dir.exists(path.inter.msa)) stop('Internal MSA directory doesn???t exist')

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


for(s.comb in pref.combinations){
  
  pokaz('* Combination', s.comb, file=file.log.main, echo=echo.main)

  # Get accessions
  file.comb = paste0(path.features.msa, aln.type.in, s.comb,'.h5')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)
  base.len = length(h5read(file.comb, paste0(gr.accs.e, accessions[1])))
  
  
  # n.rows.block = round(max.block.elemnt / n.acc)
  # seq.blocks = seq(1, base.len - n.rows.block, n.rows.block)
  
  # ---- All MAFFT results for the combination ----
  pref = paste('Gap', s.comb, sep = '_')
  
  
  files.mafft1 = list.files(path = path.mafft.out, 
                            pattern = paste0('^', pref, '.*_flank_', n.flank, '.*_aligned\\.fasta$'))
  files.mafft2 = list.files(path = path.mafft.out, 
                            pattern = paste0('^', pref, '.*_flank_', n.flank, '.*_aligned2\\.fasta$'))
  
  files.mafft1.pref = sapply(files.mafft1, function(s) strsplit(s, '_aligned')[[1]][1])
  files.mafft2.pref = sapply(files.mafft2, function(s) strsplit(s, '_aligned')[[1]][1])
  
  idx1.dup = which(files.mafft1.pref %in% files.mafft2.pref)
  if(length(idx1.dup) > 0){
    pokazAttention('duplicated alignments were found')
    files.mafft1 = files.mafft1[-idx1.dup]
  }
  
  files.mafft = c(files.mafft1,
                  files.mafft2)
  
  
  mafft.res = data.frame(file = files.mafft)
  
  if(nrow(mafft.res) > 0){
    y = lapply(mafft.res$file, function(s) strsplit(s, '_')[[1]][1:6])
    z = t(matrix(unlist(y), ncol = length(y)))
    mafft.res$comb = paste(z[,2], z[,3], sep = '_')
    mafft.res$beg = as.numeric(z[,5])
    mafft.res$end = as.numeric(z[,6])
    mafft.res$id = as.numeric(z[,3])  
  }
  
  # ---- Get long alignment positions ----
  pokaz('Read Long alignments..', file=file.log.main, echo=echo.main)
  idx.skip = c()
  mafft.aln.pos = list()
  
  pokaz('Number of mafft files', nrow(mafft.res), file=file.log.main, echo=echo.main)
  if(nrow(mafft.res) > 0){

    for(i in 1:nrow(mafft.res)){
      
      # pokaz('Aln', i, mafft.res$file[i], file=file.log.main, echo=echo.main)
      
      file.aln = paste0(path.mafft.out, mafft.res$file[i])
      
      if(!file.exists(file.aln)) {
        idx.skip = c(idx.skip, i)
        next
      }
      
      if (grepl("_aligned2", mafft.res$file[i])) {
        remove.flank = F
      } else {
        remove.flank = T
      }
      
      aln.seq = readFasta(file.aln)
      
      
      # ---
      # WITHOUT REFINEMENT
      n.aln.seq = length(aln.seq)
      name.aln.seq = names(aln.seq)
      name.acc = sapply(name.aln.seq, function(s) strsplit(s, '\\|')[[1]][1])
      pos.aln = sapply(name.aln.seq, function(s) strsplit(s, '\\|')[[1]][3:4])
      
      
      len.aln.seq <- nchar(aln.seq[1])
      # aln.mx <- matrix(, nrow = n.aln.seq, ncol = len.aln.seq)
      pos.mx <- matrix(0, nrow = n.aln.seq, ncol = len.aln.seq)
      for (i.seq in 1:length(aln.seq)) {
        tmp = strsplit(aln.seq[i.seq], "")[[1]]
        tmp.nongap = which(tmp != '-')
        
        if(remove.flank){
          tmp.nongap = tmp.nongap[-(1:(n.flank))]
          tmp.nongap <- tmp.nongap[1:(length(tmp.nongap) - n.flank)]  
        }
        
        p1 = pos.aln[1, i.seq]
        p2 = pos.aln[2, i.seq]
        pos.tmp = p1:p2
        pos.mx[i.seq, tmp.nongap] = pos.tmp
      }
      # aln.mx = aln.mx[,colSums(aln.mx != '-') != 0]
      pos.mx = pos.mx[,colSums(pos.mx != 0) != 0, drop=F]
      row.names(pos.mx) = name.acc
      
      if(min(colSums(pos.mx != 0)) == 0){
        pokaz(mafft.res$file[i])
        stop('Zeros in the alignment')
      }
      
      
      mafft.aln.pos[[mafft.res$file[i]]] = pos.mx
      # ---
      
    }
    
    if(length(idx.skip) > 0){
      mafft.aln.pos = mafft.aln.pos[-idx.skip]
      mafft.res = mafft.res[-idx.skip,,drop=F]
    }
    
    # warnings()
    mafft.res$len = unlist(lapply(mafft.aln.pos, ncol))
    mafft.res$extra = mafft.res$len - (mafft.res$end - mafft.res$beg - 1)
    
    # Skip if some are shorter than the initial aligned block
    idx.confusing = which(mafft.res$extra < 0)
    if(length(idx.confusing) > 0){
      pokazAttention('Long: Wrong lengths of alignment and gaps')
      pokazAttention("N confusing:", length(idx.confusing))
      
      mafft.res = mafft.res[-idx.confusing,,drop=F]
      mafft.aln.pos = mafft.aln.pos[-idx.confusing]
    } 
    # if(min(mafft.res$extra) < 0) stop()
    mafft.res$extra[mafft.res$extra < 0] = 0
  } else {
    mafft.res = data.frame()
  }
  
  # ---- Short alignments ----
  pokaz('Read Short alignments..', file=file.log.main, echo=echo.main)
  file.msa.res = paste0(path.inter.msa, 'aln_short_', s.comb, '.rds')
  if(file.exists(file.msa.res)){
    msa.res = readRDS(file.msa.res)
    msa.res$len = unlist(lapply(msa.res$aln, nrow))
    msa.res$extra = msa.res$len - (msa.res$ref.pos$end - msa.res$ref.pos$beg - 1)
    
    idx.confusing = which(msa.res$extra < 0)
    if(length(idx.confusing) > 0){
      pokazAttention('Short: Wrong lengths of alignment and gaps')
      pokazAttention("Confusing:", idx.confusing)
      
      msa.res$ref.pos = msa.res$ref.pos[-idx.confusing,,drop=F]
      msa.res$aln = msa.res$aln[-idx.confusing]
      msa.res$len = msa.res$len[-idx.confusing]
      msa.res$extra = msa.res$extra[-idx.confusing]
    } 
  
  } else {
    msa.res = data.frame()
  }

  # ---- Singletons alignments ----
  pokaz('Read Singletons..', file=file.log.main, echo=echo.main)
  file.single.res = paste0(path.inter.msa, 'singletons_', s.comb, '.rds')
  if(file.exists(file.single.res)){
    single.res = readRDS(file.single.res)
    
    # save(list = ls(), file = "tmp_workspace_x.RData")
    
    if(is.null(dim(single.res$pos.end))){
      single.res$pos.end = t(as.matrix(single.res$pos.end))
    }
    
    if(is.null(dim(single.res$pos.beg))){
      single.res$pos.beg = t(as.matrix(single.res$pos.beg))
    }
    
    single.res$len = rowSums(single.res$pos.end) - rowSums(single.res$pos.beg)  + 1  
    
    single.res$extra = single.res$len - (single.res$ref.pos$end - single.res$ref.pos$beg - 1) - 2
    
    idx.confusing = which((single.res$ref.pos$end - single.res$ref.pos$beg - 1) != 0)
    if(length(idx.confusing) > 0){
      pokazAttention('Singletons: Wrong lengths of alignment and gaps')
      pokazAttention("Number of confusing singletons:", length(idx.confusing))
      
      single.res$len = single.res$len[-idx.confusing]
      single.res$extra = single.res$extra[-idx.confusing]
      
      single.res$pos.beg = single.res$pos.beg[-idx.confusing, ,drop = F]
      single.res$pos.end = single.res$pos.end[-idx.confusing, ,drop = F]
      single.res$ref.pos = single.res$ref.pos[-idx.confusing, ,drop = F]
    }
      
  } else {
    single.res = data.frame()
  }

  pokaz(1)
  
  # ---- Analysis of positions ----
  # Here I wouls like fo find function of positions corresponcences between 4 things: 
  # old coordinates, long, short and singleton coordinates
  
  n.shift = rep(0, base.len)
  
  n.shift[mafft.res$end] = mafft.res$extra  # Long extra
  n.shift[msa.res$ref.pos$end] = msa.res$extra  # Short extra
  n.shift[single.res$ref.pos$end] = single.res$extra # Singletons extra
  n.shift = cumsum(n.shift)
  
  fp.main = (1:base.len) + n.shift
  
  pokaz(2)
  
  # -- 
  # Singletons
  fp.single = list()
  if(length(single.res) != 0){
    if(length(single.res$len) != 0){
      for(i in 1:length(single.res$len)){
        n.pos = single.res$len[i] - 2
        fp.single[[i]] = fp.main[single.res$ref.pos$beg[i]] + (1:n.pos)
      } 
    }
  }

  pokaz(3)
  # save(list = ls(), file = "tmp_workspace_3_oo.RData")
  
  # Short
  fp.short = list()
  if(length(msa.res) != 0){
    for(i in 1:length(msa.res$len)){
      n.pos = msa.res$len[i]
      fp.short[[i]] = fp.main[msa.res$ref.pos$beg[i]] + (1:n.pos)
    }    
    pokaz(4)
    
    # Check short
    for(i in 1:length(msa.res$len)){
      if(is.null(msa.res$aln[[i]])) next
      if(length(fp.short[[i]]) != nrow(msa.res$aln[[i]])){
        save(list = ls(), file = "tmp_workspace.RData")
        stop(paste0('Short', i)) 
      }
    }
  }

  pokaz(5)
  
  # Long
  fp.long = list()
  if(length(mafft.aln.pos) != 0){
    for(i in 1:length(mafft.aln.pos)){
      n.pos = ncol(mafft.aln.pos[[i]])
      fp.long[[i]] = fp.main[mafft.res$beg[i]] + (1:n.pos)
    }  
  }
  
  pos.beg.all = list(single.res$ref.pos$beg, msa.res$ref.pos$beg, mafft.res$beg)
  pos.end.all = list(single.res$ref.pos$end, msa.res$ref.pos$end, mafft.res$end)
  
  pos.delete.all = 0
  for(i.pos in 1:3){
    pos.beg = pos.beg.all[[i.pos]]
    pos.end = pos.end.all[[i.pos]]
    pos.delete = rep(0, base.len)
    pos.delete[pos.beg] = 1
    pos.delete[pos.end] = pos.delete[pos.end] - 1
    pos.delete = cumsum(pos.delete)
    pos.delete[pos.beg] = 0
    pos.delete[pos.end] = 0
    
    pos.delete.all = pos.delete.all + pos.delete
  }
  pos.delete = pos.delete.all

  if(sum(unlist(pos.end.all) - unlist(pos.beg.all) - 1) != sum(pos.delete)){
    save(list = ls(), file = 'tmx_workspace_step18.RData')
    stop('Wrong identification of positions to delete')
  } 
  
  pos.remain = pos.delete == 0
  
  fp.main[pos.delete != 0] = 0
  
  base.len.aln = max(fp.main)
  
  # Check-points
  fp.add = c(unlist(fp.single), unlist(fp.short), unlist(fp.long))
  if(sum(duplicated(c(fp.main[fp.main != 0], fp.add))) != 0) {
    save(list = ls(), file = paste0("tmp_workspace1_",s.comb,"_pointa.RData"))
    stop('Something is wrong with positions; Point A')
  } 
  # if(length(unique(c(fp.main, fp.add))) != (max(fp.main) + 1)) stop('Something if wrotng with positions; Point B')  # it's not trow anymore
  
  
  # ---- Define blocks before the big alignments ----
  # pos.beg = fp.main[mafft.res$beg] + 1
  # pos.beg.bins <- cut(pos.beg, breaks = c(seq.blocks, Inf), labels = FALSE)
  # pos.block.end = tapply(pos.beg, pos.beg.bins, max)
  # pos.block.end[length(pos.block.end)] = base.len
  
  file.res = paste0(path.features.msa, aln.type.out, s.comb,'.h5')
  if (file.exists(file.res)) file.remove(file.res)
  h5createFile(file.res)
  
  suppressMessages({
    h5createGroup(file.res, gr.accs.e)
  })
  
  pos.nonzero = fp.main != 0
  for(acc in accessions){
    pokaz('Accession', acc, file=file.log.main, echo=echo.main)
    v = h5read(file.comb, paste0(gr.accs.e, acc))
    v.aln = rep(0, base.len.aln)
    v.aln[fp.main[pos.nonzero]] = v[pos.nonzero]
    # v.aln[fp.main[pos.remain]] = v[pos.remain]
    
    # save(list = ls(), file = 'tmx_workspace_step18.RData')
    
    # Add singletons
    if(length(single.res$len) != 0){
      for(i in 1:length(single.res$len)){
        if(single.res$pos.beg[i, acc] != 0){
          pos = single.res$pos.beg[i, acc]:single.res$pos.end[i, acc]
          pos = pos[-c(1, length(pos))]
          
          v.aln[fp.single[[i]]] = pos
          # if(length(unique(v.aln)) != (sum(v.aln != 0) + 1)) stop('1')
        } 
      }   
      v.aln.nozero = v.aln[v.aln != 0]
      if(length(unique(v.aln.nozero)) != (sum(v.aln.nozero != 0))){
        save(list = ls(), file = "tmp_workspace.RData")
        stop('1: Duplicated positions in Singletons')
      } 
    }

    # Add short
    if(length(msa.res$len) != 0){
      for(i in 1:length(msa.res$len)){
        if(acc %in% colnames(msa.res$aln[[i]])){
          v.aln[fp.short[[i]]] = msa.res$aln[[i]][,acc]
        } 
      }
      v.aln.nozero = v.aln[v.aln != 0]
      if(length(unique(v.aln.nozero)) != (sum(v.aln.nozero != 0))){
        save(list = ls(), file = "tmp_workspace.RData")
        stop('2: Duplicated positions in short alignments')
      }       
    }

    
    # add long
    if(length(mafft.aln.pos) != 0){
      for(i in 1:length(mafft.aln.pos)){
        if(acc %in% rownames(mafft.aln.pos[[i]])){
          v.aln[fp.long[[i]]] = mafft.aln.pos[[i]][acc,]
        } 
      }      
      v.aln.nozero = v.aln[v.aln != 0]
      if(length(unique(v.aln.nozero)) != (sum(v.aln.nozero != 0))){
        save(list = ls(), file = "tmp_workspace.RData")
        stop('3: Duplicated positions in long alignments')
      } 
    }
    
    # Maybe something was overlapped by accident
    
    suppressMessages({
      h5write(v.aln, file.res, paste0(gr.accs.e, acc))
      
      stat.comb <- rbind(stat.comb, 
                         data.frame(acc = acc, 
                                    coverage = sum(v.aln != 0),
                                    comb = s.comb))
      
    })
    
  }
  H5close()
  gc()
  
}

warnings()

saveRDS(stat.comb, paste0(path.inter.msa, 'stat_coverage.rds'))

pokaz('Done.', file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- Manual testing ----



