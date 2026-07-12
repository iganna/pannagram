# Alignment-1. Remaining syntenic (major) matches

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/synteny_func.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("--path.chr"),    type = "character", default = NULL, help = "Path to query chromosome fasta files"),
  make_option(c("--path.aln"),    type = "character", default = NULL, help = "Path to the output directory with alignments"),
  make_option(c("--path.gaps"),   type = "character", default = NULL, help = "Path to the directory with gaps"),
  
  make_option(c("--ref"),         type = "character", default = NULL, help = "Name of the reference genome"),
  make_option(c("--accessions"),    type = "character", default = NULL, help = "File containing accessions to analyze"),
  make_option(c("--combinations"),  type = "character", default = NULL, help = "File containing combinations to analyze"),
  
  make_option(c("--cores"),       type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option(c("--path.log"),    type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"),   type = "character", default = NULL, help = "Level of log to be shown on the screen")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# TODO:
max.len = 10^6
len.blast = 50

# Validating that the alignment coordinates correspond to the genome sequence is
# a debug-only check (consistent with synteny_05_merge_gaps.R, where it is off by
# default). When FALSE, the whole chromosome is never exploded into a per-nucleotide
# vector and its reverse complement is never computed -- gaps are sliced with substr().
check.genomes = F


# print(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- Values of parameters ----

# Number of cores
num.cores <- opt$cores

path.chr      <- ifelse(!is.null(opt$path.chr), opt$path.chr, stop('Folder with chromosomes is not specified'))
path.aln      <- ifelse(!is.null(opt$path.aln), opt$path.aln, stop('Folder with Alignments is not specified'))
base.acc      <- ifelse(!is.null(opt$ref), opt$ref, stop('Reference genome is not specified'))
path.gaps     <- ifelse(!is.null(opt$path.gaps), opt$path.gaps, stop('Folder with Gaps is not specified'))

# Create folders for the alignment results
if(!dir.exists(path.aln)) dir.create(path.aln)
if(!dir.exists(path.gaps)) dir.create(path.gaps)

# Accessions
file.acc <- ifelse(!is.null(opt$accessions), opt$accessions, stop("File with accessions are not specified"))
tmp <- read.table(file.acc, stringsAsFactors = F)
accessions <- as.character(tmp[,1])
pokaz('Names of genomes for the analysis:', accessions, 
      file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- Combinations ----

file.combinations <- ifelse(!is.null(opt$combinations), opt$combinations, stop("File with combinations are not specified"))
if (length(readLines(file.combinations)) == 0) {
  files.maj <- list.files(path.aln, pattern = "\\maj.rds$")
  
  # Filter files that start with one of the values in accessions
  files.maj <- files.maj[sapply(files.maj, function(x) any(sapply(accessions, function(a) startsWith(x, a))))]
  
  pokaz('All blast files', length(files.maj), file=file.log.main, echo=echo.main)
} else {
  combinations = read.table(file.combinations)
  files.maj = c()
  for(acc in accessions){
    for(i.comb in 1:nrow(combinations)){
      file.tmp = paste0(acc, '_', 
                        combinations[i.comb,1], '_',
                        combinations[i.comb,2], '_maj.rds')
      if(file.exists(paste0(path.aln, file.tmp))){
        files.maj = c(files.maj, file.tmp)  
      }
    }
  }
}

pokaz(files.maj, file=file.log.main, echo=echo.main)

if(length(files.maj) == 0) stop('No alignments provided')
pokaz('Number of alignments:', length(files.maj), file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(f.maj,
                          done.set = character(0),
                          echo.loop=T){
  initial.vars <- ls()

  pref.comb <- sub("\\_maj.rds$", "", f.maj)

  # ---- Checkpoint: skip already completed items ----
  item.id <- pref.comb
  if(item.id %in% done.set){
    return(NULL)
  }

  # One log file per worker (bounded number of files); also the checkpoint ledger
  file.log.loop = initLoopLog(path.log)

  pokaz("File", f.maj, file=file.log.loop, echo=echo.loop)
  # Remove extensions
  

  # --- --- --- --- --- --- --- --- --- --- ---
  
  # Parser prefix of the result file
  info <- pref2info(pref.comb)
  query.chr <- info$query.chr
  base.chr <- info$base.chr
  acc <- info$acc
  
  pokaz(acc, query.chr, base.chr, file=file.log.loop, echo=echo.loop)
  
  # Read reference sequence (kept as a single string; gap sub-sequences are
  # sliced with substr() below, so the whole chromosome is never turned into a
  # per-nucleotide vector unless the genome-correspondence check is enabled).
  base.file = paste0(base.acc, '_chr', base.chr , '.', 'fasta', collapse = '')
  pokaz('Base:', base.file, file=file.log.loop, echo=echo.loop)
  base.str = unname(readFastaMy(paste0(path.chr, base.file))[1])
  base.len = nchar(base.str)
  pokaz('Length of base:', base.len, file=file.log.loop, echo=echo.loop)

  # Read query sequence (also kept as a string)
  query.file = paste0(acc, '_chr',query.chr, '.fasta')
  pokaz('Query:', query.file, file=file.log.loop, echo=echo.loop)

  query.str = unname(readFastaMy(paste0(path.chr, query.file))[1])
  query.len = nchar(query.str)
  pokaz('Length of query:', query.len, file=file.log.loop, echo=echo.loop)

  # Per-nucleotide vectors + reverse complement are only needed for the optional
  # genome-correspondence validation below.
  if(check.genomes){
    base.fas.fw = seq2nt(base.str)
    base.fas.bw = revCompl(base.fas.fw)
    query.fas.chr = seq2nt(query.str)
  }
  
  x = readRDS(paste0(path.aln, f.maj))
  x = cleanOverlaps(x)
  
  # save(list = ls(), file = "tmp_workspace.RData")
  
  if((nrow(x) <= 1) || (is.null(x))) {
    pokaz('No gaps', file=file.log.loop, echo=echo.loop)
    
    # ---- Checkpoint marker: item fully processed ----
    markDone(item.id, file=file.log.loop, echo=echo.loop)
    
    # Cleanup variables
    final.vars <- ls()
    new.vars <- setdiff(final.vars, initial.vars)
    rm(list = new.vars)
    gc()
    return(NULL)
  }
  
  # ---- Get gaps ----
  pokaz('Get gaps', file=file.log.loop, echo=echo.loop)
  
  # save(list = ls(), file = "tmp_workspace_get_gap.RData")
  
  if(check.genomes){
    checkCorrespToGenome(x=setDir(x, base.len = base.len),
                         query.fas = query.fas.chr,
                         base.fas.fw = base.fas.fw,
                         base.fas.bw = base.fas.bw)
  }
  
  # Find occupied positions
  pos.q.free = rep(0, query.len)  # free positions in query
  pos.b.free = rep(0, base.len)  # free positions in base

  for(irow in 1:nrow(x)){
    pos.q.free[x$V2[irow]:x$V3[irow]] <- pos.q.free[x$V2[irow]:x$V3[irow]] + 1
    pos.b.free[x$V4[irow]:x$V5[irow]] <- pos.b.free[x$V4[irow]:x$V5[irow]] + 1
  }
  
  # save(list = ls(), file = "tmp_workspace.RData")
  
  if((sum(pos.q.free > 1) != 0)) stop('Coverage query is wrong')
  if((sum(pos.b.free > 1) != 0)) stop('Coverage base is wrong')
  
  pos.q.free = -pos.q.free
  pos.b.free = -pos.b.free
  
  # ---- Masking ----
  # Accession
  file.out.masking = paste0(path.chr, 'mask_', acc, '_chr', query.chr, '.rds', collapse = '')
  if(file.exists(file.out.masking)){
    pokaz('Read masking of the accession', file.out.masking, file=file.log.loop, echo=echo.loop)
    pos.masking = readRDS(file.out.masking)
    for(irow in 1:nrow(pos.masking)){
      pos.q.free[pos.masking$beg[irow]:pos.masking$end[irow]] = -1
    }
  }
  
  file.out.masking = paste0(path.chr, 'mask_', base.acc, '_chr', base.chr, '.rds', collapse = '')
  if(file.exists(file.out.masking)){
    pokaz('Read masking of the reference', file.out.masking, file=file.log.loop, echo=echo.loop)
    pos.masking = readRDS(file.out.masking)
    for(irow in 1:nrow(pos.masking)){
      pos.b.free[pos.masking$beg[irow]:pos.masking$end[irow]] = -1
    }
  }
  
  # Reference
  
  # ---- Write gaps ----
  
  # Within non-occupied positions find those, which can be
  
  pref.comparisson = paste0('acc_', acc, '_qchr_', query.chr, '_bchr_', base.chr, '_')

  # Idempotent: gap fasta files below are written with append=T, so drop any partial
  # output left by a previously interrupted run of this item before re-creating it.
  # Covers the batched normal files (<pref>b<N>_query.fasta / _base.fasta) and residuals.
  all.gap.files = list.files(path.gaps, full.names = TRUE)
  rm.gap.files = all.gap.files[startsWith(basename(all.gap.files), pref.comparisson) &
                               (endsWith(all.gap.files, 'query.fasta') | endsWith(all.gap.files, 'base.fasta'))]
  if(length(rm.gap.files) > 0) invisible(file.remove(rm.gap.files))

  # Normal gaps are written in BATCHES of `batch.size.gaps` gaps: each batch becomes its
  # own <pref>b<N>_query.fasta / _base.fasta so step-7 blasts a query gap only against its
  # ~batch.size.gaps neighbours instead of the whole all-vs-all base. The batch target file
  # names are set per gap in the loop below.
  #
  # Larger batches => FEWER step-7 blastn invocations => less per-process startup overhead
  # (the dominant single-core cost of step 7: ~0.2s x n_batches). This is CORRECTNESS-NEUTRAL:
  # a query gap and its paired base gap are always written to the SAME batch, and synteny_05
  # keeps a hit only when pref1==pref2 (before any greedy selection), so cross-batch/cross-pref
  # hits are discarded and coverage is identical for any batch size. It is a pure timing knob:
  # bigger => step-7 faster, step-8 (merge) slightly heavier as it reads more discarded hits.
  # Kept at 100 (reverted a 500 trial: on real repeat-rich anopheles gaps the larger all-vs-all
  # per batch offset the fewer-startup savings and nudged this step slower). Sweep if you re-tune.
  batch.size.gaps = 100
  n.gap.written = 0

  # Accumulate gap sequences in memory and write each batch file ONCE after the
  # loop. Previously writeFastaMy(append=T) opened/closed the file on every gap
  # (~2 * #gaps file operations -> the dominant cost of this step).
  acc.q = list(); acc.b = list(); acc.batch = integer(0)

  # ---- Precompute gap geometry for all adjacent block pairs (vectorised) ----
  # Replaces ~4 per-gap sort(c(4 coords)) calls with element-wise 2nd/3rd order
  # statistics (pos[2]/pos[3] of 4 numbers) computed once over the whole vector.
  # For two within-block-sorted pairs, s2 = min(max(o1,o3), min(o2,o4)),
  # s3 = max(max(o1,o3), min(o2,o4)); identical to sort(c(a,b,c,d))[2:3].
  nx = nrow(x); ip = 1:(nx - 1)
  pref.gaps.v = paste0('gap_', ip, '_', ip + 1, '_')
  qo1 = pmin(x$V2[ip], x$V3[ip]);     qo2 = pmax(x$V2[ip], x$V3[ip])
  qo3 = pmin(x$V2[ip+1], x$V3[ip+1]); qo4 = pmax(x$V2[ip+1], x$V3[ip+1])
  q.s2.v = pmin(pmax(qo1, qo3), pmin(qo2, qo4))
  q.s3.v = pmax(pmax(qo1, qo3), pmin(qo2, qo4))
  d1.v = q.s3.v - q.s2.v
  bo1 = pmin(x$V4[ip], x$V5[ip]);     bo2 = pmax(x$V4[ip], x$V5[ip])
  bo3 = pmin(x$V4[ip+1], x$V5[ip+1]); bo4 = pmax(x$V4[ip+1], x$V5[ip+1])
  b.s2.v = pmin(pmax(bo1, bo3), pmin(bo2, bo4))
  b.s3.v = pmax(pmax(bo1, bo3), pmin(bo2, bo4))
  d2.v = b.s3.v - b.s2.v

  for(irow in 1:(nrow(x)-1)){

    # Common file name (precomputed)
    pref.gap = pref.gaps.v[irow]

    # If from another block - don't consider the gap
    if(x$bl[irow] != x$bl[irow+1]) next

    d1 = d1.v[irow]
    d2 = d2.v[irow]
    if(d1 * d2 == 0) next  # Glue Zero, but not necessary
    
    # Short gaps will be covered later!!!
    # if((d1 <= len.blast) & (d2 <= len.blast)){
    #   x.tmp = x[c(irow,irow+1),]
    #   x.tmp = setDir(x.tmp, base.len = base.len)
    #   x.tmp = glueByThreshold(x.tmp, 100, query.fas = query.fas.chr,
    #                       base.fas.fw = base.fas.fw,
    #                       base.fas.bw = base.fas.bw, file.log=NULL)
    #   x.tmp = getBase(x.tmp, base.len = base.len)
    #   x[irow+1,] = x.tmp
    #   idx.remove = c(idx.remove, irow)
    # }
    
    # Don't consider for blast short sequences
    if(d1 <= len.blast){
      pos.q.free[(x$V3[irow]+1):(x$V2[irow+1]-1)] = -1
    }
    if(d2 <= len.blast){
      pos.b.free[(x$V5[irow]+1):(x$V4[irow+1]-1)] = -1
    }
    
    if(!((d1 >= len.blast) & (d2 >= len.blast))) next
    
    # Create files for BLAST (gap boundaries precomputed above)
    pos.gap.q = q.s2.v[irow]:q.s3.v[irow]
    pos.gap.b = b.s2.v[irow]:b.s3.v[irow]
    
    # Remove flanking positions
    pos.gap.q = pos.gap.q[-c(1, length(pos.gap.q))]  
    pos.gap.b = pos.gap.b[-c(1, length(pos.gap.b))]
    
    # Save occupancy
    pos.q.free[pos.gap.q] = irow
    pos.b.free[pos.gap.b] = irow
    
    # if already occupied by masking - do not 
    
    # if(file.exists(file.gap.query)) next
    
    if(abs(pos.gap.q[1] - pos.gap.q[length(pos.gap.q)]) > max.len) next
    if(abs(pos.gap.b[1] - pos.gap.b[length(pos.gap.b)]) > max.len) next

    # ---- Batch index: rolls to a new file every batch.size.gaps written gaps ----
    i.batch = n.gap.written %/% batch.size.gaps

    # ---- Build query chunks (pos.gap.q is a contiguous ascending range -> substr slice) ----
    s.q = substr(query.str, pos.gap.q[1], pos.gap.q[length(pos.gap.q)])
    n.bl = 500
    len.s.q = nchar(s.q)
    if(len.s.q > n.bl){
      p.beg = seq(1, len.s.q, n.bl)
      p.end = seq(n.bl, len.s.q, n.bl)
      if(p.end[length(p.end)] != len.s.q) p.end = c(p.end, len.s.q)
      if(length(p.end) != length(p.beg)) stop('Wrong lengths of p.beg and p.end')
    } else {
      p.beg = 1
      p.end = len.s.q
    }
    s.q = splitSeq(s.q, n = n.bl)
    pref.q = paste(pref.comparisson, pref.gap,
                   'query', '|', pos.gap.q[p.beg], '|', pos.gap.q[p.end], sep = '')
    if(sum(pos.gap.q[p.beg] > pos.gap.q[p.end]) > 0) stop('Wrong boundaries of gap blocks - 2')
    if(length(s.q) != length(pref.q)) stop('Chunk lengths do not much')
    names(s.q) = pref.q

    # ---- Build base ----
    s.b = substr(base.str, pos.gap.b[1], pos.gap.b[length(pos.gap.b)])
    names(s.b) = paste(pref.comparisson, pref.gap,
                       'base', '|', pos.gap.b[1], '|', pos.gap.b[length(pos.gap.b)], sep = '')

    # ---- Accumulate; each batch file is written once, after the loop ----
    k = length(acc.q) + 1L
    acc.q[[k]] = s.q
    acc.b[[k]] = s.b
    acc.batch[k] = i.batch

    n.gap.written = n.gap.written + 1  # advances the batch index (every batch.size.gaps gaps)

  }  # irow search for gaps

  # ---- Flush accumulated gaps: ONE write per batch file (was append-per-gap) ----
  if(length(acc.q) > 0){
    for(b in sort(unique(acc.batch))){
      idx = which(acc.batch == b)
      writeFastaMy(do.call(c, acc.q[idx]),
                   paste0(path.gaps, pref.comparisson, 'b', b, '_query.fasta'))
      writeFastaMy(do.call(c, acc.b[idx]),
                   paste0(path.gaps, pref.comparisson, 'b', b, '_base.fasta'))
    }
  }

  # ---- Write remained blocks ----
  pokaz('Create fasta for the remained sequences', file=file.log.loop, echo=echo.loop)
  file.gap.query = paste0(path.gaps, pref.comparisson, 'residual_query.fasta', collapse = '')
  # Base file
  file.gap.base = paste0(path.gaps, pref.comparisson, 'residual_base.fasta', collapse = '')
  pokaz('Create gaps for', file.gap.query, file=file.log.loop, echo=echo.loop)
  pokaz('Create gaps for', file.gap.base, file=file.log.loop, echo=echo.loop)
  
  x$idx = 1:nrow(x)
  
  ## ---- Write query ----
  # Query: Zero-coverage blocks
  
  
  # file.ws = "tmp_workspace.RData"
  # all.local.objects <- ls()
  # save(list = all.local.objects, file = file.ws)
  # pokaz('Workspace is saved in', file.ws, file=file.log.loop, echo=echo.loop)
  # stop('Enough..')
  
  x = x[order(x$V2),]
  diffs = findOnes( (pos.q.free == 0) * 1)

  res.q = list()   # accumulate residual query seqs; one write after the loop
  if(nrow(diffs) > 0){
    for(irow in 1:nrow(diffs)){
      
      if((diffs$end[irow] - diffs$beg[irow]) < len.blast) next
      if((diffs$end[irow] - diffs$beg[irow]) > max.len) next
      pos.gap.q = diffs$beg[irow]:diffs$end[irow]
      
      n.prev = findInterval(diffs$beg[irow] - 1L, x$V2)  # x sorted by V2 -> count of V2 < beg == max(which(V2<beg))
      irow.next = which(x$V3 >  diffs$end[irow])

      if(n.prev == 0){
        irow.prev = 0
      } else {
        irow.prev = x$id[n.prev]
      }
      
      if(length(irow.next) == 0){
        irow.next = nrow(x) + 1
      } else {
        irow.next = x$id[min(irow.next)]
      }
      
      pref.gap = paste0('connect_', irow.prev, '_', irow.next, '_')

      # Define Chunks (pos.gap.q is a contiguous ascending range -> substr slice)
      s.q = substr(query.str, pos.gap.q[1], pos.gap.q[length(pos.gap.q)])
      n.bl = 500
      len.s.q = nchar(s.q)
      if(len.s.q > n.bl){
        p.beg = seq(1, len.s.q, n.bl)
        p.end = seq(n.bl, len.s.q, n.bl)
        if(p.end[length(p.end)] != len.s.q) p.end = c(p.end, len.s.q)
        if(length(p.end) != length(p.beg)) stop('Wrong lengths of p.beg and p.end')
      } else {
        p.beg = 1
        p.end = len.s.q
      }
      s.q = splitSeq(s.q, n = n.bl)
      
      # Standard naming (as before)
      pref.q = paste(pref.comparisson, pref.gap,
                     'resid_query', '|', pos.gap.q[p.beg], '|', pos.gap.q[p.end], sep = '')
      # pokaz('pos.gap.q', pos.gap.q, file=file.log.loop, echo=echo.loop)
      # pokaz('p.beg', p.beg, file=file.log.loop, echo=echo.loop)
      # pokaz('p.end', p.end, file=file.log.loop, echo=echo.loop)
      if(sum(pos.gap.q[p.beg] > pos.gap.q[p.end]) > 0) stop('Wrong boundaries of gap blocks - 2')
      
      if(length(s.q) != length(pref.q)) stop('Chunk lengths do not much')
      names(s.q) = pref.q
      res.q[[length(res.q) + 1L]] = s.q
    }
  }
  if(length(res.q) > 0) writeFastaMy(do.call(c, res.q), file.gap.query)

  ## ---- Write base ----
  # Query: Zero-coverage blocks
  
  # Switch V4 and V5
  idx = which(x$V4 > x$V5)
  tmp = x$V4[idx]
  x$V4[idx] = x$V5[idx]
  x$V5[idx] = tmp
  
  # Sorting
  x = x[order(x$V4),]
  diffs = findOnes( (pos.b.free == 0) * 1)

  res.b = list()   # accumulate residual base seqs; one write after the loop
  if(nrow(diffs) > 0){
    for(irow in 1:nrow(diffs)){
      
      if((diffs$end[irow] - diffs$beg[irow]) < len.blast) next
      if((diffs$end[irow] - diffs$beg[irow]) > max.len) next
      pos.gap.b = diffs$beg[irow]:diffs$end[irow]
      
      n.prev = findInterval(diffs$beg[irow] - 1L, x$V4)  # x sorted by V4 -> count of V4 < beg == max(which(V4<beg))
      irow.next = which(x$V5 >  diffs$end[irow])

      if(n.prev == 0){
        irow.prev = 0
      } else {
        irow.prev = x$id[n.prev]
      }
      
      if(length(irow.next) == 0){
        irow.next = nrow(x) + 1
      } else {
        irow.next = x$id[min(irow.next)]
      }
      
      pref.gap = paste0('connect_', irow.prev, '_', irow.next, '_')
      
      
      s.b = substr(base.str, pos.gap.b[1], pos.gap.b[length(pos.gap.b)])

      # Skip blocks that are more than half N (count N/n directly on the string)
      if((nchar(s.b) - nchar(gsub('[Nn]', '', s.b))) > nchar(s.b) / 2) next

      s.base.names = paste(pref.comparisson, pref.gap,
                           'resid_base', '|', pos.gap.b[1], '|', pos.gap.b[length(pos.gap.b)], sep = '')
      
      names(s.b) = s.base.names
      
      res.b[[length(res.b) + 1L]] = s.b
    }
  }
  if(length(res.b) > 0) writeFastaMy(do.call(c, res.b), file.gap.base)

  # ---- Checkpoint marker: item fully processed ----
  markDone(item.id, file=file.log.loop, echo=echo.loop)

  # Cleanup variables
  final.vars <- ls()
  new.vars <- setdiff(final.vars, initial.vars)
  rm(list = new.vars)
  gc()
  return(NULL)
}

# ***********************************************************************

# ---- Build the set of completed items from all worker logs (any core count) ----
done.set <- getDoneSet(path.log)
files.maj = files.maj[!(sub("\\_maj.rds$", "", files.maj) %in% done.set)]

if(length(done.set) > 0){
  pokaz('Skip already done:', length(done.set), file=file.log.main, echo=echo.main)
}

if(length(files.maj) == 0){
  pokaz('All files are done', file=file.log.main, echo=echo.main)
  quit(save = "no")
}


# ---- Loop  ----

if(num.cores == 1){
  assign('.worker.id', 1, envir = .GlobalEnv)   # stable single file core_1.log
  for(f.maj in files.maj){
    loop.function(f.maj, done.set = done.set, echo.loop=echo.loop)
  }
} else {

  batch.size <- 2 * num.cores  # Define the batch size
  
  # Initialize a temporary list to store results
  tmp <- list()
  
  # Loop through files in batches
  for (start in seq(1, length(files.maj), by = batch.size)) {
    
    end <- min(start + batch.size - 1, length(files.maj))  # Define the end of the current batch
    batch.files <- files.maj[start:end]  # Subset files for the current batch
    
    # Create and register a new cluster for the current batch
    myCluster <- makeCluster(num.cores, type = "PSOCK")
    registerDoParallel(myCluster)

    # Assign a stable worker id (1..N) to each worker -> bounded, reused file names
    parallel::clusterApply(myCluster, seq_len(num.cores),
                           function(i) assign('.worker.id', i, envir = .GlobalEnv))

    pokaz(batch.files, file=file.log.main, echo=echo.main)

    # Run parallel loop for files in the current batch
    batch.results <- foreach(f.maj = batch.files, .packages = c('crayon', 'pannagram'), .verbose = FALSE) %dopar% {
      loop.function(f.maj, done.set = done.set, echo.loop = echo.loop)
    }
    
    tmp <- c(tmp, batch.results)  # Store the batch results in the main list
    
    stopCluster(myCluster)  # Stop the cluster after completing the batch
    gc()
  }
}


pokaz('Done.', file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- Manual testing ----

if(F){
source(system.file("pangen/synteny_func.R", package = "pannagram"))
source(system.file("utils/utils.R", package = "pannagram"))

  # file.ws = "tmp_workspace.RData"
  # all.local.objects <- ls()
  # save(list = all.local.objects, file = file.ws)
  # pokaz('Workspace is saved in', file.ws, file=file.log.loop, echo=echo.loop)
  # stop('Enough..')
}


