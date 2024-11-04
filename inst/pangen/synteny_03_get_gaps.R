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
  files.maj <- files.blast[sapply(files.maj, function(x) any(sapply(accessions, function(a) startsWith(x, a))))]
  
  pokaz('All blast files', length(files.blast), file=file.log.main, echo=echo.main)
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

pokaz(files.maj)

if(length(files.maj) == 0) stop('No alignments provided')
pokaz('Number of alignments:', length(files.maj), file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(f.maj, 
                          echo.loop=T){
  
  # Remove extensions
  pref.comb <- sub("\\_maj.rds$", "", f.maj)
  
  # Log files
  file.log.loop = paste0(path.log, 'loop_file_', 
                         pref.comb,
                         '.log')
  if(!file.exists(file.log.loop)){
    invisible(file.create(file.log.loop))
  }
  
  # ---- Check log Done ----
  if(checkDone(file.log.loop)){
    return()
  }

  # --- --- --- --- --- --- --- --- --- --- ---
  
  # Parser prefix of the result file
  info <- pref2info(pref.comb)
  query.chr <- info$query.chr
  base.chr <- info$base.chr
  acc <- info$acc
  
  pokaz(acc, query.chr, base.chr, file=file.log.loop, echo=echo.loop)
  
  # Read reference sequences
  base.file = paste0(base.acc, '_chr', base.chr , '.', 'fasta', collapse = '')
  pokaz('Base:', base.file, file=file.log.loop, echo=echo.loop)
  base.fas.fw = readFastaMy(paste0(path.chr, base.file))
  base.fas.fw = seq2nt(base.fas.fw)
  base.fas.bw = revCompl(base.fas.fw)
  base.len = length(base.fas.bw)
  pokaz('Length of base:', base.len, file=file.log.loop, echo=echo.loop)
  
  # Read query sequences
  query.file = paste0(acc, '_chr',query.chr, '.fasta')
  pokaz('Query:', query.file, file=file.log.loop, echo=echo.loop)
  
  query.fas.chr = readFastaMy(paste0(path.chr, query.file))
  query.fas.chr = seq2nt(query.fas.chr)
  query.len = length(query.fas.chr)
  pokaz('Length of query:', query.len, file=file.log.loop, echo=echo.loop)
  
  x = readRDS(paste0(path.aln, f.maj))
  
  if((nrow(x) <= 1) || (is.null(x))) {
    pokaz('No gaps', file=file.log.loop, echo=echo.loop)
    
    rmSafe(x)
    rmSafe(x.major)
    rmSafe(base.fas.bw)
    rmSafe(base.fas.fw)
    rmSafe(query.fas.chr)
    gc()
    
    pokaz('Done.', file=file.log.loop, echo=echo.loop)
    return(NULL)
  }
  
  # ---- Get gaps ----
  pokaz('Get gaps', file=file.log.loop, echo=echo.loop)
  
  checkCorrespToGenome(x=setDir(x, base.len = base.len),
                       query.fas = query.fas.chr, 
                       base.fas.fw = base.fas.fw, 
                       base.fas.bw = base.fas.bw)
  
  # Find occupied positions
  pos.q.free = rep(0, query.len)  # free positions in query
  pos.b.free = rep(0, base.len)  # free positions in base
  for(irow in 1:nrow(x)){
    pos.q.free[x$V2[irow]:x$V3[irow]] <- pos.q.free[x$V2[irow]:x$V3[irow]] + 1
    pos.b.free[x$V4[irow]:x$V5[irow]] <- pos.b.free[x$V4[irow]:x$V5[irow]] + 1
  }
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
  # Query-file
  file.gap.query = paste0(path.gaps, pref.comparisson, 'query.fasta', collapse = '')
  # Base file
  file.gap.base = paste0(path.gaps, pref.comparisson, 'base.fasta', collapse = '')
  pokaz('Create gaps for', file.gap.query, file=file.log.loop, echo=echo.loop)
  pokaz('Create gaps for', file.gap.base, file=file.log.loop, echo=echo.loop)
  
  for(irow in 1:(nrow(x)-1)){
    
    # Common file name
    pref.gap = paste0('gap_', irow, '_', irow + 1, '_')
    
    # If from another block - don't consider the gap
    if(x$bl[irow] != x$bl[irow+1]) next
    
    pos = sort(c(x$V2[irow], x$V3[irow], x$V2[irow+1], x$V3[irow+1]))
    d1 = pos[3] - pos[2]
    
    pos = sort(c(x$V4[irow], x$V5[irow], x$V4[irow+1], x$V5[irow+1]))
    d2 = pos[3] - pos[2]
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
    
    # Create files for BLAST
    pos = sort(c(x$V2[irow], x$V3[irow], x$V2[irow+1], x$V3[irow+1]))
    pos.gap.q = pos[2]:pos[3]
    pos = sort(c(x$V4[irow], x$V5[irow], x$V4[irow+1], x$V5[irow+1]))
    pos.gap.b = pos[2]:pos[3]
    
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
    
    # ---- Write query ----
    # Define Chunks
    s.q = query.fas.chr[pos.gap.q]
    s.q = nt2seq(s.q)
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
    
    # Standsrd naming (as before)
    pref.q = paste(pref.comparisson, pref.gap,
                   'query', '|', pos.gap.q[p.beg], '|', pos.gap.q[p.end], sep = '')
    
    if(sum(pos.gap.q[p.beg] > pos.gap.q[p.end]) > 0) stop('Wrong boundaries of gap blocks - 2')
    
    if(length(s.q) != length(pref.q)) stop('Chunk lengths do not much')
    names(s.q) = pref.q
    writeFastaMy(s.q, file.gap.query, append = T)
    
    # ---- Write base ----
    s.b = base.fas.fw[pos.gap.b]
    s.b = nt2seq(s.b)
    
    s.base.names = paste(pref.comparisson, pref.gap,
                         'base', '|', pos.gap.b[1], '|', pos.gap.b[length(pos.gap.b)], sep = '')
    
    names(s.b) = s.base.names
    
    writeFastaMy(s.b, file.gap.base, append = T)
    
  }  # irow search for gaps
  
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
  
  if(nrow(diffs) > 0){
    for(irow in 1:nrow(diffs)){
      
      if((diffs$end[irow] - diffs$beg[irow]) < len.blast) next
      if((diffs$end[irow] - diffs$beg[irow]) > max.len) next
      pos.gap.q = diffs$beg[irow]:diffs$end[irow]
      
      irow.prev = which(x$V2 <  diffs$beg[irow])  # x is sorted by V2
      irow.next = which(x$V3 >  diffs$end[irow])
      
      if(length(irow.prev) == 0){
        irow.prev = 0
      } else {
        irow.prev = x$id[max(irow.prev)]
      }
      
      if(length(irow.next) == 0){
        irow.next = nrow(x) + 1
      } else {
        irow.next = x$id[min(irow.next)]
      }
      
      pref.gap = paste0('connect_', irow.prev, '_', irow.next, '_')
      
      # Define Chunks
      s.q = query.fas.chr[pos.gap.q]
      s.q = nt2seq(s.q)
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
      writeFastaMy(s.q, file.gap.query, append = T)
    }
  }
  
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
  
  
  if(nrow(diffs) > 0){
    for(irow in 1:nrow(diffs)){
      
      if((diffs$end[irow] - diffs$beg[irow]) < len.blast) next
      if((diffs$end[irow] - diffs$beg[irow]) > max.len) next
      pos.gap.b = diffs$beg[irow]:diffs$end[irow]
      
      irow.prev = which(x$V4 <  diffs$beg[irow])  # x is sorted by V2
      irow.next = which(x$V5 >  diffs$end[irow])
      
      if(length(irow.prev) == 0){
        irow.prev = 0
      } else {
        irow.prev = x$id[max(irow.prev)]
      }
      
      if(length(irow.next) == 0){
        irow.next = nrow(x) + 1
      } else {
        irow.next = x$id[min(irow.next)]
      }
      
      pref.gap = paste0('connect_', irow.prev, '_', irow.next, '_')
      
      
      s.b = base.fas.fw[pos.gap.b]
      
      if((sum(s.b == 'N') + sum(s.b == 'n')) > length(s.b) / 2) next
      
      s.b = nt2seq(s.b)
      
      s.base.names = paste(pref.comparisson, pref.gap,
                           'resid_base', '|', pos.gap.b[1], '|', pos.gap.b[length(pos.gap.b)], sep = '')
      
      names(s.b) = s.base.names
      
      writeFastaMy(s.b, file.gap.base, append = T)
    }
  }
  
  rmSafe(x)
  rmSafe(x.major)
  rmSafe(base.fas.bw)
  rmSafe(base.fas.fw)
  rmSafe(query.fas.chr)
  gc()
  
  pokaz('Done.', file=file.log.loop, echo=echo.loop)
  return(NULL)
}

# ***********************************************************************
# ---- Loop  ----


if(num.cores == 1){
  for(f.maj in files.maj){
    loop.function(f.maj, echo.loop=echo.loop)
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
    
    # Run parallel loop for files in the current batch
    batch.results <- foreach(f.maj = batch.files, .packages = c('crayon'), .verbose = FALSE) %dopar% {
      loop.function(f.maj, echo.loop = echo.loop)
    }
    
    tmp <- c(tmp, batch.results)  # Store the batch results in the main list
    
    stopCluster(myCluster)  # Stop the cluster after completing the batch
  }
}

warnings()

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


