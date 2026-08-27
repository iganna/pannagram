suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
  library(pannagram)
})

source(system.file("pangen/comb_func.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.features.msa", type = "character", default = NULL, help = "Path to msa directory (features)"),
  make_option("--path.inter.msa",    type = "character", default = NULL, help = "Path to msa directory (internal)"),
  make_option("--path.chromosomes",  type = "character", default = NULL, help = "Path to directory with chromosomes"),

  make_option("--len.min.merge",     type = "integer",   default = NULL, help = "Min length of both breaks to be considered; defaults to len.short"),
  make_option("--gap.max.merge",     type = "integer",   default = 5000, help = "Max distance between two breaks"),
  make_option("--len.max.merge",     type = "integer",   default = 25000, help = "Max pangenome span of a merged break (same limit solveLong2 enforces in comb_04)"),
  make_option("--dist.tol",          type = "double",    default = 0.5,   help = "Max gap between two breaks as a share of their length; same meaning as dist.tol in mergeOverlapsTolerance()"),
  make_option("--len.max.sv",        type = "integer",   default = 15000, help = "Breaks longer than this are not merged; same meaning as len.max.sv in mergeOverlapsTolerance()"),
  make_option("--p.ident.merge",     type = "double",    default = 90,   help = "Min percent of identity between two breaks"),
  make_option("--cover.merge",       type = "double",    default = 0.5,  help = "Min coverage of the shorter break"),
  make_option("--member.tol",        type = "double",    default = 0.5,  help = "Min share of the break length for an accession to count"),

  make_option("--cores",             type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",          type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",         type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser, args = args)

#TODO: SHOULD BE PARAMATERS
len.short = 50

# ***********************************************************************
# ---- Logging ----
source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----
source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

aln.type.in = aln.type.clean
aln.type.in = paste0(aln.type.in, '_')

# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores = opt$cores
if(is.null(num.cores)) stop('Wrong number of cores: NULL')
pokaz('Number of cores', num.cores, file=file.log.main, echo=echo.main)
if(num.cores > 1){
  myCluster <- makeCluster(num.cores, type = "PSOCK")
  registerDoParallel(myCluster)
}

# Everything above len.short goes to mafft as a 'long' locus in the next step,
# and that is exactly the set worth a sequence check here. Tying the two
# together keeps the boundary meaningful instead of an arbitrary round number:
# on chr1 a threshold of 300 skipped 789 of 1921 candidate pairs, and 43 of
# those turned out to be real merges spread evenly over 51..299 bp.
len.min.merge = ifelse(is.null(opt$len.min.merge), len.short, opt$len.min.merge)
gap.max.merge = opt$gap.max.merge
# Nothing used to bound the span of a merge here, so a chain of iterations could
# undo what solveLong2 achieved in comb_04 and hand mafft an arbitrarily long locus.
len.max.merge = opt$len.max.merge
# Same distance rule mergeOverlapsTolerance() uses in comb_04: gap.max.merge alone
# lets two short breaks 5 kb apart be joined, which the relative test rejects.
dist.tol      = opt$dist.tol
len.max.sv    = opt$len.max.sv
p.ident.merge = opt$p.ident.merge
cover.merge   = opt$cover.merge
member.tol    = opt$member.tol

# Path with the MSA output (features)
path.features.msa <- opt$path.features.msa
path.inter.msa <- opt$path.inter.msa

if (is.null(path.features.msa) || is.null(path.inter.msa)) {
  stop("Error: both --path.features.msa and --path.inter.msa must be provided")
}

if (!dir.exists(path.features.msa)) stop('Features MSA directory does not exist')
if (!dir.exists(path.inter.msa)) stop('Internal MSA directory does not exist')

# Path to chromosomes
if (is.null(opt$path.chromosomes)) {
  stop("Error: --path.chromosomes must be provided")
}
path.chromosomes <- opt$path.chromosomes

# Every path below is glued with paste0, so a missing trailing slash would turn
# <dir> + 'breaks_ws_' into a sibling file name instead of a file inside <dir>.
addSlash <- function(s.path) {
  if (substr(s.path, nchar(s.path), nchar(s.path)) != '/') s.path <- paste0(s.path, '/')
  return(s.path)
}
path.features.msa <- addSlash(path.features.msa)
path.inter.msa    <- addSlash(path.inter.msa)
path.chromosomes  <- addSlash(path.chromosomes)
path.log          <- addSlash(path.log)

# ***************************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type.in, ".*")
files <- list.files(path = path.features.msa, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.in, "", files)
pref.combinations <- sub("\\.h5$", "", pref.combinations)

if(length(pref.combinations) == 0) {
  stop('No files with the ref-based alignments are found')
}

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ***************************************************************************
# ---- Functions ----

# Length of every accession's own piece inside a break.
getLenAcc <- function(v.beg, v.end) {
  len.acc <- abs(v.end) - abs(v.beg) - 1
  len.acc[(v.beg == 0) | (v.end == 0)] <- 0
  len.acc[len.acc < 0] <- 0
  return(len.acc)
}

# Accessions that substantially occupy a break. An accession contributing a
# couple of nucleotides while all of its sequence sits in the neighbouring break
# is alignment noise, not a second event, so it does not count as a member.
getMembers <- function(len.row, accessions, member.tol) {
  len.max <- max(len.row)
  if (len.max == 0) return(character(0))
  return(accessions[len.row >= member.tol * len.max])
}

# The sequence of a break in the accession that carries the longest piece.
getSeqOfBreak <- function(i.break, acc, v.beg, v.end, genome) {
  p1 <- abs(v.beg[i.break, acc]) + 1
  p2 <- abs(v.end[i.break, acc]) - 1
  if (p2 < p1) return('')
  s <- substring(genome, p1, p2)
  if (v.beg[i.break, acc] < 0) {
    s <- stringi::stri_reverse(chartr("ACGTRYSWKMBDHVNXUacgtryswkmbdhvnxu",
                                      "TGCAYRSWMKVHDBNNATGCAYRSWMKVHDBNNA", s))
  }
  return(as.character(s))
}

# Do two sequences look like the same thing? One blastn per pair is what makes
# this step expensive, so the pairs are handed out in chunks and a worker keeps
# reusing the same two temporary files instead of creating a pair per call.
# Every HSP is tested, not only the longest one: the longest hit is not
# necessarily the one that passes both thresholds.
isSameSeq <- function(s1, s2, f1, f2, p.ident.merge, cover.merge) {
  writeLines(c('>s1', s1), f1)
  writeLines(c('>s2', s2), f2)

  res <- system2('blastn',
                 args = c('-query', f1, '-subject', f2,
                          '-outfmt', shQuote('6 pident length'),
                          '-perc_identity', p.ident.merge),
                 stdout = TRUE, stderr = FALSE)

  if (length(res) == 0) return(FALSE)
  hits <- do.call(rbind, strsplit(res, '\t'))
  hits <- matrix(as.numeric(hits), ncol = 2)
  len.min <- min(nchar(s1), nchar(s2))
  return(any(hits[, 2] >= cover.merge * len.min))
}

# ***********************************************************************
# ---- MAIN program body ----

# Rebuild the set of completed combinations from all worker logs (any core count)
done.set <- getDoneSet(path.log)
if(length(done.set) > 0){
  pokaz('Skip already done:', length(done.set), file=file.log.main, echo=echo.main)
}
assign('.worker.id', 1, envir = .GlobalEnv)   # sequential outer loop -> single core_1.log

for(s.comb in pref.combinations){

  # ---- Checkpoint: skip already completed combinations ----
  s.comb.id <- s.comb
  if(s.comb.id %in% done.set) next

  # One log file per worker (bounded number of files); also the checkpoint ledger
  file.log.loop = initLoopLog(path.log)

  pokaz('* Combination', s.comb, file=file.log.loop, echo=echo.loop)
  q.chr = strsplit(s.comb, '_')[[1]][1]

  # Load variables from the previous step
  file.ws.pre = paste0(path.inter.msa, 'breaks_ws_pre_', s.comb, '.RData')
  file.ws     = paste0(path.inter.msa, 'breaks_ws_', s.comb, '.RData')
  if(!file.exists(file.ws.pre)){
    # NOT marked as done on purpose: a missing input is a problem to look at, and
    # marking it would make the combination unreachable on the next run.
    pokazAttention('File', file.ws.pre, 'does not exist, skip')
    next
  }
  load(file.ws.pre)

  # Merging changes which breaks are neighbours, so a chain A-B-C needs more than
  # one pass: the first pass joins A-B, the second sees the new A+B next to C.
  # Iterating here instead of relying on a re-run of the pipeline, which would
  # never come -- the combination is marked done at the end of this loop.
  n.iter.max <- 10
  n.merged.total <- 0

for(i.iter in 1:n.iter.max){

  # Only 'long' breaks are worth the sequence check: the short ones are indels
  # that abpoa aligns anyway, the singletons have nothing to compare.
  is.long <- (breaks$single != 1) & (breaks$len.acc > len.short)

  idx.ord <- order(breaks$idx.beg)
  idx.long <- idx.ord[is.long[idx.ord]]

  pokaz('Long breaks', length(idx.long), 'of', nrow(breaks), file=file.log.loop, echo=echo.loop)

  # Lengths only for the rows that can take part: the full breaks x accessions
  # matrix is never needed.
  len.acc.long <- getLenAcc(v.beg[idx.long, , drop=FALSE],
                            v.end[idx.long, , drop=FALSE])
  len.max.long <- apply(len.acc.long, 1, max)

  # ---- Candidate pairs: neighbouring, close, and with disjoint members ----
  # Disjoint membership on its own means nothing (most neighbouring breaks are
  # disjoint); it only marks the pairs whose sequences are worth comparing.
  n.long <- length(idx.long)
  is.cand <- rep(FALSE, max(n.long - 1, 0))
  if(n.long > 1){
    for(t in 1:(n.long-1)){
      i.br <- idx.long[t]
      j.br <- idx.long[t+1]

      # The distance rule is the one mergeOverlapsTolerance() applies in comb_04:
      # an absolute cap on the gap, a relative cap on the gap against the size of
      # what is being joined, and no merging of breaks that are already large.
      s.gap <- breaks$idx.beg[j.br] - breaks$idx.end[i.br] + 1
      if((s.gap < 0) || (s.gap > gap.max.merge)) next
      if(min(len.max.long[t], len.max.long[t+1]) < len.min.merge) next
      if(max(len.max.long[t], len.max.long[t+1]) > len.max.sv) next
      if((s.gap / max(len.max.long[t], len.max.long[t+1])) > dist.tol) next

      # comb_04 gets its span bound from solveLong2; this step runs after it and
      # iterates, so without this the chain could undo that bound.
      if((breaks$idx.end[j.br] - breaks$idx.beg[i.br] + 1) > len.max.merge) next

      acc.i <- getMembers(len.acc.long[t,],   accessions, member.tol)
      acc.j <- getMembers(len.acc.long[t+1,], accessions, member.tol)
      if(length(acc.i) == 0 || length(acc.j) == 0) next
      if(length(intersect(acc.i, acc.j)) > 0) next

      is.cand[t] <- TRUE
    }
  }

  t.cand <- which(is.cand)
  cand.i <- idx.long[t.cand]
  cand.j <- idx.long[t.cand + 1]

  pokaz('Iteration', i.iter, ': candidate pairs', length(cand.i),
        file=file.log.loop, echo=echo.loop)
  if(length(cand.i) == 0) break

  # ---- Compare the sequences of the candidates ----
  idx.merge <- integer(0)
  if(length(cand.i) > 0){

    # The representative of a break is the accession with the longest piece;
    # reading one chromosome per accession keeps this to one pass over the disk.
    acc.i <- accessions[apply(len.acc.long[t.cand, , drop=FALSE],   1, which.max)]
    acc.j <- accessions[apply(len.acc.long[t.cand+1, , drop=FALSE], 1, which.max)]

    seq.i <- rep('', length(cand.i))
    seq.j <- rep('', length(cand.j))

    for(acc in unique(c(acc.i, acc.j))){
      file.chromosome = paste(path.chromosomes, acc, '_chr', q.chr, '.fasta', sep = '')
      if(!file.exists(file.chromosome)){
        pokaz('File', file.chromosome, 'does not exist')
        stop()
      }
      genome = readFasta(file.chromosome)

      for(t in which(acc.i == acc)) seq.i[t] <- getSeqOfBreak(cand.i[t], acc, v.beg, v.end, genome)
      for(t in which(acc.j == acc)) seq.j[t] <- getSeqOfBreak(cand.j[t], acc, v.beg, v.end, genome)

      rmSafe(genome)
    }

    # ---- The blastn calls, spread over the cores ----
    # This is the whole cost of the step, so it is the thing that gets the cores.
    idx.test <- which((nchar(seq.i) > 0) & (nchar(seq.j) > 0))
    if(length(idx.test) > 0){
      n.chunk <- min(num.cores, length(idx.test))
      chunks <- split(idx.test, cut(seq_along(idx.test), n.chunk, labels = FALSE))

      if(num.cores > 1){
        # No .export for isSameSeq: foreach already picks up the globals the
        # expression uses, and naming it again only produces a warning per call.
        res.chunks <- foreach(i.chunk = seq_along(chunks),
                              .combine = c) %dopar% {
          f1 <- paste0(path.log, 'tmp_', s.comb, '_', i.chunk, '_1.fasta')
          f2 <- paste0(path.log, 'tmp_', s.comb, '_', i.chunk, '_2.fasta')
          ok <- integer(0)
          for(t in chunks[[i.chunk]]){
            if(isSameSeq(seq.i[t], seq.j[t], f1, f2, p.ident.merge, cover.merge)) ok <- c(ok, t)
          }
          if(file.exists(f1)) file.remove(f1)
          if(file.exists(f2)) file.remove(f2)
          ok
        }
        idx.merge <- sort(res.chunks)
      } else {
        f1 <- paste0(path.log, 'tmp_', s.comb, '_1.fasta')
        f2 <- paste0(path.log, 'tmp_', s.comb, '_2.fasta')
        for(t in idx.test){
          if(isSameSeq(seq.i[t], seq.j[t], f1, f2, p.ident.merge, cover.merge)) idx.merge <- c(idx.merge, t)
        }
        if(file.exists(f1)) file.remove(f1)
        if(file.exists(f2)) file.remove(f2)
      }
    }
  }

  if(length(idx.merge) == 0) break

  # ---- Merge the pairs ----
  # The two breaks become one: the interval spans both, and every accession keeps
  # the longest of the pieces available inside that span.
  if(length(idx.merge) > 0){

    i.keep <- cand.i[idx.merge]
    i.drop <- cand.j[idx.merge]

    # Candidates are neighbours among the LONG breaks, so short breaks and
    # singletons may sit between them. They are swallowed by the merged interval
    # and must be merged in as well -- leaving them as separate rows would nest a
    # break inside another one, which the coordinate map of comb_10 cannot
    # represent (it adds every break's `extra` at a single point and takes a
    # cumsum, so a nested break shifts its container's end but not its beginning).
    mid.list <- vector('list', length(i.keep))
    for(k in seq_along(i.keep)){
      mid.list[[k]] <- which((breaks$idx.beg > breaks$idx.end[i.keep[k]]) &
                             (breaks$idx.end < breaks$idx.beg[i.drop[k]]))
    }

    # A break can belong to two pairs at once (A-B and B-C). Take them greedily
    # within this pass so that no break is merged twice; what is left over is
    # picked up by the next iteration, when A+B and C are neighbours. The
    # swallowed breaks take part in the bookkeeping on equal terms.
    idx.used <- integer(0)
    idx.take.pair <- rep(FALSE, length(i.keep))
    for(k in seq_along(i.keep)){
      grp <- c(i.keep[k], mid.list[[k]], i.drop[k])
      if(any(grp %in% idx.used)) next
      idx.take.pair[k] <- TRUE
      idx.used <- c(idx.used, grp)
    }
    if(sum(!idx.take.pair) > 0){
      pokaz('Chained pairs left for the next iteration', sum(!idx.take.pair),
            file=file.log.loop, echo=echo.loop)
    }
    i.keep   <- i.keep[idx.take.pair]
    i.drop   <- i.drop[idx.take.pair]
    mid.list <- mid.list[idx.take.pair]

    n.swallowed <- sum(lengths(mid.list))
    if(n.swallowed > 0){
      pokaz('Breaks swallowed by the merged intervals', n.swallowed,
            file=file.log.loop, echo=echo.loop)
    }

    rows.drop <- integer(0)
    for(t in seq_along(i.keep)){
      i.br <- i.keep[t]
      j.br <- i.drop[t]
      grp  <- c(i.br, mid.list[[t]], j.br)

      # Positions: whoever has the longest piece gives it. Membership was checked
      # with a tolerance, so an accession may hold a few noise nucleotides in one
      # break and its real sequence in another; comparing lengths keeps the real
      # one. With `i` first, ties keep the same row the pairwise version kept.
      len.grp <- getLenAcc(v.beg[grp, , drop=FALSE], v.end[grp, , drop=FALSE])
      sel     <- grp[max.col(t(len.grp), ties.method = "first")]
      i.acc   <- seq_along(accessions)

      v.beg[i.br, ] <- v.beg[cbind(sel, i.acc)]
      v.end[i.br, ] <- v.end[cbind(sel, i.acc)]

      len.new <- getLenAcc(v.beg[i.br, , drop=FALSE], v.end[i.br, , drop=FALSE])

      breaks$idx.end[i.br]  <- breaks$idx.end[j.br]
      breaks$len.acc[i.br]  <- max(len.new)
      breaks$cnt[i.br]      <- sum(breaks$cnt[grp])
      breaks$len[i.br]      <- breaks$idx.end[i.br] - breaks$idx.beg[i.br] + 1
      breaks$single[i.br]   <- sum(len.new > 0)
      breaks$len.mean[i.br] <- mean(len.new)

      rows.drop <- c(rows.drop, mid.list[[t]], j.br)
    }

    breaks <- breaks[-rows.drop,]
    v.beg  <- v.beg[-rows.drop, , drop=FALSE]
    v.end  <- v.end[-rows.drop, , drop=FALSE]

    n.merged.total <- n.merged.total + length(i.keep)
    pokaz('Iteration', i.iter, ': merged', length(i.keep),
          ', breaks left', nrow(breaks), file=file.log.loop, echo=echo.loop)
  }

}   # end of the merging iterations

  pokaz('Merged in total', n.merged.total, file=file.log.loop, echo=echo.loop)

  # Nested/overlapping breaks are unrepresentable downstream: comb_10 would only
  # notice three steps later, as 'Wrong mapping: length mismatch'. Same guard
  # mergeOverlapsTolerance() already has.
  ord <- order(breaks$idx.beg)
  if(nrow(breaks) > 1 &&
     any(breaks$idx.beg[ord][-1] <= breaks$idx.end[ord][-nrow(breaks)])){
    stop('comb_04b: overlapping breaks after merging')
  }

  # ---- Save the workspace for the next step ----
  # Always written, even when nothing was merged: breaks_ws_<comb>.RData is what
  # the next step reads, and it has to exist for every combination.
  all.local.objects <- c('breaks', 'v.beg', 'v.end', 'accessions')
  save(list = all.local.objects, file = file.ws)

  # ---- Checkpoint marker: combination fully processed ----
  markDone(s.comb.id, file=file.log.loop, echo=echo.loop)
}

if(num.cores > 1){
  stopCluster(myCluster)
}
