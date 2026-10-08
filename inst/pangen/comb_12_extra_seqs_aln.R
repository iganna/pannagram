# Align extra-long fragments with the step-16 aligner (refineAlignment).
#
# Same CLI as before, same contract with comb_11 (reads <id>_group.fasta / _group.txt /
# breaks_info_<comb>.RData) and with comb_13 (writes <id>_out.RData with idx.new/msa.new
# and <id>_len.RData with len.aln), so pannagram.sh needs no change.
#
# What changes: instead of adding the accessions' break sequences to the consensus ONE
# BY ONE (each add = one mafftAdd or one makeblastdb+blastn, plus a full cbind of the
# growing matrices), every accession's own sequence across the break window is assembled
# once and the whole locus goes through refineAlignment() - the same clustering +
# alignTwoLong/chainHsps merge that comb_08 (step 16) uses on bad-MAFFT loci. That
# removes the 15 kb fork, the per-accession BLAST, the order dependence and the
# quadratic bookkeeping.
#
# Operational fixes ported from comb_08_mafft2.R:
#   - largest loci first (a big break starting last would hold one core to the end)
#   - per-break scratch cleanup (comb_12 left every tempfile in one flat tmp/ dir)
#   - gc(reset) + peak-memory logging per break
#   - breaks that are skipped are RECORDED (skipped.tsv), not dropped silently

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
  library(pannagram)
  library(igraph)
  library(R.utils)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func.R", package = "pannagram"))
source(system.file("pangen/synteny_func.R", package = "pannagram"))
source(system.file("pangen/comb_func_mafft_refine2.R", package = "pannagram"))

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.cons"),        type = "character", default = NULL, help = "path to consensus directory"),
  make_option(c("--path.chromosomes"), type = "character", default = NULL, help = "path to directory with chromosomes"),
  make_option(c("--path.extra"),       type = "character", default = NULL, help = "path to directory with the break files"),
  make_option(c("-c", "--cores"),      type = "integer",   default = 1,    help = "number of cores to use for parallel processing"),
  make_option(c("--path.log"),         type = "character", default = NULL, help = "Path for log files"),
  make_option(c("--log.level"),        type = "character", default = NULL, help = "Level of log to be shown on the screen"),
  make_option(c("--aln.type.in"),      type = "character", default = NULL, help = "Alignment type"),
  make_option(c("--path.work"),        type = "character", default = NULL, help = "scratch dir for mafft/blast (default: <path.extra>/tmp/)"),
  make_option(c("--max.total.bp"),     type = "double",    default = 5e6,  help = "skip a break whose total input exceeds this")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram"))

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram"))

# ***********************************************************************

aln.type.in <- opt$aln.type.in
if(is.null(aln.type.in)) aln.type.in = aln.type.add

num.cores = opt$cores
if(is.null(num.cores)) stop('Wrong number of cores: NULL')
pokaz('Number of cores', num.cores)

if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop('Consensus folder does not exist')
if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes
if (!is.null(opt$path.extra)) path.extra <- opt$path.extra

# Scratch for mafft/blast. Kept separate from path.extra so it can be put on a local
# disk: this step creates and removes a few files per merge, and on GPFS that metadata
# traffic - not the alignment - is what the step spends its time on.
path.work = if(is.null(opt$path.work)) paste0(path.extra, 'tmp/') else opt$path.work
if(!endsWith(path.work, '/')) path.work = paste0(path.work, '/')
if (!dir.exists(path.work)) dir.create(path.work, recursive = TRUE)

max.total.bp = opt$max.total.bp
file.skip = paste0(path.log, 'skipped.tsv')

# **************************************************************************
# ---- Combinations ----

aln.type.in = paste0(aln.type.in, '_')
s.pattern <- paste0("^", aln.type.in, "\\d+_\\d+\\.h5$")
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.in, "", files)
pref.combinations <- sub("\\.h5$", "", pref.combinations)
if(length(pref.combinations) == 0) stop('No files alignments are found')

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

echo = T

done.set <- getDoneSet(path.log)
if(length(done.set) > 0){
  pokaz('Skip already done:', length(done.set), file=file.log.main, echo=echo.main)
}
assign('.worker.id', 1, envir = .GlobalEnv)

# Strip the "acc_<ACC>_break_..." wrapper that comb_11 puts on every sequence name
accFromName <- function(s){
  sub("^acc_", "", sapply(strsplit(s, '_break_'), `[`, 1))
}

for(s.comb in pref.combinations){

  if(echo) pokaz('* Combination', s.comb)

  file.breaks.info = paste0(path.extra, "breaks_info_",s.comb,".RData")
  if(!file.exists(file.breaks.info)){
    if(echo) pokaz('No breaks for combination', s.comb, '- nothing to align')
    next
  }
  load(file.breaks.info)  # "breaks.init", "breaks"

  breaks$len.true = breaks$idx.end - breaks$idx.beg - 1

  file.comb = paste0(path.cons, aln.type.in, s.comb,'.h5')
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]

  # ---- Cost proxy, largest first ----------------------------------------------
  # The pool hands out breaks in list order; in index order a megabase break can start
  # last and keep one core busy long after the rest of the chromosome is done.
  add.bp = tapply(breaks.init$len.acc, factor(breaks.init$gr, levels = seq_len(nrow(breaks))), sum)
  add.bp[is.na(add.bp)] = 0
  breaks$cost = breaks$len.true * length(accessions) + as.numeric(add.bp)
  order.b = order(-breaks$cost)

  # ---- Additional alignments ----

  loop.function <- function(i.b, breaks, done.set = character(0), echo.loop=T){

    item.id <- breaks$id.s[i.b]
    if(item.id %in% done.set) return(NULL)

    file.log.loop = initLoopLog(path.log)

    file.br.group = paste0(path.extra, item.id, '_group.fasta')
    file.br.idx   = paste0(path.extra, item.id, '_group.txt')
    file.br.out   = paste0(path.extra, item.id, '_out.RData')
    file.br.len   = paste0(path.extra, item.id, '_len.RData')

    pos.b <- breaks$idx.beg[i.b]
    pos.e <- breaks$idx.end[i.b]
    len.true <- breaks$len.true[i.b]

    # Individual breaks that live inside this merged group
    breaks.tmp <- breaks.init[(breaks.init$idx.beg >= pos.b) & (breaks.init$idx.end <= pos.e),,drop=F]

    # ---- The already aligned interior of the group (may be absent) ----
    has.group = file.exists(file.br.idx)
    if(has.group){
      indexes = read.table(file.br.idx, row.names = 1, header = F)
      rownames(indexes) = accFromName(rownames(indexes))
      aln = aln2mx(readFasta(file.br.group))
      rownames(aln) = accFromName(rownames(aln))
      if(ncol(aln) != len.true)  stop('Group alignment width does not match the break')
      if(ncol(indexes) != len.true) stop('Group index width does not match the break')
    }

    # ---- Assemble, per accession, its own contiguous sequence over the window ----
    # Column j of the group is pangenome column pos.b + j. A break sits immediately
    # after its left flank column idx.beg, so ordering by (major = pangenome column,
    # minor = offset inside the break) puts every base in its own-coordinate order.
    seqs = character(0)
    own.list = list()
    for(acc in accessions){
      maj = integer(0); mnr = integer(0); own = numeric(0); nts = character(0)

      if(has.group && (acc %in% rownames(indexes))){
        vb = as.numeric(indexes[acc,])
        ab = aln[acc,]
        nz = which(vb != 0)
        if(length(nz) > 0){
          if(any(ab[nz] == '-')) stop('Aligned position carries a gap')
          maj = c(maj, pos.b + nz); mnr = c(mnr, rep(0L, length(nz)))
          own = c(own, vb[nz]);     nts = c(nts, ab[nz])
        }
      }

      br = breaks.tmp[breaks.tmp$acc == acc,,drop=F]
      if(nrow(br) > 0){
        for(j in seq_len(nrow(br))){
          poses = (br$val.beg[j]+1):(br$val.end[j]-1)
          s.j = seq2nt(br$seq[j])
          if(length(s.j) != length(poses)) stop('Mismatch sequence and positions')
          maj = c(maj, rep(br$idx.beg[j], length(poses)))
          mnr = c(mnr, seq_along(poses))
          own = c(own, poses); nts = c(nts, s.j)
        }
      }

      if(length(own) == 0) next
      o = order(maj, mnr)
      # as.numeric: an accession whose sequence comes only from breaks would keep the
      # integer type of the ':' sequence, and the checks below compare against a double
      # matrix.
      own.list[[acc]] = as.numeric(own[o])
      seqs[acc] = nt2seq(toupper(nts[o]))
    }

    if(length(seqs) == 0) stop('No sequences in the break')

    total.bp = sum(nchar(seqs))
    if(total.bp > max.total.bp){
      # Recorded, not silent: comb_13 drops breaks without _out/_len, and until now
      # nothing said which ones or why.
      cat(paste(item.id, length(seqs), total.bp, 'too_big', sep='\t'), '\n',
          sep='', file=file.skip, append=TRUE)
      pokazAttention('Break too big, skipped:', item.id, total.bp,
                     file=file.log.loop, echo=echo.loop)
      return(NULL)
    }

    # ---- Align ----
    if(length(seqs) == 1){
      acc1 = names(seqs)[1]
      aln.mx = matrix(seq2nt(seqs[1]), nrow = 1, dimnames = list(acc1, NULL))
      pos.mx = matrix(seq_len(nchar(seqs[1])), nrow = 1, dimnames = list(acc1, NULL))
    } else {
      prefix.work = paste0(item.id, '_')
      path.work.pref = paste0(path.work, prefix.work)
      on.exit({
        f.tmp = list.files(path.work)
        f.tmp = f.tmp[startsWith(f.tmp, prefix.work)]
        if(length(f.tmp) > 0) unlink(file.path(path.work, f.tmp))
      }, add = TRUE)

      res = refineAlignment(seqs, path.work.pref, keep.pos = TRUE)
      aln.mx = res$aln[[length(res$aln)]]
      pos.mx = res$pos[[length(res$pos)]]
      rm(res)
    }

    if(!setequal(rownames(aln.mx), names(seqs))) stop('Aligner lost a sequence')

    # ---- Back to the comb_13 contract ----
    msa.new = matrix('-', nrow = length(accessions), ncol = ncol(aln.mx),
                     dimnames = list(accessions, NULL))
    idx.new = matrix(0, nrow = length(accessions), ncol = ncol(aln.mx),
                     dimnames = list(accessions, NULL))

    for(acc in rownames(aln.mx)){
      msa.new[acc, ] = aln.mx[acc, ]
      pp = pos.mx[acc, ]
      nzc = which(pp != 0)
      idx.new[acc, nzc] = own.list[[acc]][pp[nzc]]
      # Every base of the input must come back exactly once, and with the same own
      # coordinate: comb_13 checks uniqueness genome-wide and would only fail much later.
      got = sort(abs(idx.new[acc, nzc])); want = sort(abs(own.list[[acc]]))
      if(length(got) != length(want) || any(got != want)){
        stop(paste('Positions of', acc, 'were not preserved'))
      }
    }

    # Empty columns (all gaps) - the same cleanup the original did
    idx.zero = which(colSums(idx.new != 0) == 0)
    if(length(idx.zero) > 0){
      idx.gap = which(colSums(msa.new != '-') == 0)
      if(!setequal(idx.zero, idx.gap)) stop('Wrong zeros')
      idx.new = idx.new[, -idx.zero, drop = F]
      msa.new = msa.new[, -idx.gap, drop = F]
    }

    # comb_13 recomputes every downstream coordinate from len.aln - len.true and stops
    # on a negative gain. A fresh alignment can be NARROWER than the block it replaces
    # (the old one had gap-only columns), so pad with empty columns rather than lose
    # the break.
    if(ncol(idx.new) < len.true){
      n.pad = len.true - ncol(idx.new)
      idx.new = cbind(idx.new, matrix(0, nrow = nrow(idx.new), ncol = n.pad,
                                      dimnames = list(accessions, NULL)))
      msa.new = cbind(msa.new, matrix('-', nrow = nrow(msa.new), ncol = n.pad,
                                      dimnames = list(accessions, NULL)))
    }

    if(ncol(idx.new) != ncol(msa.new)) stop('Lengths are different')
    len.aln = ncol(idx.new)
    if(len.aln < len.true) stop(paste('The length has decreased', i.b))

    save(list = c("len.aln"), file = file.br.len)
    save(list = c("idx.new", "msa.new"), file = file.br.out)

    mem = gc(reset = TRUE)
    pokaz('break', item.id, 'n.seq', length(seqs), 'in.bp', total.bp,
          'cols', len.aln, 'mem.peak.mb', round(sum(mem[, 6])),
          file=file.log.loop, echo=echo.loop)

    markDone(item.id, file=file.log.loop, echo=echo.loop)
    return(NULL)
  }

  if(num.cores == 1){
    for(i.b in order.b){
      tryCatch({
        loop.function(i.b, breaks, done.set = done.set, echo.loop = echo.loop)
      },
      error = function(err) {
        cat(paste(breaks$id.s[i.b], NA, NA, paste0('error:', err$message), sep='\t'), '\n',
            sep='', file=file.skip, append=TRUE)
        pokazAttention("An error occurred for i.b =", i.b, ":", err$message)
      })
    }
  } else {
    myCluster <- makeCluster(num.cores, type = "PSOCK")
    registerDoParallel(myCluster)
    parallel::clusterApply(myCluster, seq_len(num.cores),
                           function(i) assign('.worker.id', i, envir = .GlobalEnv))

    tmp = foreach(i.b = order.b,
                  .packages=c('rhdf5', 'crayon', 'igraph', 'pannagram', 'R.utils'))  %dopar% {
                    tryCatch({
                      loop.function(i.b, breaks, done.set = done.set, echo.loop = echo.loop)
                    },
                    error = function(err) {
                      cat(paste(breaks$id.s[i.b], NA, NA, paste0('error:', err$message), sep='\t'), '\n',
                          sep='', file=file.skip, append=TRUE)
                      pokazAttention("An error occurred for i.b =", i.b, ":", err$message)
                      return(NULL)
                    })
                  }
    stopCluster(myCluster)
  }

  H5close()
  gc()
}

pokaz('Done.', file=file.log.main, echo=echo.main)
