# Combine all alignments together into the final one

suppressMessages({
library('foreach')
library(doParallel)
library("optparse")
source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func_mafft_refine2.R", package = "pannagram"))
source(system.file("pangen/synteny_func.R", package = "pannagram"))
})

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.mafft.in"), type="character", default=NULL, 
              help="path to directory, where to combine fasta files for mafft runs", metavar="character"),
  make_option(c("--path.mafft.out"), type="character", default=NULL, 
              help="path to directory, where to mafft results are", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("--path.log"), type = "character", default = NULL,
              help = "Path for log files", metavar = "character"),
  make_option(c("--log.level"), type = "character", default = NULL,
              help = "Level of log to be shown on the screen", metavar = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging


# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing.
# No per-step cap: the memory of this step is bounded per LOCUS (see the notes in
# refineAlignment() about releasing merged parents and about keep.pos), not by
# starving the cores. A capped core count only lowers the probability of an OOM,
# it never bounds it - the peak depends on which loci happen to be in flight.
num.cores = opt$cores

if (!is.null(opt$path.mafft.in)) path.mafft.in <- opt$path.mafft.in
if (!is.null(opt$path.mafft.out)) path.mafft.out <- opt$path.mafft.out

# ***********************************************************************
# ---- Preparation ----

n.flank = 30

files.extra <- list.files(path = path.mafft.in, pattern = "\\.fasta$", full.names = F)
files.extra = files.extra[!grepl('_aligned', files.extra)]

if(length(files.extra) == 0){
  pokaz('Number of files for extra alignment')
  quit(status = 0)
} else {
  pokaz('Number of files to align', length(files.extra))
}

path.mafft.in.tmp = paste0(path.mafft.in, 'tmp/')

if (!file.exists(path.mafft.in.tmp)) {
  dir.create(path.mafft.in.tmp)
}

# ***********************************************************************
# ---- MAIN program body ----

loop.function <- function(f.in,
                          done.set = character(0),
                          echo.loop=T){
  pokaz(f.in)

  # ---- Checkpoint: skip already completed items ----
  item.id <- sub("\\.[^.]*$", "", basename(f.in))
  if(item.id %in% done.set){
    return(NULL)
  }

  # One log file per worker (bounded number of files); also the checkpoint ledger
  file.log.loop = initLoopLog(path.log)
  
  seqs = readFasta(paste0(path.mafft.in, f.in))
  seqs.clean = seq2clean(seqs, n.flank)
  
  seqs.clean = seqs.clean[nchar(seqs.clean) > 7]
  if(length(seqs.clean) < 2) {
    pokazAttention('Not enough sequences to align', file=file.log.loop, echo=echo.loop)
    markDone(item.id, file=file.log.loop, echo=echo.loop)
    return()
  }
  
  # Proportion of non-N nucleotides: drop N-heavy sequences INDIVIDUALLY
  # (a single N-rich sequence must not discard the whole locus).
  # Counted on the string: seq2nt() would expand every sequence of the locus into a
  # character vector (8 bytes per base against 1 in the string) just to be discarded.
  n.n = 1 - (nchar(seqs.clean) - nchar(gsub('[Nn]', '', seqs.clean))) / nchar(seqs.clean)
  seqs.clean = seqs.clean[n.n >= 0.5]

  if(length(seqs.clean) < 2){
    pokazAttention('Too many N to align', file=file.log.loop, echo=echo.loop)
    markDone(item.id, file=file.log.loop, echo=echo.loop)
    return()
  }
  
  path.work = paste0(path.mafft.in.tmp, sub('\\.fasta', '', basename(f.in)), '_')
  pokaz(path.work)

  # Drop the scratch of this locus as soon as it is finished, and also when it fails.
  # refineAlignmentRouted() leaves ~17 files per locus (seqs_*, aln_*, seg*, tbl, blast)
  # in the single per-chromosome tmp/ directory; kept until the end of the chromosome
  # they pile up to tens of thousands (53908 on a 3144-locus chromosome), and the
  # metadata load of `cores` workers creating and removing files in one flat directory
  # produces transient I/O failures. With per-locus cleanup the peak is ~17 * cores.
  prefix.work = basename(path.work)
  on.exit({
    files.tmp = list.files(path.mafft.in.tmp)
    files.tmp = files.tmp[startsWith(files.tmp, prefix.work)]   # the trailing '_' of the
    if(length(files.tmp) > 0){                                  # prefix excludes locus_1 vs locus_10
      unlink(file.path(path.mafft.in.tmp, files.tmp))
    }
  }, add = TRUE)
  # keep.pos = FALSE: only the final alignment is written out, so the parallel list of
  # position matrices (the same size as the alignment matrices) is never built.
  res = refineAlignmentRouted(seqs.clean, path.work, keep.pos = FALSE)
  
  alignment = res$aln[[length(res$aln)]]
  rm(res)
  
  alignment.seq = mx2aln(alignment)
  rm(alignment)
  
  file.out = paste0(path.mafft.out, sub('\\.fasta', '', basename(f.in)), "_aligned.fasta")
  writeFasta(alignment.seq, file.out)
  rm(seqs, seqs.clean, alignment.seq)

  # Release the locus back to the OS before taking the next one: a PSOCK worker lives
  # for the whole step, so without this its heap only ever grows to the largest locus
  # it happened to see. gc(reset) also reports the peak, logged next to the input size
  # so that the real cost of a locus can be read off the worker logs.
  mem = gc(reset = TRUE)
  pokaz('mem.peak.mb', round(sum(mem[, 6])), 'input.mb',
        round(file.size(paste0(path.mafft.in, f.in)) / 2^20, 1),
        file=file.log.loop, echo=echo.loop)

  # ---- Checkpoint marker: item fully processed ----
  markDone(item.id, file=file.log.loop, echo=echo.loop)

}


# ***********************************************************************
# ---- Loop  ----

# Rebuild the set of completed items from all worker logs (any core count)
done.set <- getDoneSet(path.log)
if(length(done.set) > 0){
  pokaz('Skip already done:', length(done.set), file=file.log.main, echo=echo.main)
}

if(num.cores == 1){
  assign('.worker.id', 1, envir = .GlobalEnv)   # stable single file core_1.log
  for(f.in in files.extra){
    loop.function(f.in,
                  done.set = done.set,
                  echo.loop=echo.loop)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK")
  registerDoParallel(myCluster)

  # Assign a stable worker id (1..N) to each worker -> bounded, reused file names
  parallel::clusterApply(myCluster, seq_len(num.cores),
                         function(i) assign('.worker.id', i, envir = .GlobalEnv))

  tmp = foreach(f.in = files.extra,
                .packages=c('crayon'),
                .verbose = F)  %dopar% {
                  loop.function(f.in,
                                done.set = done.set,
                                echo.loop=echo.loop)
                }
  stopCluster(myCluster)
}


pokaz('Done.', file=file.log.main, echo=echo.main)

