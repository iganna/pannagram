# Exact boundaries of the consensus sequences of SV families by their copies.
#
# The consensus of sv_08_family_cons.R is built from the SVs themselves, and the
# boundaries of an SV do not have to be the boundaries of the element. So every
# consensus is searched in the pangenome consensus by simsearch, the full-length
# copies are taken with short flanks and aligned, and the element is the block
# where the copies agree: every copy sits in its own place of the genome, so the
# flanks do not agree. Optionally, the boundaries are shifted to the strongest
# target site duplication shared by the copies.
#
# In the end every consensus gets the measures of its quality: how many copies
# support it, how close they are to it, and how sharp its boundaries are.

suppressMessages({
  library(optparse)
  library(pannagram)
  library(parallel)
  library(data.table)
})

args = commandArgs(trailingOnly = TRUE)

option_list = list(
  make_option("--path.out",     type = "character", default = NULL, help = "Path for the results, the same as of sv_08_family_cons.R"),
  make_option("--file.cons",    type = "character", default = NULL, help = "Consensus sequences, default <path.out>/sv_families_cons_<sim>_<cov>.fasta"),
  make_option("--path.seq",     type = "character", default = NULL, help = "Directory with the pangenome consensus seq_cons_*.fasta"),
  make_option("--seq.suff",     type = "character", default = "",   help = "Suffix of the seq_cons files, e.g. _extra1"),
  make_option("--similarity",   type = "integer",   default = 85,   help = "Similarity of the families"),
  make_option("--coverage",     type = "integer",   default = 85,   help = "Coverage of the families"),
  make_option("--cores",        type = "integer",   default = 1,    help = "Number of cores"),
  make_option("--sim.copy",     type = "integer",   default = 85,   help = "Minimum similarity of a copy to the consensus, simsearch -sim"),
  make_option("--cov.copy",     type = "integer",   default = 90,   help = "Minimum coverage of the consensus by a full-length copy, simsearch -cov"),
  make_option("--max.len",      type = "double",    default = 3,    help = "Maximum length of a full-length copy relative to the consensus"),
  make_option("--min.copies",   type = "integer",   default = 5,    help = "Minimum number of full-length copies"),
  make_option("--n.copies",     type = "integer",   default = 20,   help = "Maximum number of copies to align"),
  make_option("--flank",        type = "integer",   default = 50,   help = "Length of the flanks"),
  make_option("--bnd.cov",      type = "integer",   default = 3,    help = "Minimum coverage of a column inside the element"),
  make_option("--bnd.cov.frac", type = "double",    default = 0.5,  help = "Minimum share of the copies covering a column inside the element"),
  make_option("--bnd.agree",    type = "double",    default = 0.75, help = "Minimum share of the dominant nucleotide in a column inside the element"),
  make_option("--bnd.wnd",      type = "integer",   default = 10,   help = "Window at the boundary"),
  make_option("--bnd.good",     type = "integer",   default = 8,    help = "Minimum number of good columns in the window"),
  make_option("--tsd",          type = "logical",   default = FALSE, help = "Shift the boundaries to the target site duplication"),
  make_option("--tsd.min",      type = "integer",   default = 3,    help = "Minimum length of the TSD"),
  make_option("--tsd.max",      type = "integer",   default = 20,   help = "Maximum length of the TSD"),
  make_option("--tsd.shift",    type = "integer",   default = 3,    help = "Maximum shift of a boundary to the TSD"),
  make_option("--q.wnd",        type = "integer",   default = 20,   help = "Window for the sharpness of the boundaries"),
  make_option("--q.ident",      type = "double",    default = 0.8,  help = "Minimum median identity of the copies for a good consensus"),
  make_option("--q.contrast",   type = "double",    default = 0.3,  help = "Minimum drop of the agreement across a sharp boundary"),
  make_option("--min.out",      type = "integer",   default = 20,   help = "Flank that has to remain outside of the element, otherwise the end is open"),
  make_option("--mafft.args",   type = "character", default = "--op 5 --ep 0.2 --quiet", help = "Arguments of MAFFT"),
  make_option("--path.code",    type = "character", default = "",   help = "Directory with the sources, if the package is not rebuilt yet")
)

opt = parse_args(OptionParser(option_list = option_list), args = args)

# The package may be not rebuilt yet: take the functions from the repository
if(nzchar(opt$path.code)){
  for(f in c('analys/sv_aln_func.R')){
    p = file.path(opt$path.code, f)
    if(file.exists(p)) source(p)
  }
}

path.out = opt$path.out
if(is.null(path.out)) stop('path.out is not provided')
suff = paste0(opt$similarity, '_', opt$coverage)
file.cons = if(is.null(opt$file.cons)) file.path(path.out, paste0('sv_families_cons_', suff, '.fasta')) else opt$file.cons
if(!file.exists(file.cons)) stop(paste0('No consensus sequences: ', file.cons))

if(is.null(opt$path.seq) || !dir.exists(opt$path.seq)) stop(paste0('No seq dir: ', opt$path.seq))
files.seq = list.files(opt$path.seq, pattern = paste0('^seq_cons_[0-9]+_[0-9]+', opt$seq.suff, '\\.fasta$'),
                       full.names = TRUE)
if(length(files.seq) == 0) stop(paste0('No seq_cons files in ', opt$path.seq))

path.ext = file.path(path.out, 'ext')
dir.create(path.ext, showWarnings = FALSE, recursive = TRUE)
path.sim = file.path(path.ext, 'simsearch')
dir.create(path.sim, showWarnings = FALSE, recursive = TRUE)

# ***********************************************************************
# ---- Copies ----

cons.old = readFasta(file.cons)
fam.ids = sub('^fam_', '', sub('\\|.*', '', names(cons.old)))
names(cons.old) = fam.ids
pokaz('Consensus sequences:', length(cons.old))

file.query = file.path(path.ext, 'query.fasta')
writeFasta(setNames(cons.old, paste0('fam_', fam.ids)), file.query)

suff.sim = paste0(opt$sim.copy, '_', opt$cov.copy)
copies = NULL
for(f in files.seq){
  pref = sub('\\..*', '', basename(f))
  file.tab = file.path(path.sim, paste0(pref, '_', suff.sim, '.table'))
  if(!file.exists(file.tab)){
    status = system(paste('simsearch -in_seq', shQuote(file.query), '-on_genome', shQuote(f),
                          '-out', shQuote(path.sim), '-sim', opt$sim.copy, '-cov', opt$cov.copy,
                          '-cores', opt$cores))
    if(status != 0) stop(paste0('simsearch has failed on ', f))
  }
  # no copies at all in this sequence: simsearch writes nothing
  if(!file.exists(file.tab)) next
  x = read.table(file.tab, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  x$file = f
  copies = rbind(copies, x)
}
if(is.null(copies)) stop('No copies found')
copies$family = sub('^fam_', '', copies$sequence)
pokaz('Copies:', nrow(copies))

# Full-length copies only: simsearch has already required the coverage of the
# consensus. A copy may still be much longer than the consensus: when the SVs
# carry only the ends of the element, every copy in the genome has the middle,
# and the alignment of the copies brings it back. Only a copy so long that it is
# rather two copies chained together is left out.
len.copy = copies$end.genome - copies$beg.genome + 1
copies = copies[len.copy <= opt$max.len * copies$len.seq, ]
copies = copies[order(copies$family, -copies$coverage, -copies$similarity), ]
n.full = table(factor(copies$family, levels = fam.ids))
copies = copies[ave(seq_len(nrow(copies)), copies$family, FUN = seq_along) <= opt$n.copies, ]
copies$name = paste0(copies$genome.chr, ':', copies$beg.genome, '-', copies$end.genome, ':', copies$strand)
pokaz('Full-length copies taken:', nrow(copies))

# The sequences are cut one chromosome at a time, so that the whole pangenome
# is never in memory.
pokaz('Cutting the copies with the flanks of', opt$flank, 'bp...')
seq.copy = character(nrow(copies))
for(f in unique(copies$file)){
  s.chr = readFasta(f)
  idx = which(copies$file == f)
  for(i in idx){
    s = s.chr[[copies$genome.chr[i]]]
    b = max(1, copies$beg.genome[i] - opt$flank)
    e = min(nchar(s), copies$end.genome[i] + opt$flank)
    x = substr(s, b, e)
    seq.copy[i] = if(copies$strand[i] == '-') revCompl(x) else x
  }
  rm(s.chr)
}
seq.copy = setNames(toupper(seq.copy), copies$name)
copies.by.fam = split(seq_len(nrow(copies)), copies$family)

# ***********************************************************************
# ---- One family ----

oneFamily <- function(f.id){

  s.cons = cons.old[[f.id]]
  res = data.frame(family = f.id, len.old = nchar(s.cons), len.cons = NA,
                   n.full = as.integer(n.full[[f.id]]), n.aln = 0,
                   open5 = NA, open3 = NA, status = '', stringsAsFactors = FALSE)
  cons = NULL

  tryCatch({
    idx = copies.by.fam[[f.id]]
    if(length(idx) < opt$min.copies) stop('few full-length copies')
    res$n.aln = length(idx)
    s.cp = seq.copy[idx]

    # ---- 1. alignment of the copies with their flanks ----
    mx = svAlnMafft(s.cp, path.work = tempdir(), args = opt$mafft.args)
    bn = svAlnBounds(mx, min.cov = opt$bnd.cov, min.cov.frac = opt$bnd.cov.frac, min.agree = opt$bnd.agree,
                     wnd = opt$bnd.wnd, wnd.good = opt$bnd.good)
    if(is.null(bn)) stop('no block of agreement')

    # an end is open when the block of agreement reaches the end of the flanks:
    # the element goes on further than the consensus has shown
    cov.ok = (bn$cov >= bn$min.cov)
    res$open5 = sum(cov.ok[seq_len(bn$beg - 1)]) < opt$min.out
    res$open3 = sum(cov.ok[seq_along(cov.ok) > bn$end]) < opt$min.out

    # ---- 2. boundaries ----
    cons.col = svAlnCons(mx)$cons
    col.nt = which(cons.col != '-')
    beg = bn$beg
    end = bn$end

    # A shift is measured in the nucleotides of the copies, and it is applied to
    # the nucleotides of the consensus: in the few columns around a boundary the two
    # do not differ.
    shiftCol <- function(col, shift){
      i = which.min(abs(col.nt - col))
      col.nt[min(max(i + shift, 1), length(col.nt))]
    }

    res$tsd.len = 0; res$tsd.seq = ''; res$tsd.n = 0; res$tsd.shift5 = 0; res$tsd.shift3 = 0
    if(opt$tsd){
      rb = svAlnRowBounds(mx, beg, end)
      tsd = svFindTSD(s.cp[rownames(mx)], rb$beg, rb$end, k.min = opt$tsd.min, k.max = opt$tsd.max,
                      max.shift = opt$tsd.shift)
      if(!is.null(tsd)){
        res$tsd.len = tsd$k; res$tsd.seq = tsd$tsd; res$tsd.n = tsd$n.tsd
        res$tsd.shift5 = tsd$shift5; res$tsd.shift3 = tsd$shift3
        beg = shiftCol(beg, tsd$shift5)
        end = shiftCol(end, tsd$shift3)
      }
    }

    # ---- 3. consensus ----
    s.new = nt2seq(cons.col[beg:end][cons.col[beg:end] != '-'])
    res$len.cons = nchar(s.new)

    # ---- 4. quality ----
    mx.el = rbind(cons = cons.col[beg:end], mx[, beg:end, drop = FALSE])
    id = svAlnIdent(mx.el, 'cons')
    res$ident.med = round(median(id$ident), 3)
    res$ident.q10 = round(unname(quantile(id$ident, 0.1)), 3)
    res$cov.med = round(median(id$cov), 3)

    # Sharpness: the copies agree inside and do not agree outside. The drop is what
    # matters and not the level outside: the aligner pulls unrelated flanks together,
    # and with a few copies they agree at the level of 0.5 just by that. The windows
    # are counted in the nucleotides of the consensus: a column of an insertion in
    # one copy agrees on nothing and says nothing about the boundary.
    w = opt$q.wnd
    meanAgree <- function(i) if(length(i) == 0) NA else round(mean(bn$agree[i]), 3)
    res$agree.in5  = meanAgree(head(col.nt[(col.nt >= beg) & (col.nt <= end)], w))
    res$agree.out5 = meanAgree(tail(col.nt[col.nt < beg], w))
    res$agree.in3  = meanAgree(tail(col.nt[(col.nt >= beg) & (col.nt <= end)], w))
    res$agree.out3 = meanAgree(head(col.nt[col.nt > end], w))

    # how much the consensus has changed
    mx.old = svAlnMafft(c(new = s.new, old = s.cons), path.work = tempdir(), args = opt$mafft.args)
    id.old = svAlnIdent(mx.old, 'new')
    res$old.ident = round(id.old$ident, 3)
    res$old.cov = round(id.old$cov, 3)

    flags = c(if(res$open5) 'open5', if(res$open3) 'open3',
              if(res$ident.med < opt$q.ident) 'low_ident',
              if(is.na(res$agree.out5) || (res$agree.in5 - res$agree.out5 < opt$q.contrast)) 'blurred5',
              if(is.na(res$agree.out3) || (res$agree.in3 - res$agree.out3 < opt$q.contrast)) 'blurred3')
    res$quality = if(length(flags) == 0) 'good' else paste(flags, collapse = ',')
    res$status = 'ok'

    cons = s.new
    saveRDS(list(copies = copies[idx, ], mx = mx, bounds = c(beg = beg, end = end),
                 agree = bn$agree, cov = bn$cov, ident = id, res = res),
            file.path(path.ext, paste0('fam_', f.id, '.rds')))

  }, error = function(e){
    res$status <<- paste0('error: ', conditionMessage(e))
  })

  list(res = res, cons = cons)
}

# ***********************************************************************
# ---- All the families ----

pokaz('Processing', length(fam.ids), 'families on', opt$cores, 'cores...')
t0 = Sys.time()
out = mclapply(fam.ids, oneFamily, mc.cores = opt$cores, mc.preschedule = FALSE)
pokaz('Done in', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), 'minutes')

tab = as.data.frame(rbindlist(lapply(out, function(x) x$res), fill = TRUE))

# a family that could not be refined keeps its consensus
cons = cons.old
for(x in out) if(!is.null(x$cons)) cons[[x$res$family]] = x$cons
tab$cons.src = ifelse(tab$status == 'ok', 'ext', 'old')

file.tab = file.path(path.out, paste0('sv_families_cons_ext_', suff, '.txt'))
write.table(tab, file.tab, sep = '\t', quote = FALSE, row.names = FALSE)
file.fasta = file.path(path.out, paste0('sv_families_cons_ext_', suff, '.fasta'))
writeFasta(setNames(cons, paste0('fam_', names(cons), '|', nchar(cons))), file.fasta)

pokaz('Table:', file.tab)
pokaz('Consensus sequences:', file.fasta)
pokaz('Refined:', sum(tab$cons.src == 'ext'), 'of', nrow(tab))
if('quality' %in% colnames(tab)){
  cnt = table(tab$quality[tab$cons.src == 'ext'])
  for(s in names(cnt)) pokaz('  ', s, ':', cnt[[s]])
}
n.err = sum(tab$status != 'ok')
if(n.err > 0){
  pokazAttention('Families not refined:', n.err)
  print(head(sort(table(sub('^error: ', '', tab$status[tab$status != 'ok'])), decreasing = TRUE), 5))
}
