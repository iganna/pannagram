# Alignment, consensus and structural order of every family of SVs.
#
# For every family: the members are brought to one strand, the (nearly)
# full-length ones are aligned by MAFFT, the edges of the alignment are trimmed
# by the coverage of the columns, and the consensus is checked for the terminal
# features of mobile elements.
#
# Every alignment is stored, so that the consensus and the classification can be
# recomputed with other thresholds without running MAFFT again.

suppressMessages({
  library(optparse)
  library(pannagram)
  library(parallel)
  library(data.table)
})

args = commandArgs(trailingOnly = TRUE)

option_list = list(
  make_option("--path.sv",     type = "character", default = NULL, help = "Path to the sv dir of the project"),
  make_option("--path.out",    type = "character", default = NULL, help = "Path for the results"),
  make_option("--similarity",  type = "integer",   default = 85,   help = "Similarity of the families"),
  make_option("--coverage",    type = "integer",   default = 85,   help = "Coverage of the families"),
  make_option("--cores",       type = "integer",   default = 1,    help = "Number of cores"),
  make_option("--min.members", type = "integer",   default = 3,    help = "Minimum number of members in a family"),
  make_option("--n.repr",      type = "integer",   default = 20,   help = "Maximum number of members to align"),
  make_option("--min.votes",   type = "integer",   default = 2,    help = "Members that have to agree for the class to be given"),
  make_option("--path.graphs", type = "character", default = NULL, help = "Directory with the per-family pieces, default <path.sv>/graphs_<sim>_<cov>/"),
  make_option("--trim.cov",    type = "integer",   default = 3,    help = "Minimum coverage of a column"),
  make_option("--trim.wnd",    type = "integer",   default = 5,    help = "Size of the clean window"),
  make_option("--trim.agree",  type = "double",    default = 0.7,  help = "Minimum share of the dominant nucleotide in a column"),
  make_option("--use.gap",     type = "logical",   default = TRUE, help = "Let a gap compete for the column"),
  make_option("--mafft.args",  type = "character", default = "--op 5 --ep 0.2 --quiet", help = "Arguments of MAFFT"),
  make_option("--path.code",   type = "character", default = "",   help = "Directory with the sources, if the package is not rebuilt yet")
)

opt = parse_args(OptionParser(option_list = option_list), args = args)

# The package may be not rebuilt yet: take the functions from the repository
if(nzchar(opt$path.code)){
  for(f in c('analys/graph_func.R', 'analys/sv_te_func.R', 'analys/sv_aln_func.R')){
    p = file.path(opt$path.code, f)
    if(file.exists(p)) source(p)
  }
}

path.sv = opt$path.sv
if(is.null(path.sv) || !dir.exists(path.sv)) stop(paste0('No SV dir: ', path.sv))
path.out = opt$path.out
if(is.null(path.out)) stop('path.out is not provided')
path.aln = file.path(path.out, 'aln')
dir.create(path.aln, showWarnings = FALSE, recursive = TRUE)
path.vote = file.path(path.out, 'votes')
dir.create(path.vote, showWarnings = FALSE, recursive = TRUE)

suff = paste0(opt$similarity, '_', opt$coverage)
file.seqs = file.path(path.sv, 'seq_sv_large.fasta')
file.fam  = file.path(path.sv, paste0('sv_families_', suff, '.txt'))
for(f in c(file.seqs, file.fam)) if(!file.exists(f)) stop(paste0('No file: ', f))

# ***********************************************************************
# ---- Reading ----

pokaz('Reading the sequences...')
seqs = readFasta(file.seqs)
pokaz('Sequences:', length(seqs))

fam = fread(file.fam, sep = '\t', header = FALSE, data.table = FALSE)
fam[[2]] = as.character(fam[[2]])
fam = fam[fam[[1]] %in% names(seqs), ]
n.mem = table(fam[[2]])
fam.ids = names(n.mem)[n.mem >= opt$min.members]
pokaz('Families:', length(fam.ids), 'of', length(n.mem), 'with at least', opt$min.members, 'members')

path.graphs = if(is.null(opt$path.graphs)) file.path(path.sv, paste0('graphs_', suff)) else opt$path.graphs
if(!dir.exists(path.graphs)) stop(paste0('No per-family pieces: ', path.graphs,
                                         '. Run sv_03b_family_graphs.R first.'))
pokaz('Per-family pieces:', path.graphs)

sv.by.fam = split(fam[[1]], fam[[2]])

# ***********************************************************************
# ---- One family ----

oneFamily <- function(f.id){

  res = data.frame(family = f.id, n.members = length(sv.by.fam[[f.id]]),
                   n.aln = 0, len.cons = 0, aln.len = 0, trim.beg = NA, trim.end = NA,
                   cons.src = '', ori.by = '', te.class = 'failed', status = '', stringsAsFactors = FALSE)
  cons = NULL

  tryCatch({
    sv = sv.by.fam[[f.id]]
    file.aln = file.path(path.aln, paste0('fam_', f.id, '.rds'))

    # ---- 1. the alignment, stored ----
    if(file.exists(file.aln)){
      st = readRDS(file.aln)
      res$ori.by = st$ori.by
    } else {
      file.g = file.path(path.graphs, paste0('fam_', f.id, '.txt'))
      mutual = NULL
      if(file.exists(file.g)){
        g = read.table(file.g, sep = '\t', header = TRUE, stringsAsFactors = FALSE,
                       quote = '', comment.char = '')
        # a mutual pair is an edge present in both directions in the final graph
        i.both = g$in.graph & g$in.graph.rev
        sv.both = unique(c(g$name.query[i.both], g$name.target[i.both]))
        mutual = setNames(sv %in% sv.both, sv)
      }

      seqs.ori = NULL
      ori.by = 'graph'
      if(file.exists(file.g)){
        seqs.ori = tryCatch(getComponentSequences(sv, seqs, path.graphs = path.graphs, family = f.id),
                            error = function(e) NULL)
      }
      if(is.null(seqs.ori)){
        seqs.ori = svFamilyOrient(seqs[sv])$seqs
        ori.by = 'kmer'
      }

      idx = svFamilyRepr(nchar(seqs.ori), n.repr = opt$n.repr,
                         mutual = if(is.null(mutual)) NULL else mutual[names(seqs.ori)])
      if(length(idx) < 2) stop('less than two members to align')

      # the rows go from the longest one down, so that the picture of the alignment
      # can be read: the length of a member is where its own story begins
      idx = idx[order(nchar(seqs.ori[idx]), decreasing = TRUE)]

      mx = svAlnMafft(seqs.ori[idx], path.work = tempdir(), args = opt$mafft.args)
      st = list(mx = mx, ori.by = ori.by, n.mutual = if(is.null(mutual)) NA else sum(mutual),
                strand = ifelse(toupper(seqs.ori[idx]) == toupper(seqs[names(seqs.ori)[idx]]), '+', '-'))
      saveRDS(st, file.aln)
      res$ori.by = ori.by
    }

    res$n.aln = nrow(st$mx)
    res$aln.len = ncol(st$mx)

    # ---- 2. the consensus, for the sequence itself ----
    tr = svAlnTrim(st$mx, min.cov = opt$trim.cov, wnd = opt$trim.wnd,
                   min.agree = opt$trim.agree, use.gap = opt$use.gap)
    s.cons = if(is.null(tr$mx)) svAlnCons(st$mx)$seq else svAlnCons(tr$mx)$seq
    if(nchar(s.cons) < 40) s.cons = svAlnCons(st$mx)$seq
    if(nchar(s.cons) < 40) stop('the consensus is too short')
    res$len.cons = nchar(s.cons)
    if(!is.null(tr$mx)){ res$trim.beg = tr$beg; res$trim.end = tr$end }
    cons = setNames(s.cons, paste0('fam_', f.id, '|', res$len.cons))

    # ---- 3. the structural order: every member of the family votes ----
    # A consensus can break where the copies are clean -- two length outliers
    # stretch the alignment and the terminal repeat that every copy carries is
    # lost -- and it can also invent a repeat out of two flanks that only one
    # member each covers. So the verdict belongs to the members themselves, and
    # all of them vote, not only the ones picked for the alignment. The detectors
    # look for repeats inside a sequence, so the strand does not matter here.
    lv = c('LTR', 'TIR', 'Helitron', 'LINE_like')
    file.v = file.path(path.vote, paste0('fam_', f.id, '.rds'))
    # The whole verdict of every member is kept, not just its class, so that a
    # change in how the evidence is weighed does not cost another full pass.
    if(file.exists(file.v)){
      v.m = readRDS(file.v)
    } else {
      v.m = rbindlist(lapply(names(seqs[sv]), function(nm){
        x = toupper(seqs[[nm]])
        if(nchar(x) < 40) data.frame(sv = nm, len = nchar(x), te.class = 'unknown')
        else svTeOrder(x, nm)
      }), fill = TRUE)
      saveRDS(v.m, file.v)
    }
    cl.m = v.m$te.class
    n.cl = table(factor(cl.m, levels = c(lv, 'unknown')))
    res$n.vote = length(cl.m)

    # The consensus is read in both variants -- the trimming cuts the flanks off
    # one family and the very termini off another -- and the stronger one is kept.
    rank.cls = c('LTR', 'TIR', 'Helitron', 'LINE_like', 'unknown')
    cand = list(full = svAlnCons(st$mx)$seq)
    if(!is.null(tr$mx)) cand$trim = s.cons
    o.cons = NULL
    for(nm in names(cand)){
      if(nchar(cand[[nm]]) < 40) next
      o = svTeOrder(cand[[nm]], name = paste0('fam_', f.id))
      if(is.null(o.cons) || match(o$te.class, rank.cls) < match(o.cons$te.class, rank.cls)){
        o.cons = o; res$cons.src = nm
      }
    }
    if(is.null(o.cons)) o.cons = svTeOrder(s.cons, name = paste0('fam_', f.id))
    res = cbind(res[, setdiff(colnames(res), 'te.class')], o.cons[, -(1:2)])
    colnames(res)[colnames(res) == 'te.class'] = 'cons.class'

    res$n.LTR = n.cl[['LTR']]; res$n.TIR = n.cl[['TIR']]
    res$n.Helitron = n.cl[['Helitron']]; res$n.LINE = n.cl[['LINE_like']]
    res$n.unknown = n.cl[['unknown']]

    # A family where nobody carries any structure has no ballots to count at all,
    # and that is a different answer from a family whose members disagree. Two
    # agreeing members are enough to name the class, but they have to be the
    # dominant one: a tie with another class leaves the family unnamed.
    n.cls = sum(n.cl[lv])
    v = sort(as.integer(n.cl[lv]), decreasing = TRUE)
    i.top = which.max(n.cl[lv])
    res$n.cls = n.cls
    res$te.class = if(n.cls == 0) 'cant_vote' else
                   if(v[1] >= opt$min.votes && v[1] > v[2]) lv[i.top] else 'unknown'
    res$support = if(res$te.class %in% lv) round(v[1] / n.cls, 3) else 0
    res$status = 'ok'

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

# the rows of the failed families have fewer columns, so they are filled in
tab = as.data.frame(rbindlist(lapply(out, function(x) x$res), fill = TRUE))
cons = unlist(lapply(out, function(x) x$cons))

file.tab = file.path(path.out, paste0('sv_te_order_cons_', suff, '.txt'))
write.table(tab, file.tab, sep = '\t', quote = FALSE, row.names = FALSE)
if(length(cons) > 0) writeFasta(cons, file.path(path.out, paste0('sv_families_cons_', suff, '.fasta')))

pokaz('Table:', file.tab)
pokaz('Consensus sequences:', length(cons))
if('te.class' %in% colnames(tab)){
  cnt = table(tab$te.class)
  for(s in names(cnt)) pokaz('  ', s, ':', cnt[[s]])
}
n.err = sum(tab$status != 'ok')
if(n.err > 0){
  pokazAttention('Families with errors:', n.err)
  print(head(sort(table(sub(':.*', '', tab$status[tab$status != 'ok'])), decreasing = TRUE), 5))
}
