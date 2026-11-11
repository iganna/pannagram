# Subfamilies of every family of SVs and their refined consensus sequences.
#
# Takes the alignments stored by sv_08_family_cons.R. For every family: the
# aligned representatives are split into subfamilies, a consensus is called for
# every subfamily, and then all the members of the family are aligned to these
# consensus sequences, assigned to the closest one, and the consensus is called
# again from its members until it stops changing.

suppressMessages({
  library(optparse)
  library(pannagram)
  library(parallel)
  library(data.table)
})

args = commandArgs(trailingOnly = TRUE)

option_list = list(
  make_option("--path.sv",     type = "character", default = NULL, help = "Path to the sv dir of the project"),
  make_option("--path.out",    type = "character", default = NULL, help = "Path for the results, the same as of sv_08_family_cons.R"),
  make_option("--similarity",  type = "integer",   default = 85,   help = "Similarity of the families"),
  make_option("--coverage",    type = "integer",   default = 85,   help = "Coverage of the families"),
  make_option("--cores",       type = "integer",   default = 1,    help = "Number of cores"),
  make_option("--min.members", type = "integer",   default = 3,    help = "Minimum number of members in a subfamily"),
  make_option("--split.dist",  type = "double",    default = 0.1,  help = "Distance at which a family is split into subfamilies"),
  make_option("--min.overlap", type = "integer",   default = 50,   help = "Minimum overlap of two members to be compared"),
  make_option("--min.ident",   type = "double",    default = 0.7,  help = "Minimum identity of a member to the consensus of its subfamily"),
  make_option("--n.iter",      type = "integer",   default = 3,    help = "Maximum number of rounds of the refinement"),
  make_option("--n.max",       type = "integer",   default = 100,  help = "Maximum number of members voting for the consensus"),
  make_option("--min.len",     type = "integer",   default = 40,   help = "Minimum length of the consensus"),
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
  for(f in c('analys/sv_aln_func.R')){
    p = file.path(opt$path.code, f)
    if(file.exists(p)) source(p)
  }
}

path.sv = opt$path.sv
if(is.null(path.sv) || !dir.exists(path.sv)) stop(paste0('No SV dir: ', path.sv))
path.out = opt$path.out
if(is.null(path.out)) stop('path.out is not provided')
path.aln = file.path(path.out, 'aln')
if(!dir.exists(path.aln)) stop(paste0('No alignments: ', path.aln, '. Run sv_08_family_cons.R first.'))
path.refine = file.path(path.out, 'refine')
dir.create(path.refine, showWarnings = FALSE, recursive = TRUE)

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
sv.by.fam = split(fam[[1]], fam[[2]])

fam.ids = sub('^fam_', '', sub('\\.rds$', '', list.files(path.aln, pattern = '^fam_.*\\.rds$')))
fam.ids = fam.ids[fam.ids %in% names(sv.by.fam)]
pokaz('Families with an alignment:', length(fam.ids))

# ***********************************************************************
# ---- One family ----

oneFamily <- function(f.id){

  res = data.frame(family = f.id, subfam = NA, n.members = length(sv.by.fam[[f.id]]),
                   n.repr = 0, n.assigned = 0, len.cons = 0, ident.med = NA, cov.med = NA,
                   n.iter = 0, converged = FALSE, status = '', stringsAsFactors = FALSE)
  cons = NULL
  mem = NULL

  tryCatch({
    file.res = file.path(path.refine, paste0('fam_', f.id, '.rds'))
    if(file.exists(file.res)){
      rf = readRDS(file.res)
    } else {
      st = readRDS(file.path(path.aln, paste0('fam_', f.id, '.rds')))
      mx = st$mx

      # ---- 1. subfamilies of the representatives ----
      cl = svFamilySplit(mx, max.dist = opt$split.dist, min.size = 2, min.overlap = opt$min.overlap)
      cons.0 = character(0)
      for(k in sort(unique(cl[!is.na(cl)]))){
        cons.0[as.character(k)] = svAlnConsTrim(mx[which(cl == k), , drop = FALSE], min.len = opt$min.len,
                                                trim.cov = opt$trim.cov, trim.wnd = opt$trim.wnd,
                                                trim.agree = opt$trim.agree, use.gap = opt$use.gap)
      }
      if(length(cons.0) == 0) stop('no subfamily')

      # ---- 2. all the members on the strand of the alignment ----
      # The consensus of the whole alignment is the reference: it is on the
      # strand of the representatives, and every subfamily is on the same one.
      s.ref = svAlnCons(mx)$seq
      seqs.ori = svFamilyOrient(seqs[sv.by.fam[[f.id]]], ref = s.ref)

      # ---- 3. refinement by all the members ----
      rf = svFamilyRefine(seqs.ori$seqs, cons.0, n.iter = opt$n.iter, n.max = opt$n.max,
                          min.members = opt$min.members, min.ident = opt$min.ident, min.len = opt$min.len,
                          trim.cov = opt$trim.cov, trim.wnd = opt$trim.wnd,
                          trim.agree = opt$trim.agree, use.gap = opt$use.gap,
                          path.work = tempdir(), mafft.args = opt$mafft.args)
      rf$best$strand = seqs.ori$strand[rf$best$sv]
      rf$n.repr = table(factor(cl, levels = names(cons.0)))
      saveRDS(rf, file.res)
    }

    # ---- 4. the summary by subfamilies ----
    best = rf$best
    mem = data.frame(family = f.id, best, stringsAsFactors = FALSE)

    res = do.call(rbind, lapply(names(rf$cons), function(k){
      b = best[which(best$subfam == k), ]
      data.frame(family = f.id, subfam = k, n.members = nrow(best),
                 n.repr = if(k %in% names(rf$n.repr)) as.integer(rf$n.repr[[k]]) else 0,
                 n.assigned = nrow(b), len.cons = nchar(rf$cons[[k]]),
                 ident.med = round(median(b$ident), 3), cov.med = round(median(b$cov), 3),
                 n.iter = rf$n.iter, converged = rf$converged, status = 'ok',
                 stringsAsFactors = FALSE)
    }))
    cons = setNames(rf$cons, paste0('fam_', f.id, '.', names(rf$cons), '|', nchar(rf$cons)))

  }, error = function(e){
    res$status <<- paste0('error: ', conditionMessage(e))
  })

  list(res = res, cons = cons, mem = mem)
}

# ***********************************************************************
# ---- All the families ----

pokaz('Processing', length(fam.ids), 'families on', opt$cores, 'cores...')
t0 = Sys.time()
out = mclapply(fam.ids, oneFamily, mc.cores = opt$cores, mc.preschedule = FALSE)
pokaz('Done in', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), 'minutes')

tab = as.data.frame(rbindlist(lapply(out, function(x) x$res), fill = TRUE))
mem = as.data.frame(rbindlist(lapply(out, function(x) x$mem), fill = TRUE))
cons = unlist(lapply(out, function(x) x$cons))

file.tab = file.path(path.out, paste0('sv_subfam_', suff, '.txt'))
file.mem = file.path(path.out, paste0('sv_subfam_members_', suff, '.txt'))
write.table(tab, file.tab, sep = '\t', quote = FALSE, row.names = FALSE)
write.table(mem, file.mem, sep = '\t', quote = FALSE, row.names = FALSE)
if(length(cons) > 0) writeFasta(cons, file.path(path.out, paste0('sv_subfam_cons_', suff, '.fasta')))

pokaz('Table:', file.tab)
pokaz('Members:', file.mem)
pokaz('Subfamilies:', sum(tab$status == 'ok'), 'in', length(unique(tab$family[tab$status == 'ok'])), 'families')
if(nrow(mem) > 0) pokaz('Unassigned members:', sum(is.na(mem$subfam)), 'of', nrow(mem))
n.err = sum(tab$status != 'ok')
if(n.err > 0){
  pokazAttention('Families with errors:', n.err)
  print(head(sort(table(sub(':.*', '', tab$status[tab$status != 'ok'])), decreasing = TRUE), 5))
}
