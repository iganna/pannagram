# Annotation groups
#
# INTERVAL FORMAT (_n).
# Input:  GFF files of the accessions (<acc>.gff; not every accession needs one)
#         and the pangenome alignment.
# Output: annogroups.txt     - groups of genes that occupy the same place of the pangenome
#         correspondence.txt - group x accession -> IDs of the genes
#         lost_genes.txt     - genes that were not placed, with the reason
#         projection.txt     - models projected to the accessions that lack a gene of a group
#         gff/<acc>.gff      - annotation of every accession: own genes + projected ones, common IDs
#         pangen.gff         - all the annotations in the pangenome coordinates

suppressMessages({
  library(rhdf5)
  library(optparse)
  library(pannagram)
  library(crayon)
})

source(system.file("utils/chunk_hdf5.R", package = "pannagram"))
source(system.file("utils/interval_func.R", package = "pannagram")) # interval codec + hdf5 layer
source(system.file("analys/annogroup_func_n.R", package = "pannagram"))

# ***********************************************************************
args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.features.msa", type = "character", default = NULL, help = "Path to msa dir (features)"),
  make_option("--path.annot",   type = "character", default = NULL,   help = "Path to the GFF files, <acc>.gff"),
  make_option("--path.res",     type = "character", default = NULL,   help = "Path to the results"),
  make_option("--path.chr",     type = "character", default = NULL,   help = "Path to the chromosomes (to check the projected models)"),
  make_option("--ref",          type = "character", default = NULL,   help = "prefix of the reference file"),
  make_option("--aln.type",     type = "character", default = aln.type.msa, help = "type of alignment ('pan', 'ref', etc.)"),
  make_option("--s.chr",        type = "character", default = '_Chr', help = "Chromosome name format in GFF: '<anything><s.chr><number>'"),
  make_option("--feature",      type = "character", default = 'all',  help = "Types of the top-level features to group, comma-separated ('all' by default)"),
  make_option("--id.prefix",    type = "character", default = 'AG',   help = "Prefix of the group IDs"),
  make_option("--overlap",      type = "numeric",   default = 0.5,    help = "Genes are in one group if they overlap by at least this share of the longer one"),
  make_option("--snap",         action = "store_true", default = FALSE, help = "Move an unaligned end of a gene inwards instead of losing the gene"),
  make_option("--path.log",     type = "character", default = NULL,   help = "Path for log files"),
  make_option("--log.level",    type = "character", default = NULL,   help = "Level of log to be shown on the screen")
)

opt = parse_args(OptionParser(option_list = option_list), args = args)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ***********************************************************************
# ---- Paths and parameters ----

isNull <- function(x) is.null(x) || (x == 'NULL') || (x == '')
addSlash <- function(x) paste0(sub('/+$', '', x), '/')

if(isNull(opt$path.features.msa)) stop("'path.features.msa' must be specified.")
if(isNull(opt$path.annot)) stop("'path.annot' must be specified.")
if(isNull(opt$path.res)) stop("'path.res' must be specified.")

path.features.msa <- addSlash(opt$path.features.msa)
path.annot <- addSlash(opt$path.annot)
path.res <- addSlash(opt$path.res)
if(!dir.exists(path.features.msa)) stop(paste('No alignment directory found:', path.features.msa))
if(!dir.exists(path.annot)) stop(paste('No annotation directory found:', path.annot))

path.gff <- paste0(path.res, 'gff/')
if(!dir.exists(path.gff)) dir.create(path.gff, recursive = TRUE)

path.chr <- if(isNull(opt$path.chr)) NULL else addSlash(opt$path.chr)
if(!is.null(path.chr) && !dir.exists(path.chr)){
  pokazAttention('No directory with chromosomes:', path.chr)
  path.chr <- NULL
}
if(is.null(path.chr)) pokazAttention('Chromosomes are not provided, coding sequences of the projected models are not checked.')

s.chr <- opt$s.chr
id.prefix <- opt$id.prefix
min.overlap <- opt$overlap
snap <- isTRUE(opt$snap)
features <- strsplit(opt$feature, ',', fixed = TRUE)[[1]]

# ***********************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

aln.type <- opt$aln.type

ref.name <- opt[['ref']]  # not `$`: it partially matches longer option names
if(is.null(ref.name) || ref.name == "NULL") ref.name <- ''

# Common code for aln.pref, ref.suffix and s.combinations
source(system.file("utils/chunk_combinations.R", package = "pannagram"))

s.combinations <- s.combinations[order(sapply(s.combinations, comb2ref))]

# ***********************************************************************
# ---- Accessions and annotations ----

accessions <- c()
for(s.comb in s.combinations){
  accessions <- union(accessions, h5AccNames(paste0(path.features.msa, aln.pref, s.comb, ref.suff, '.h5')))
}

files.gff <- list.files(path = path.annot, pattern = "\\.(gff|gff3)$", full.names = FALSE)
acc.files <- sub("\\.(gff|gff3)$", "", files.gff)
if(any(!(acc.files %in% accessions))){
  pokazAttention('Annotations of accessions that are not in the alignment are skipped:',
                 acc.files[!(acc.files %in% accessions)])
}
files.gff <- files.gff[acc.files %in% accessions]
acc.files <- acc.files[acc.files %in% accessions]
if(length(files.gff) == 0){
  pokazAttention('No GFF files for the aligned accessions in', path.annot, '(expected <accession>.gff). Nothing to do.')
  quit(save = "no", status = 0)
}

pokaz('Accessions in the alignment:', length(accessions), '; with annotation:', length(acc.files))

gff.list <- list()
for(i.file in seq_along(files.gff)){
  acc <- acc.files[i.file]
  gff <- annoReadGFF(paste0(path.annot, files.gff[i.file]))
  if(nrow(gff) == 0){
    pokazAttention('No features in the annotation of', acc)
    next
  }
  gff <- annoHierarchy(gff)
  gff$chr <- annoChr(gff$V1, s.chr)
  if(all(is.na(gff$chr))){
    pokazAttention('Sequence names of', acc, 'do not match the format', paste0('"*', s.chr, 'N"'),
                   ', for example:', gff$V1[1], '- the annotation is skipped. Check the chromosome format.')
    next
  }

  gff$is.root <- gff$root == seq_len(nrow(gff))
  is.used <- gff$is.root & (('all' %in% features) | (gff$V3 %in% features))
  gff$use <- is.used[gff$root]

  # Place in the pangenome, filled in chromosome by chromosome
  gff$pan.chr <- 0L
  gff$pan.beg <- 0L
  gff$pan.end <- 0L
  gff$delta <- 0L
  gff$snapped <- FALSE
  gff$reason <- ifelse(is.na(gff$chr), 'chr_format', 'chr_not_aligned')

  # Groups of the root: new ID, the primary group and its order, all the groups
  gff$root.id <- NA_character_
  gff$group <- NA_character_
  gff$group.ord <- NA_integer_
  gff$groups <- NA_character_

  gff.list[[acc]] <- gff
  pokaz('Accession', acc, ': features', nrow(gff), ', top-level', sum(is.used))
}
acc.gff <- names(gff.list)
if(length(acc.gff) == 0){
  pokazAttention('No annotations to work with. Nothing to do.')
  quit(save = "no", status = 0)
}

# Names of the chromosomes of an accession, as they are in its annotation
chrName <- function(acc, i.chr){
  gff <- gff.list[[acc]]
  if(!is.null(gff)){
    s <- gff$V1[which(gff$chr == i.chr)]
    if(length(s) > 0) return(names(which.max(table(s))))
  }
  paste0(acc, '_Chr', i.chr)
}

# ***********************************************************************
# ---- Main loop: chromosomes ----

groups.all <- list()
corresp.all <- list()
proj.tbl.all <- list()
proj.gff.all <- list()
n.groups <- 0

for(s.comb in s.combinations){

  q.chr <- as.integer(strsplit(s.comb, '_')[[1]][1])
  p.chr <- comb2ref(s.comb)
  pokaz('Chromosome', s.comb)

  file.comb <- paste0(path.features.msa, aln.pref, s.comb, ref.suff, '.h5')
  accs.comb <- h5AccNames(file.comb)
  ivl.list <- list()
  for(acc in accs.comb) ivl.list[[acc]] <- annoIvlPrep(h5IvlRead(file.comb, acc))

  # ---- Genes -> pangenome ----
  roots <- list()
  rows.chr <- list()
  for(acc in acc.gff){
    gff <- gff.list[[acc]]
    idx <- which((gff$chr == q.chr) & gff$use)
    if(length(idx) == 0) next
    if(!(acc %in% accs.comb)) next

    # Top-level features define the groups: exact ends, unless --snap.
    # The rest is for drawing and projection: ends are moved inwards.
    x <- annoLift(ivl.list[[acc]], gff$V4[idx], gff$V5[idx], snap = TRUE)
    if(!snap){
      is.r <- gff$is.root[idx]
      x[is.r, ] <- annoLift(ivl.list[[acc]], gff$V4[idx[is.r]], gff$V5[idx[is.r]], snap = FALSE)
    }
    gff$pan.chr[idx] <- p.chr
    gff$pan.beg[idx] <- x$pan.beg
    gff$pan.end[idx] <- x$pan.end
    gff$delta[idx] <- x$delta
    gff$snapped[idx] <- x$snapped
    gff$reason[idx] <- x$reason
    gff.list[[acc]] <- gff

    idx <- idx[gff$reason[gff$root[idx]] == '']  # children of the lost genes are out
    rows.chr[[acc]] <- idx
    idx.r <- idx[gff$is.root[idx]]
    if(length(idx.r) == 0) next
    roots[[acc]] <- data.frame(acc = acc, row = idx.r,
                               beg = gff$pan.beg[idx.r], end = gff$pan.end[idx.r],
                               strand = annoStrand(gff$V7[idx.r], gff$delta[idx.r]),
                               type = gff$V3[idx.r],
                               id = ifelse(is.na(gff$id[idx.r]),
                                           paste0(gff$V3[idx.r], ':', gff$V1[idx.r], ':', gff$V4[idx.r], '-', gff$V5[idx.r]),
                                           gff$id[idx.r]),
                               stringsAsFactors = FALSE)
  }
  if(length(roots) == 0){
    pokazAttention('No genes were placed on the chromosome', s.comb)
    next
  }
  roots <- do.call(rbind, roots)
  rownames(roots) <- NULL

  # ---- Groups ----
  groups <- list()
  map <- list()
  n.prev <- 0
  for(s.s in unique(roots$strand)){
    idx.s <- which(roots$strand == s.s)
    x <- annoGroups(roots$beg[idx.s], roots$end[idx.s], cutoff = min.overlap)
    a <- annoAssign(roots$beg[idx.s], roots$end[idx.s], x$group, x$groups)
    a$gene <- idx.s[a$gene]
    a$group <- a$group + n.prev
    n.prev <- n.prev + nrow(x$groups)
    x$groups$strand <- s.s
    groups[[s.s]] <- x$groups
    map[[s.s]] <- a
  }
  groups <- do.call(rbind, groups)
  map <- do.call(rbind, map)

  # Groups are numbered along the chromosome
  ord <- order(groups$beg, groups$end)
  map$group <- match(map$group, ord)
  groups <- groups[ord, , drop = FALSE]
  rownames(groups) <- NULL
  n.gr <- nrow(groups)
  groups$group <- sprintf('%s%dG%06d', id.prefix, p.chr, seq_len(n.gr))
  groups$chr <- p.chr
  groups$ord <- n.groups + seq_len(n.gr)
  n.groups <- n.groups + n.gr

  # A gene is a primary member of its own group and a secondary one of the
  # groups it covers; the primary members go first
  map <- map[order(map$group, !map$primary, roots$beg[map$gene]), , drop = FALSE]
  map$acc <- roots$acc[map$gene]
  map.pr <- map[map$primary, , drop = FALSE]
  s.key <- paste(map.pr$group, map.pr$acc)
  map.pr$multi <- s.key %in% s.key[duplicated(s.key)]

  n.cl <- table(groups$cluster, groups$strand)
  groups$n.split <- as.integer(n.cl[cbind(as.character(groups$cluster), groups$strand)])
  f.gr <- factor(map.pr$group, levels = seq_len(n.gr))
  groups$n.acc <- as.integer(tapply(map.pr$acc, f.gr, function(x) length(unique(x))))
  groups$n.genes <- as.integer(table(f.gr))
  groups$n.acc.multi <- as.integer(tapply(map.pr$acc[map.pr$multi], f.gr[map.pr$multi], function(x) length(unique(x))))
  groups$n.acc.multi[is.na(groups$n.acc.multi)] <- 0L
  groups$n.shared <- as.integer(table(factor(map$group[!map$primary], levels = seq_len(n.gr))))
  groups$types <- as.character(tapply(roots$type[map.pr$gene], f.gr, function(x) paste(sort(unique(x)), collapse = ',')))

  # Secondary members are in brackets
  corresp <- tapply(ifelse(map$primary, roots$id[map$gene], paste0('(', roots$id[map$gene], ')')),
                    list(factor(map$group, levels = seq_len(n.gr)), factor(map$acc, levels = acc.gff)),
                    paste, collapse = ',')
  corresp[is.na(corresp)] <- '-'
  corresp.all[[s.comb]] <- data.frame(group = groups$group, corresp, check.names = FALSE, stringsAsFactors = FALSE)

  # ---- Common IDs of the genes ----
  o <- order(map$gene, !map$primary)
  gene.groups <- tapply(groups$group[map$group[o]], map$gene[o], paste, collapse = ',')
  map.pr$id.new <- paste0(groups$group[map.pr$group], '.', map.pr$acc)
  map.pr$id.new[map.pr$multi] <- paste0(map.pr$id.new[map.pr$multi], '.', annoRunNum(s.key)[map.pr$multi])
  for(acc in unique(map.pr$acc)){
    k <- which(map.pr$acc == acc)
    r <- roots$row[map.pr$gene[k]]
    gff.list[[acc]]$root.id[r] <- map.pr$id.new[k]
    gff.list[[acc]]$group[r] <- groups$group[map.pr$group[k]]
    gff.list[[acc]]$group.ord[r] <- groups$ord[map.pr$group[k]]
    gff.list[[acc]]$groups[r] <- gene.groups[as.character(map.pr$gene[k])]
  }

  # ---- Models to project ----
  # A model of a group: the one with the most frequent ends, taken from an
  # accession that has exactly one gene in the group.
  cand <- map.pr
  cand$len <- roots$end[cand$gene] - roots$beg[cand$gene]
  s.key <- paste(cand$group, roots$beg[cand$gene], roots$end[cand$gene])
  cand$freq <- as.integer(table(s.key)[s.key])
  cand <- cand[order(cand$group,
                     -(!cand$multi),
                     -cand$freq, -cand$len,
                     match(cand$acc, accessions)), , drop = FALSE]
  models <- cand[!duplicated(cand$group), , drop = FALSE]
  models$row <- roots$row[models$gene]

  # ---- Projection: coordinates ----
  proj.list <- list()
  for(acc.t in accs.comb){
    need <- models[!(models$group %in% map$group[map$acc == acc.t]), , drop = FALSE]
    if(nrow(need) == 0) next

    proj <- list()
    for(acc.s in unique(need$acc)){
      gff <- gff.list[[acc.s]]
      r <- need$row[need$acc == acc.s]
      idx <- rows.chr[[acc.s]]
      idx <- idx[gff$root[idx] %in% r]
      proj[[acc.s]] <- data.frame(src.acc = acc.s, src.row = idx, src.root = gff$root[idx],
                                  src.parent = gff$parent.idx[idx],
                                  group = need$group[need$acc == acc.s][match(gff$root[idx], r)],
                                  gff[idx, c('V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'id',
                                             'pan.beg', 'pan.end', 'delta', 'snapped', 'reason')],
                                  stringsAsFactors = FALSE)
    }
    proj <- do.call(rbind, proj)
    rownames(proj) <- NULL

    is.lifted <- proj$reason == ''
    x <- annoProject(ivl.list[[acc.t]], proj$pan.beg[is.lifted], proj$pan.end[is.lifted])
    proj$beg <- 0L; proj$end <- 0L; proj$n.cov <- 0; proj$state <- 'absent'; proj$strand <- proj$V7
    proj$beg[is.lifted] <- x$beg
    proj$end[is.lifted] <- x$end
    proj$n.cov[is.lifted] <- x$n.cov
    proj$state[is.lifted] <- x$state
    proj$strand[is.lifted] <- annoStrand(annoStrand(proj$V7[is.lifted], proj$delta[is.lifted]), x$delta)

    # Transcript of a CDS piece
    proj$is.root <- proj$src.row == proj$src.root
    proj$tr <- paste(proj$src.acc, ifelse(is.na(proj$src.parent), proj$src.row, proj$src.parent), sep = '|')
    proj$complete <- (proj$state == 'ok') & !proj$snapped
    proj.list[[acc.t]] <- proj
  }

  # ---- Projection: coding sequences ----
  cds.src <- list()
  cds.trg <- list()
  if(!is.null(path.chr)){
    for(acc in accs.comb){
      # Source models of this accession and models projected to it
      src <- do.call(rbind, lapply(proj.list, function(p) p[(p$src.acc == acc) & (p$V3 == 'CDS'),
                                                           c('tr', 'V4', 'V5', 'V7', 'V8'), drop = FALSE]))
      if(!is.null(src)) src <- src[!duplicated(paste(src$tr, src$V4)), , drop = FALSE]
      trg <- proj.list[[acc]]
      if(!is.null(trg)){
        tr.ok <- tapply(trg$complete[trg$V3 == 'CDS'], trg$tr[trg$V3 == 'CDS'], all)
        trg <- trg[(trg$V3 == 'CDS') & (trg$tr %in% names(tr.ok)[tr.ok]), , drop = FALSE]
      }
      if((is.null(src) || (nrow(src) == 0)) && (is.null(trg) || (nrow(trg) == 0))) next

      file.chr <- paste0(path.chr, acc, '_chr', q.chr, '.fasta')
      if(!file.exists(file.chr)){
        pokazAttention('No chromosome file', file.chr)
        next
      }
      s.seq <- readFasta(file.chr)[1]

      if(!is.null(src) && (nrow(src) > 0)){
        cds.src[[acc]] <- annoCdsProps(s.seq, src$V4, src$V5, src$V7, src$V8, src$tr)
      }
      if(!is.null(trg) && (nrow(trg) > 0)){
        cds.trg[[acc]] <- annoCdsProps(s.seq, trg$beg, trg$end, trg$strand, trg$V8, trg$tr)
      }
      rm(s.seq)
    }
  }
  cds.src <- do.call(rbind, cds.src)

  # ---- Projection: status and GFF ----
  for(acc.t in names(proj.list)){
    proj <- proj.list[[acc.t]]
    is.cds <- proj$V3 == 'CDS'

    # Transcripts
    tr <- data.frame(tr = unique(proj$tr[is.cds]), stringsAsFactors = FALSE)
    tr$group <- proj$group[is.cds][match(tr$tr, proj$tr[is.cds])]
    tr$status <- rep('truncated', nrow(tr))
    tr.ok <- tapply(proj$complete[is.cds], proj$tr[is.cds], all)
    is.ok <- tr$tr %in% names(tr.ok)[tr.ok]
    tr$status[is.ok] <- 'unchecked'
    p.trg <- cds.trg[[acc.t]]
    if(!is.null(p.trg) && !is.null(cds.src)){
      k <- which(is.ok & (tr$tr %in% p.trg$tr) & (tr$tr %in% cds.src$tr))
      tr$status[k] <- annoCdsStatus(cds.src[match(tr$tr[k], cds.src$tr), , drop = FALSE],
                                    p.trg[match(tr$tr[k], p.trg$tr), , drop = FALSE])
    }

    # Models: the best transcript; without CDS - whether everything is in place
    md <- proj[proj$is.root, c('group', 'src.acc', 'id', 'V3', 'V4', 'V5', 'V7',
                               'beg', 'end', 'strand', 'n.cov', 'state'), drop = FALSE]
    md$status <- ifelse(as.logical(tapply(proj$complete, proj$group, all)[as.character(md$group)]),
                        'intact', 'truncated')
    if(nrow(tr) > 0){
      st <- tapply(match(tr$status, anno.status.rank), tr$group, min)
      k <- as.character(md$group) %in% names(st)
      md$status[k] <- anno.status.rank[st[as.character(md$group[k])]]
    }
    md$status[md$state %in% c('absent', 'rearranged')] <- md$state[md$state %in% c('absent', 'rearranged')]
    md$cov <- round(pmin(md$n.cov / (md$V5 - md$V4 + 1), 1), 3)

    s.name <- chrName(acc.t, q.chr)
    proj.tbl.all[[paste(s.comb, acc.t)]] <-
      data.frame(group = groups$group[md$group], acc = acc.t,
                 src.acc = md$src.acc,
                 src.id = ifelse(is.na(md$id), '-', md$id),
                 chr = s.name, beg = md$beg, end = md$end,
                 strand = ifelse(md$beg > 0, md$strand, '.'),
                 status = md$status, cov = md$cov, stringsAsFactors = FALSE)

    # Features to write: the gene is in place, the feature and its parent too
    keep <- (proj$state %in% c('ok', 'partial')) &
            (proj$group %in% md$group[!(md$status %in% c('absent', 'rearranged'))])
    idx.parent <- match(paste(proj$src.acc, proj$src.parent), paste(proj$src.acc, proj$src.row))
    for(i in 1:5){
      drop <- keep & !is.na(idx.parent) & !keep[idx.parent]
      if(!any(drop)) break
      keep[drop] <- FALSE
    }
    if(!any(keep)) next

    idx.new <- cumsum(keep)
    idx.new[!keep] <- NA
    gff.p <- data.frame(V1 = s.name, V2 = 'pannagram', V3 = proj$V3,
                        V4 = proj$beg, V5 = proj$end, V6 = proj$V6,
                        V7 = proj$strand, V8 = proj$V8, V9 = proj$V9,
                        id = proj$id, parent.idx = idx.new[idx.parent],
                        stringsAsFactors = FALSE)[keep, , drop = FALSE]
    gff.p$root <- idx.new[match(paste(proj$src.acc, proj$src.root), paste(proj$src.acc, proj$src.row))][keep]
    is.root <- proj$is.root[keep]
    i.md <- match(proj$group[keep], md$group)
    root.id <- ifelse(is.root, paste0(groups$group[proj$group[keep]], '.', acc.t), NA)
    root.attr <- annoAttrJoin(paste0('annogroup=', groups$group[proj$group[keep]]),
                              paste0('projected_from=', proj$src.acc[keep]),
                              ifelse(is.na(gff.p$id), NA, paste0('orig_id=', gff.p$id)),
                              paste0('status=', md$status[i.md]),
                              paste0('cov=', md$cov[i.md]))
    gff.p <- annoRename(gff.p, root.id, root.attr)
    gff.p$ord.chr <- q.chr
    gff.p$ord.beg <- gff.p$V4[gff.p$root]
    gff.p$ord.root <- paste('p', s.comb, gff.p$root)
    proj.gff.all[[acc.t]][[s.comb]] <- gff.p[, c(paste0('V', 1:9), 'ord.chr', 'ord.beg', 'ord.root')]
  }

  n.proj <- table(factor(unlist(lapply(proj.tbl.all[paste(s.comb, names(proj.list))],
                                       function(x) x$group[!(x$status %in% c('absent', 'rearranged'))])),
                         levels = groups$group))
  groups$n.projected <- as.integer(n.proj)
  groups.all[[s.comb]] <- groups

  pokaz('  groups:', n.gr, '; in clusters of several groups:', sum(groups$n.split > 1),
        '; with several genes per accession:', sum(groups$n.acc.multi > 0),
        '; covered by genes of other groups:', sum(groups$n.shared > 0))
  rm(ivl.list, proj.list, cds.src, cds.trg)
  rhdf5::h5closeAll()
  gc()
}

if(length(groups.all) == 0){
  pokazAttention('No annotation groups were found.')
  quit(save = "no", status = 0)
}

# ***********************************************************************
# ---- Tables ----

writeTable <- function(x, file){
  write.table(x, file, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
}

groups <- do.call(rbind, groups.all)
groups$len <- groups$end - groups$beg + 1
writeTable(groups[, c('group', 'chr', 'strand', 'beg', 'end', 'len', 'n.acc', 'n.genes',
                      'n.acc.multi', 'n.shared', 'n.split', 'n.projected', 'types')],
           paste0(path.res, 'annogroups.txt'))

writeTable(do.call(rbind, corresp.all), paste0(path.res, 'correspondence.txt'))

lost <- do.call(rbind, lapply(acc.gff, function(acc){
  gff <- gff.list[[acc]]
  idx <- which(gff$is.root & gff$use & (gff$reason != ''))
  data.frame(acc = rep(acc, length(idx)), id = ifelse(is.na(gff$id[idx]), '-', gff$id[idx]),
             type = gff$V3[idx], chr = gff$V1[idx], beg = gff$V4[idx], end = gff$V5[idx],
             strand = gff$V7[idx], reason = gff$reason[idx], stringsAsFactors = FALSE)
}))
writeTable(lost, paste0(path.res, 'lost_genes.txt'))

proj.tbl <- do.call(rbind, proj.tbl.all)
if(is.null(proj.tbl)){
  proj.tbl <- data.frame(group = character(0), acc = character(0), src.acc = character(0),
                         src.id = character(0), chr = character(0), beg = integer(0), end = integer(0),
                         strand = character(0), status = character(0), cov = numeric(0))
}
writeTable(proj.tbl, paste0(path.res, 'projection.txt'))

pokaz('Groups:', nrow(groups), '; lost genes:', nrow(lost))
if(nrow(lost) > 0) pokaz('  lost by reason:', paste(names(table(lost$reason)), table(lost$reason), collapse = ', '))
if(nrow(proj.tbl) > 0) pokaz('  projected models by status:', paste(names(table(proj.tbl$status)), table(proj.tbl$status), collapse = ', '))

# ***********************************************************************
# ---- GFF of the accessions and of the pangenome ----

pangen.all <- list()
for(acc in accessions){

  gff.own <- NULL
  if(acc %in% acc.gff){
    gff <- gff.list[[acc]]
    root.attr <- annoAttrJoin(ifelse(is.na(gff$groups), NA, paste0('annogroup=', gff$groups)),
                              ifelse(is.na(gff$id), NA, paste0('orig_id=', gff$id)))
    gff <- annoRename(gff, gff$root.id, root.attr)

    # Pangenome coordinates of the grouped genes
    idx <- which(!is.na(gff$root.id[gff$root]) & (gff$reason == ''))
    if(length(idx) > 0){
      gff.pan <- gff[idx, paste0('V', 1:9)]
      gff.pan$V1 <- paste0('PanGen_Chr', gff$pan.chr[idx])
      gff.pan$V4 <- gff$pan.beg[idx]
      gff.pan$V5 <- gff$pan.end[idx]
      gff.pan$V7 <- annoStrand(gff$V7[idx], gff$delta[idx])
      is.r <- gff$is.root[idx]
      gff.pan$V9[is.r] <- annoAttrJoin(paste0('ID=', gff$id.new[idx[is.r]]),
                                       paste0('Parent=', gff$group[idx[is.r]]),
                                       annoAttrDrop(gff.pan$V9[is.r], 'ID'))
      gff.pan$ord.group <- gff$group.ord[gff$root[idx]]
      gff.pan$ord.acc <- match(acc, accessions)
      pangen.all[[acc]] <- gff.pan
    }

    gff.own <- gff[, paste0('V', 1:9)]
    gff.own$ord.chr <- gff$chr
    gff.own$ord.beg <- gff$V4[gff$root]
    gff.own$ord.root <- paste('n', gff$root)
  }

  gff.proj <- if(is.null(proj.gff.all[[acc]])) NULL else do.call(rbind, proj.gff.all[[acc]])
  gff.acc <- rbind(gff.own, gff.proj)
  if(is.null(gff.acc) || (nrow(gff.acc) == 0)) next
  gff.acc <- gff.acc[order(gff.acc$ord.chr, gff.acc$ord.beg, gff.acc$ord.root), , drop = FALSE]
  writeGFF(gff.acc, paste0(path.gff, acc, '.gff'))
}

gff.groups <- data.frame(V1 = paste0('PanGen_Chr', groups$chr), V2 = 'pannagram', V3 = 'annogroup',
                         V4 = groups$beg, V5 = groups$end, V6 = '.', V7 = groups$strand, V8 = '.',
                         V9 = paste0('ID=', groups$group, ';n_acc=', groups$n.acc, ';n_genes=', groups$n.genes),
                         ord.group = groups$ord, ord.acc = 0, stringsAsFactors = FALSE)
gff.pan <- rbind(gff.groups, do.call(rbind, pangen.all))
gff.pan <- gff.pan[order(gff.pan$ord.group, gff.pan$ord.acc), , drop = FALSE]
writeGFF(gff.pan, paste0(path.res, 'pangen.gff'))

pokaz('Results are in', path.res)
