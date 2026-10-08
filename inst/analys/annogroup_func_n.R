# ============================================================================
#  Annotation groups: functions behind analys_05_annogroups_n.R
#
#  INTERVAL FORMAT (_n): the alignment is read as interval tables
#  (see utils/interval_func.R), nothing is expanded to the pangenome length.
#
#  The file is sourced by the script and is not a part of the library.
# ============================================================================

# Top-level GFF records that describe a sequence and not a feature
anno.skip.types <- c('region', 'chromosome', 'contig', 'scaffold', 'supercontig', 'biological_region')

anno.stops  <- c('TAA', 'TAG', 'TGA')
anno.starts <- c('ATG', 'GTG', 'TTG')

# ---------------------------------------------------------------------------
# ---- GFF ------------------------------------------------------------------
# ---------------------------------------------------------------------------

#' Read a GFF file
#'
#' Tolerant to '#' and quotes in the attributes, to comment lines and to the
#' '##FASTA' section. Returns the 9 columns V1..V9, coordinates are integers.
annoReadGFF <- function(file, skip.types = anno.skip.types){
  x <- readLines(file, warn = FALSE)
  i.fasta <- which(startsWith(x, '##FASTA'))
  if(length(i.fasta) > 0) x <- x[seq_len(i.fasta[1] - 1)]
  x <- x[!startsWith(x, '#') & (nchar(x) > 0)]
  x <- strsplit(x, '\t', fixed = TRUE)
  x <- x[lengths(x) == 9]

  m <- matrix(as.character(unlist(x, use.names = FALSE)), ncol = 9, byrow = TRUE)
  gff <- data.frame(V1 = m[,1], V2 = m[,2], V3 = m[,3],
                    V4 = suppressWarnings(as.integer(m[,4])),
                    V5 = suppressWarnings(as.integer(m[,5])),
                    V6 = m[,6], V7 = m[,7], V8 = m[,8], V9 = m[,9],
                    stringsAsFactors = FALSE)
  gff <- gff[!is.na(gff$V4) & !is.na(gff$V5) & (gff$V4 > 0) & (gff$V5 > 0), , drop = FALSE]
  gff <- gff[!(tolower(gff$V3) %in% skip.types), , drop = FALSE]

  idx.swap <- which(gff$V4 > gff$V5)
  if(length(idx.swap) > 0){
    tmp <- gff$V4[idx.swap]
    gff$V4[idx.swap] <- gff$V5[idx.swap]
    gff$V5[idx.swap] <- tmp
  }
  rownames(gff) <- NULL
  gff
}

#' Value of an attribute of the 9th GFF column (NA if absent)
annoAttr <- function(v9, key){
  res <- rep(NA_character_, length(v9))
  m <- regexpr(paste0('(^|;)\\s*', key, '=[^;]*'), v9, perl = TRUE)
  ok <- m > 0
  res[ok] <- sub(paste0('^;?\\s*', key, '='), '', regmatches(v9, m), perl = TRUE)
  res
}

#' Remove attributes from the 9th GFF column
annoAttrDrop <- function(v9, keys){
  v9 <- gsub(paste0('(^|;)\\s*(', paste(keys, collapse = '|'), ')=[^;]*'), '', v9, perl = TRUE)
  sub('^;+', '', v9)
}

#' Join attribute strings with ';', empty and NA ones are skipped
annoAttrJoin <- function(...){
  parts <- list(...)
  res <- NULL
  for(p in parts){
    p[is.na(p)] <- ''
    res <- if(is.null(res)) p else paste0(res, ifelse((res != '') & (p != ''), ';', ''), p)
  }
  res
}

#' Parent-child structure of a GFF
#'
#' Adds: `id`, `parent` (the first one), `parent.idx` (row of the parent or NA)
#' and `root` (row of the top-level feature). A feature without a parent is its
#' own root, whatever its type is: `gene` in most annotations, `CDS`/`tRNA` in Prokka.
annoHierarchy <- function(gff){
  n <- nrow(gff)
  gff$id <- annoAttr(gff$V9, 'ID')
  gff$parent <- sub(',.*', '', annoAttr(gff$V9, 'Parent'))
  pidx <- match(gff$parent, gff$id, incomparables = NA)
  pidx[!is.na(pidx) & (pidx == seq_len(n))] <- NA
  gff$parent.idx <- pidx

  r <- ifelse(is.na(pidx), seq_len(n), pidx)
  for(i in 1:20){
    r.next <- r[r]
    if(all(r.next == r)) break
    r <- r.next
  }
  gff$root <- r
  gff
}

#' Chromosome number from the sequence name: "<anything><s.chr><number>"
annoChr <- function(seqid, s.chr){
  chr <- rep(NA_integer_, length(seqid))
  ok <- grepl(paste0(s.chr, '\\d+$'), seqid)
  chr[ok] <- as.integer(sub(paste0('.*', s.chr, '(\\d+)$'), '\\1', seqid[ok]))
  chr
}

# ---------------------------------------------------------------------------
# ---- Accession <-> pangenome ----------------------------------------------
# ---------------------------------------------------------------------------

#' Interval table with the columns needed for the lookups below
annoIvlPrep <- function(ivl){
  if((nrow(ivl) > 0) && all(ivl$block == 0)){
    ivl <- ivl[order(ivl$pan.beg), , drop = FALSE]
    ivl$block <- as.integer(ivlBlocks(ivl))
  }
  acc.end <- ivl$acc.beg + ivl$delta * (ivl$len - 1L)
  ivl$lo <- pmin(ivl$acc.beg, acc.end)
  ivl$hi <- pmax(ivl$acc.beg, acc.end)
  ivl$pan.end <- ivl$pan.beg + ivl$len - 1L
  ivl
}

#' Features of an accession -> pangenome coordinates
#'
#' @param ivl Interval table after annoIvlPrep().
#' @param beg,end Coordinates in the accession.
#' @param snap If TRUE, an unaligned end is moved inwards to the nearest aligned position.
#' @return data.frame: pan.beg, pan.end, delta (+1/-1, strand of the block),
#'   snapped, reason ('' - lifted, 'unaligned', 'diff_blocks').
annoLift <- function(ivl, beg, end, snap = FALSE){
  n <- length(beg)
  res <- data.frame(pan.beg = integer(n), pan.end = integer(n), delta = integer(n),
                    snapped = logical(n), reason = rep('unaligned', n),
                    stringsAsFactors = FALSE)
  if((n == 0) || (nrow(ivl) == 0)) return(res)

  a <- ivl[order(ivl$lo), , drop = FALSE]
  m <- nrow(a)
  i1 <- findInterval(beg, a$lo)
  i2 <- findInterval(end, a$lo)
  in1 <- i1 > 0; in1[in1] <- beg[in1] <= a$hi[i1[in1]]
  in2 <- i2 > 0; in2[in2] <- end[in2] <= a$hi[i2[in2]]

  q1 <- beg; q2 <- end
  ok <- in1 & in2
  if(snap){
    s1 <- !in1 & (i1 < m)
    i1[s1] <- i1[s1] + 1L
    q1[s1] <- a$lo[i1[s1]]
    s2 <- !in2 & (i2 > 0)
    q2[s2] <- a$hi[i2[s2]]
    ok <- (in1 | s1) & (in2 | s2) & (q1 <= q2)
    res$snapped <- ok & !(in1 & in2)
  }

  k <- which(ok)
  if(length(k) == 0) return(res)
  j1 <- i1[k]; j2 <- i2[k]
  p1 <- a$pan.beg[j1] + a$delta[j1] * (q1[k] - a$acc.beg[j1])
  p2 <- a$pan.beg[j2] + a$delta[j2] * (q2[k] - a$acc.beg[j2])
  same <- (a$block[j1] == a$block[j2]) & (a$delta[j1] == a$delta[j2])

  res$pan.beg[k] <- pmin(p1, p2)
  res$pan.end[k] <- pmax(p1, p2)
  res$delta[k] <- a$delta[j1]
  res$reason[k] <- ifelse(same, '', 'diff_blocks')
  res
}

#' Pangenome intervals -> coordinates of an accession
#'
#' @param ivl Interval table after annoIvlPrep().
#' @param p1,p2 Begin and end in the pangenome.
#' @return data.frame: beg, end, delta, n.cov (number of positions of the accession
#'   aligned within the interval) and state:
#'   'ok' - both ends are aligned, 'partial' - an end is moved inwards,
#'   'absent' - nothing is aligned, 'rearranged' - the ends are in different blocks.
annoProject <- function(ivl, p1, p2){
  n <- length(p1)
  res <- data.frame(beg = integer(n), end = integer(n), delta = integer(n),
                    n.cov = numeric(n), state = rep('absent', n), stringsAsFactors = FALSE)
  if((n == 0) || (nrow(ivl) == 0)) return(res)

  a <- ivl[order(ivl$pan.beg), , drop = FALSE]
  m <- nrow(a)
  cum <- c(0, cumsum(as.numeric(a$len)))   # aligned positions before the interval

  i1 <- findInterval(p1, a$pan.beg)
  i2 <- findInterval(p2, a$pan.beg)
  in1 <- i1 > 0; in1[in1] <- p1[in1] <= a$pan.end[i1[in1]]

  # The first aligned position >= p1 and the last one <= p2
  f <- p1
  s1 <- !in1
  i1[s1] <- i1[s1] + 1L
  ok <- (i1 <= m) & (i2 > 0)
  f[s1 & ok] <- a$pan.beg[i1[s1 & ok]]
  l <- p2
  l[ok] <- pmin(p2[ok], a$pan.end[i2[ok]])
  ok[ok] <- f[ok] <= l[ok]

  k <- which(ok)
  if(length(k) == 0) return(res)
  j1 <- i1[k]; j2 <- i2[k]
  t1 <- a$acc.beg[j1] + a$delta[j1] * (f[k] - a$pan.beg[j1])
  t2 <- a$acc.beg[j2] + a$delta[j2] * (l[k] - a$pan.beg[j2])
  same <- (a$block[j1] == a$block[j2]) & (a$delta[j1] == a$delta[j2])

  res$beg[k] <- pmin(t1, t2)
  res$end[k] <- pmax(t1, t2)
  res$delta[k] <- a$delta[j1]
  res$n.cov[k] <- (cum[j2] + (l[k] - a$pan.beg[j2] + 1)) - (cum[j1] + (f[k] - a$pan.beg[j1]))
  res$state[k] <- ifelse(!same, 'rearranged',
                         ifelse((f[k] == p1[k]) & (l[k] == p2[k]), 'ok', 'partial'))
  res
}

anno.strand.flip <- c('+' = '-', '-' = '+')

#' Strand after going through a block with the orientation `delta`
annoStrand <- function(strand, delta){
  flip <- (delta < 0) & (strand %in% names(anno.strand.flip))
  strand[flip] <- anno.strand.flip[strand[flip]]
  strand
}

# ---------------------------------------------------------------------------
# ---- Groups ---------------------------------------------------------------
# ---------------------------------------------------------------------------

#' Annotation groups of one strand of one pangenome chromosome
#'
#' Two genes are the same gene if they overlap by at least `cutoff` of the longer
#' one; a group is a connected component of this relation. So neighbouring genes
#' that overlap by a few nucleotides stay apart, and a gene annotated as one
#' model in some accessions and as two in the others does not glue the two.
#'
#' Most accessions carry identical intervals, therefore the comparison is done
#' on the unique intervals of a cluster of overlapping genes, not on the genes.
#'
#' @param beg,end Pangenome coordinates of the genes.
#' @param cutoff Minimal overlap, as a share of the longer gene.
#' @return list: `group` - index of the group of every gene;
#'   `groups` - data.frame: beg, end (span of the group), cluster,
#'   sorted along the chromosome.
annoGroups <- function(beg, end, cutoff = 0.5){
  n <- length(beg)
  if(n == 0){
    return(list(group = integer(0),
                groups = data.frame(beg = integer(0), end = integer(0), cluster = integer(0))))
  }

  # Unique intervals, sorted
  s.key <- paste(beg, end)
  is.u <- !duplicated(s.key)
  o <- order(beg[is.u], end[is.u])
  u.beg <- beg[is.u][o]
  u.end <- end[is.u][o]
  gene.u <- match(s.key, s.key[is.u][o])
  m <- length(u.beg)

  # Clusters of overlapping intervals
  cm <- cummax(u.end)
  u.cl <- cumsum(c(TRUE, u.beg[-1] > cm[-m]))

  # Components within the clusters
  u.lab <- rep(1L, m)
  n.u <- tabulate(u.cl)
  for(idx in split(seq_len(m), u.cl)[n.u > 1]){
    b <- u.beg[idx]; e <- u.end[idx]; len <- e - b + 1
    ov <- outer(e, e, pmin) - outer(b, b, pmax) + 1
    adj <- (ov / outer(len, len, pmax)) >= cutoff
    lab <- seq_along(idx)
    repeat{
      lab.new <- apply(adj, 1, function(x) min(lab[x]))
      if(all(lab.new == lab)) break
      lab <- lab.new
    }
    u.lab[idx] <- lab
  }

  # Groups along the chromosome
  s.gr <- paste(u.cl, u.lab)
  is.first <- !duplicated(s.gr)
  u.gr <- match(s.gr, s.gr[is.first])
  groups <- data.frame(beg = u.beg[is.first],
                       end = as.integer(tapply(u.end, u.gr, max)),
                       cluster = u.cl[is.first])
  ord <- order(groups$beg, groups$end)
  groups <- groups[ord, , drop = FALSE]
  rownames(groups) <- NULL
  list(group = match(u.gr, ord)[gene.u], groups = groups)
}

#' Genes that cover groups of other genes
#'
#' A gene annotated as one model where the other accessions have two (or the
#' other way round) belongs to one group, but occupies the place of another one.
#' Such a gene is reported as a secondary member of that group: it is shown in
#' the correspondence table, and nothing is projected to its accession.
#'
#' @param beg,end Pangenome coordinates of the genes.
#' @param group Group of every gene, annoGroups().
#' @param groups Table of groups, annoGroups().
#' @param min.overlap Overlap, as a share of the gene or of the group.
#' @return data.frame: gene (index), group (index), primary.
annoAssign <- function(beg, end, group, groups, min.overlap = 0.5){
  res <- data.frame(gene = seq_along(beg), group = group, primary = rep(TRUE, length(beg)))

  # Only the clusters with several groups can have such genes
  cl <- groups$cluster
  idx.gr <- which(cl %in% cl[duplicated(cl)])
  if(length(idx.gr) == 0) return(res)
  idx.gene <- which(group %in% idx.gr)

  res.sec <- list()
  gr.cl <- split(idx.gr, cl[idx.gr])
  gene.cl <- split(idx.gene, cl[group[idx.gene]])
  for(s.cl in names(gr.cl)){
    g <- gr.cl[[s.cl]]
    i <- gene.cl[[s.cl]]
    ov <- pmax(outer(end[i], groups$end[g], pmin) - outer(beg[i], groups$beg[g], pmax) + 1, 0)
    is.over <- (ov >= min.overlap * (end[i] - beg[i] + 1)) |
               t(t(ov) >= min.overlap * (groups$end[g] - groups$beg[g] + 1))
    is.over[cbind(seq_along(i), match(group[i], g))] <- FALSE
    k <- which(is.over, arr.ind = TRUE)
    if(nrow(k) == 0) next
    res.sec[[s.cl]] <- data.frame(gene = i[k[, 1]], group = g[k[, 2]], primary = FALSE)
  }
  if(length(res.sec) > 0) res <- rbind(res, do.call(rbind, res.sec))
  rownames(res) <- NULL
  res
}

# ---------------------------------------------------------------------------
# ---- Coding sequences -----------------------------------------------------
# ---------------------------------------------------------------------------

annoRevCompl <- function(s){
  s <- chartr('ACGTRYKMBDHVN', 'TGCAYRMKVHDBN', s)
  vapply(s, function(x) intToUtf8(rev(utf8ToInt(x))), '', USE.NAMES = FALSE)
}

#' Properties of the coding sequences of transcripts
#'
#' @param s.chr Chromosome as one string.
#' @param beg,end,strand,phase CDS pieces.
#' @param tr Transcript of every piece.
#' @return data.frame, a row per transcript: tr, len (without the phase offset),
#'   first and last codons, stop.inner (a stop codon before the last codon).
annoCdsProps <- function(s.chr, beg, end, strand, phase, tr){
  if(length(beg) == 0){
    return(data.frame(tr = character(0), len = integer(0), first = character(0),
                      last = character(0), stop.inner = logical(0), stringsAsFactors = FALSE))
  }
  tr <- as.character(tr)
  o <- order(tr, beg)
  beg <- beg[o]; end <- end[o]; strand <- strand[o]; tr <- tr[o]
  phase <- suppressWarnings(as.integer(phase[o]))
  phase[is.na(phase)] <- 0L

  pieces <- toupper(substring(s.chr, beg, end))
  i.first <- which(!duplicated(tr))
  cds <- vapply(split(pieces, tr), paste, '', collapse = '')[tr[i.first]]
  i.last <- c(i.first[-1] - 1L, length(tr))
  is.minus <- strand[i.first] == '-'
  cds[is.minus] <- annoRevCompl(cds[is.minus])

  # The phase of the piece the translation starts from
  ph <- ifelse(is.minus, phase[i.last], phase[i.first])
  cds <- substring(cds, ph + 1L)

  n.cod <- nchar(cds) %/% 3L
  stop.inner <- vapply(seq_along(cds), function(i){
    if(n.cod[i] < 2) return(FALSE)
    p <- seq(1L, by = 3L, length.out = n.cod[i] - 1L)
    any(substring(cds[i], p, p + 2L) %in% anno.stops)
  }, logical(1))

  data.frame(tr = tr[i.first], len = nchar(cds),
             first = substr(cds, 1, 3),
             last = substr(cds, 3 * n.cod - 2, 3 * n.cod),
             stop.inner = stop.inner, stringsAsFactors = FALSE)
}

#' Status of a projected coding model against the model it was projected from
#'
#' @param src,trg Rows of annoCdsProps() for the same transcripts.
#' @return Character vector: 'intact', 'frameshift', 'stop', 'start_lost', 'stop_lost'.
annoCdsStatus <- function(src, trg){
  status <- rep('intact', nrow(src))
  status[(src$last %in% anno.stops) & !(trg$last %in% anno.stops)] <- 'stop_lost'
  status[(src$first %in% anno.starts) & !(trg$first %in% anno.starts)] <- 'start_lost'
  status[trg$stop.inner] <- 'stop'
  status[((trg$len - src$len) %% 3) != 0] <- 'frameshift'
  status
}

# The best transcript defines the status of a gene
anno.status.rank <- c('intact', 'unchecked', 'stop_lost', 'start_lost', 'stop', 'frameshift', 'truncated')

# ---------------------------------------------------------------------------
# ---- Common IDs -----------------------------------------------------------
# ---------------------------------------------------------------------------

#' Running number of an element within its key, in the order of appearance
annoRunNum <- function(key){
  o <- order(key)
  num <- integer(length(key))
  num[o] <- sequence(rle(key[o])$lengths)
  num
}

#' Rename the features of the grouped genes
#'
#' The root gets the ID "<group>.<accession>", its descendants -
#' "<root ID>.<type><number>"; `Parent` is rewritten accordingly. The pieces of
#' one feature that share an ID (multi-line CDS) keep sharing it.
#'
#' @param gff GFF with the columns of annoHierarchy().
#' @param root.id New IDs of the roots: a vector along the rows of `gff`,
#'   NA for the rows whose root is not renamed.
#' @param root.attr Attributes to add to the renamed roots (a vector along the rows).
#' @return `gff` with the new V9 and the column `id.new`.
annoRename <- function(gff, root.id, root.attr = NULL){
  n <- nrow(gff)
  gff$id.new <- NA_character_
  if(n == 0) return(gff)
  is.root <- gff$root == seq_len(n)
  id.root <- root.id[gff$root]
  idx <- which(!is.na(id.root))
  if(length(idx) == 0) return(gff)

  # Roots
  idx.r <- idx[is.root[idx]]
  gff$id.new[idx.r] <- id.root[idx.r]

  # Descendants with an ID: numbered within the root and the type
  idx.c <- idx[!is.root[idx] & !is.na(gff$id[idx])]
  if(length(idx.c) > 0){
    s.set <- paste(gff$root[idx.c], gff$V3[idx.c], sep = '|')
    s.key <- paste(s.set, gff$id[idx.c], sep = '|')
    is.new <- !duplicated(s.key)
    num <- annoRunNum(s.set[is.new])[match(s.key, s.key[is.new])]
    gff$id.new[idx.c] <- paste0(id.root[idx.c], '.', gff$V3[idx.c], num)
  }

  # Parents
  parent.new <- gff$id.new[gff$parent.idx]
  s.parent <- annoAttr(gff$V9[idx], 'Parent')
  idx.multi <- which(grepl(',', s.parent, fixed = TRUE))
  for(i in idx.multi){
    p <- strsplit(s.parent[i], ',', fixed = TRUE)[[1]]
    p.new <- gff$id.new[match(p, gff$id)]
    p.new[is.na(p.new)] <- p[is.na(p.new)]
    parent.new[idx[i]] <- paste(p.new, collapse = ',')
  }
  k <- is.na(parent.new[idx])  # a parent that is not in the file: the reference is kept as it is
  parent.new[idx[k]] <- s.parent[k]

  s.attr <- if(is.null(root.attr)) rep('', n) else root.attr
  s.attr[!is.root] <- ''
  gff$V9[idx] <- annoAttrJoin(ifelse(is.na(gff$id.new[idx]), NA, paste0('ID=', gff$id.new[idx])),
                              ifelse(is.na(parent.new[idx]), NA, paste0('Parent=', parent.new[idx])),
                              annoAttrDrop(gff$V9[idx], c('ID', 'Parent')),
                              s.attr[idx])
  gff
}
