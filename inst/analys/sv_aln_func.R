# Alignment of the members of an SV family and the consensus of the family.
#
# The orientation of an SV in the pangenome is arbitrary, and its boundaries do
# not have to coincide with the boundaries of the mobile element: some copies
# carry a piece of the flanking sequence, some are truncated. So the members are
# first brought to one strand, then aligned, and then the alignment is trimmed by
# the coverage of its columns - the step of the coverage shows where the element
# actually begins and ends.


#' Bring the members of a family to one strand by k-mers
#'
#' A fallback for \code{\link{getComponentSequences}}, which is the main way to
#' orient the members: it takes the strand from the nestedness table, where it has
#' already been computed by the similarity search, and propagates it along the
#' edges of the component. This function needs nothing but the sequences
#' themselves, and is used when the nestedness table is not available or when the
#' members of the family are not connected in it.
#'
#' The orientation is decided by the number of shared k-mers with the reference,
#' counted over the whole sequence and not over the termini: truncated copies
#' may have no termini at all.
#'
#' @param seqs Named vector of sequences.
#' @param ref A single reference sequence; by default the longest member.
#' @param k Length of the k-mer.
#'
#' @return A list with `seqs` (the sequences on one strand) and `strand`.
#'
#' @author Anna A. Igolkina
#' @export
svFamilyOrient <- function(seqs, ref = NULL, k = 12){

  seqs = toupper(seqs)
  if(is.null(ref)) ref = seqs[[which.max(nchar(seqs))]]
  ref = toupper(ref)

  kmers <- function(s){
    n = nchar(s)
    if(n < k) return(character(0))
    unique(substring(s, 1:(n - k + 1), k:n))
  }
  km.ref = kmers(ref)

  strand = rep('+', length(seqs))
  for(i in seq_along(seqs)){
    s.rc = revCompl(seqs[[i]])
    if(sum(kmers(s.rc) %in% km.ref) > sum(kmers(seqs[[i]]) %in% km.ref)){
      seqs[[i]] = s.rc
      strand[i] = '-'
    }
  }
  names(strand) = names(seqs)

  return(list(seqs = seqs, strand = strand))
}


#' Align a set of sequences with MAFFT
#'
#' @param seqs Named vector of sequences.
#' @param path.work Directory for the temporary files.
#' @param args Arguments of MAFFT. The default gap penalties are stricter than
#' those of MAFFT itself: with a cheap gap the members of a family get torn into
#' pieces scattered over the whole alignment, because a distant piece can always
#' add a few matches. On a family of 16 members the default `--op 1.53` gives 241
#' blocks of sequence, `--op 5` gives 166 and `--op 20` gives 101, but the last
#' one starts to eat the edges of the alignment.
#' @param bin Name or path of the MAFFT executable.
#'
#' @return A character matrix of the alignment, one row per sequence.
#'
#' @author Anna A. Igolkina
#' @export
svAlnMafft <- function(seqs, path.work = tempdir(), args = '--op 5 --ep 0.2 --quiet', bin = 'mafft'){

  if(length(seqs) < 2) stop('At least two sequences are needed for the alignment')
  if(is.null(names(seqs))) names(seqs) = paste0('seq_', seq_along(seqs))

  pref = tempfile(pattern = 'sv_aln_', tmpdir = path.work)
  file.in  = paste0(pref, '_in.fasta')
  file.out = paste0(pref, '_out.fasta')
  on.exit(unlink(c(file.in, file.out)), add = TRUE)

  writeFasta(seqs, file.in)
  status = system(paste(bin, args, shQuote(file.in), '>', shQuote(file.out)))
  if((status != 0) || !file.exists(file.out) || (file.size(file.out) == 0)){
    stop(paste0('MAFFT has failed, check that "', bin, '" is available'))
  }

  mx = aln2mx(readFasta(file.out))
  return(mx[names(seqs), , drop = FALSE])
}


#' Coverage of the columns of an alignment
#'
#' A sequence covers everything between its own beginning and its own end,
#' the gaps inside it do not break the coverage: an internal deletion is a part
#' of the sequence, while a gap at the edge only means that the sequence
#' does not reach that far.
#'
#' @param mx A character matrix of the alignment.
#'
#' @return An integer vector: the number of sequences covering every column.
#'
#' @author Anna A. Igolkina
#' @export
svAlnCover <- function(mx){

  cov = integer(ncol(mx))
  for(i in 1:nrow(mx)){
    idx = which(mx[i, ] != '-')
    if(length(idx) == 0) next
    span = idx[1]:idx[length(idx)]
    cov[span] = cov[span] + 1L
  }

  return(cov)
}


#' Mask of the columns covered by every sequence
#'
#' @param mx A character matrix of the alignment.
#'
#' @return A logical matrix of the size of `mx`: `TRUE` where the sequence
#' reaches the column, see \code{\link{svAlnCover}}.
#'
#' @author Anna A. Igolkina
#' @export
svAlnSpan <- function(mx){

  span = matrix(FALSE, nrow(mx), ncol(mx), dimnames = dimnames(mx))
  for(i in 1:nrow(mx)){
    idx = which(mx[i, ] != '-')
    if(length(idx) == 0) next
    span[i, idx[1]:idx[length(idx)]] = TRUE
  }

  return(span)
}


#' Agreement of the sequences in the columns of an alignment
#'
#' Only the sequences that cover the column are taken into account,
#' see \code{\link{svAlnCover}}.
#'
#' The diversity of a column is the share of its dominant *nucleotide*: an internal
#' gap cannot win the column, it only tells that this sequence has a deletion here.
#'
#' @param mx A character matrix of the alignment.
#' @param use.gap Let an internal gap compete for the column together with the
#' nucleotides; `FALSE` by default.
#'
#' @return A numeric vector: the share of the dominant nucleotide in every column.
#'
#' @author Anna A. Igolkina
#' @export
svAlnAgree <- function(mx, use.gap = FALSE){

  mx = toupper(mx)
  mx[!svAlnSpan(mx)] = ''                   # outside of the sequence: not a variant

  cnt = mx2profile(mx, gap.flag = use.gap)  # counts of A, C, G, T (and of the gap)
  n = colSums(cnt)
  agree = ifelse(n > 0, apply(cnt, 2, max) / n, 0)

  return(agree)
}


#' Trim the edges of an alignment
#'
#' The edges of an alignment of the members of a family are ragged: the boundaries
#' of an SV do not coincide with the boundaries of the element, so the first and
#' the last columns are covered by a few members and disagree with each other.
#' Both edges are trimmed inwards until a position that is covered by at least
#' `min.cov` sequences and starts a window of `wnd` columns where every column
#' agrees at the level of `min.agree`.
#'
#' @param mx A character matrix of the alignment.
#' @param min.cov Minimum number of sequences covering a column.
#' @param wnd Size of the window that has to be clean.
#' @param min.agree Minimum agreement in every column of the window.
#' @param use.gap Parameter of \code{\link{svAlnAgree}}.
#'
#' @return A list with the trimmed matrix `mx`, the profiles `cov` and `agree`
#' of the initial alignment, and the boundaries `beg` and `end` of the trimmed part.
#'
#' @author Anna A. Igolkina
#' @export
svAlnTrim <- function(mx, min.cov = 3, wnd = 5, min.agree = 0.9, use.gap = FALSE){

  cov = svAlnCover(mx)
  agree = svAlnAgree(mx, use.gap = use.gap)

  # an alignment of two sequences cannot be covered by three of them
  min.cov = min(min.cov, nrow(mx))

  res = list(mx = NULL, cov = cov, agree = agree, beg = NA, end = NA)

  n = ncol(mx)
  if(n < wnd) return(res)

  good = (cov >= min.cov) & (agree >= min.agree)

  # windows of wnd columns, all of them good
  cs = cumsum(c(0, as.integer(good)))
  wnd.good = which((cs[(wnd + 1):(n + 1)] - cs[1:(n - wnd + 1)]) == wnd)
  if(length(wnd.good) == 0) return(res)

  res$beg = wnd.good[1]
  res$end = wnd.good[length(wnd.good)] + wnd - 1
  res$mx = mx[, res$beg:res$end, drop = FALSE]

  return(res)
}


#' Consensus of an alignment
#'
#' In every column the most represented symbol among the sequences that cover
#' the column wins, see \code{\link{svAlnCover}}. The gap votes together with
#' the nucleotides: if most of the covering sequences have a deletion here, the
#' column is a deletion, and it leaves the consensus afterwards. The sequences
#' that do not reach the column do not vote at all - a leading or a trailing gap
#' says nothing about what should be in the middle of the element.
#'
#' @param mx A character matrix of the alignment.
#'
#' @return A list with the consensus `seq` (a string, gaps removed), `cons`
#' (by columns, gaps kept), the coverage `cov` and the fraction `freq` of the
#' winning symbol in every column.
#'
#' @author Anna A. Igolkina
#' @export
svAlnCons <- function(mx){

  mx = toupper(mx)
  s.val = c('A', 'C', 'G', 'T', '-')

  mx[!svAlnSpan(mx)] = ''                   # only the covering sequences vote
  cnt = mx2profile(mx, gap.flag = TRUE)

  n = colSums(cnt)
  i.max = apply(cnt, 2, which.max)
  cons = ifelse(n > 0, s.val[i.max], '-')
  freq = ifelse(n > 0, apply(cnt, 2, max) / n, 0)

  return(list(seq  = nt2seq(cons[cons != '-']),
              cons = cons,
              cov  = n / nrow(mx),
              freq = freq))
}


#' Choose the representatives of a family
#'
#' Only the members connected to another member by an edge in both directions are
#' taken: a one-sided edge means that one sequence is nested into the other, that
#' is, they share a fragment rather than the whole element, and a set of such
#' members does not align. Among those, the longest ones are taken.
#'
#' Length alone is a bad criterion: it does not agree with how well a member is
#' connected in the graph of the family, and the longest SVs are often composite.
#' Mutual edges cut those off, and the length then picks the most complete copies
#' among the remaining ones.
#'
#' @param len Lengths of the members of the family.
#' @param n.repr Number of representatives to take.
#' @param mutual Logical vector in the order of `len`: does the member have an
#' edge in both directions; `NULL` if the graph is not available.
#'
#' @return Indices of the chosen members in `len`, sorted by decreasing length.
#'
#' @author Anna A. Igolkina
#' @export
svFamilyRepr <- function(len, n.repr = 20, mutual = NULL){

  keep = seq_along(len)
  if(!is.null(mutual)){
    k = which(mutual)
    if(length(k) >= 2) keep = k     # without mutual edges the family is left as it is
  }

  ord = keep[order(len[keep], decreasing = TRUE)]

  return(ord[1:min(length(ord), n.repr)])
}


#' Alignment and consensus of one SV family
#'
#' Brings the members to one strand, aligns the (nearly) full-length ones,
#' trims the alignment by the coverage and returns the consensus.
#'
#' @param seqs Named vector of the sequences of the members of one family.
#' @param nestedness Table of the similarity search with the `strand` column,
#' used by \code{\link{getComponentSequences}} to orient the members; without it
#' the members are oriented by \code{\link{svFamilyOrient}}.
#' @param n.repr Maximum number of members to align.
#' @param mutual Named logical vector, see \code{\link{svFamilyRepr}}.
#' @param trim.cov,trim.wnd,trim.agree Parameters of \code{\link{svAlnTrim}}.
#' @param path.work,mafft.args,mafft.bin Parameters of \code{\link{svAlnMafft}}.
#'
#' @return A list with the trimmed alignment `mx`, the consensus `cons`,
#' the names `sv` of the aligned members, their strands and the coverage profile.
#'
#' @author Anna A. Igolkina
#' @export
svFamilyAln <- function(seqs, nestedness = NULL, n.repr = 20, mutual = NULL,
                        trim.cov = 3, trim.wnd = 5, trim.agree = 0.9,
                        path.work = tempdir(), mafft.args = '--op 5 --ep 0.2 --quiet', mafft.bin = 'mafft'){

  if(length(seqs) < 2) return(NULL)

  # ---- one strand ----
  # The strand of every pair is already known from the similarity search, so the
  # graph of the component is the first choice; k-mers are the fallback. The whole
  # family is oriented and not only the members to be aligned: the longest members
  # are often connected to each other only through the shorter ones, and a subset
  # of them falls apart into several components.
  seqs.ori = NULL
  method.ori = 'graph'
  if(!is.null(nestedness)){
    seqs.ori = tryCatch(getComponentSequences(names(seqs), seqs, nestedness),
                        error = function(e) NULL)
  }
  if(is.null(seqs.ori)){
    seqs.ori = svFamilyOrient(seqs)$seqs
    method.ori = 'kmer'
  }
  strand = ifelse(toupper(seqs.ori) == toupper(seqs[names(seqs.ori)]), '+', '-')
  names(strand) = names(seqs.ori)

  # ---- the members to align ----
  idx = svFamilyRepr(nchar(seqs.ori), n.repr = n.repr, mutual = mutual[names(seqs.ori)])
  if(length(idx) < 2) return(NULL)
  seqs.ori = seqs.ori[idx]
  strand = strand[names(seqs.ori)]

  mx.full = svAlnMafft(seqs.ori, path.work = path.work, args = mafft.args, bin = mafft.bin)

  tr = svAlnTrim(mx.full, min.cov = trim.cov, wnd = trim.wnd, min.agree = trim.agree)
  if(is.null(tr$mx)) return(NULL)

  cons = svAlnCons(tr$mx)

  res = list(mx      = tr$mx,
             cons    = cons,
             sv      = names(seqs.ori),
             strand  = strand,
             ori.by  = method.ori,
             cov     = tr$cov,
             agree   = tr$agree,
             trim    = c(beg = tr$beg, end = tr$end),
             n.aln   = length(idx))

  return(res)
}


#' Support of a terminal repeat by the members of a family
#'
#' A copy of an element is often truncated and carries only one of the two arms
#' of the terminal repeat. Asking every copy to show the whole repeat throws such
#' copies away, although each of them is evidence for the family. And the other
#' way round: the termini of the consensus of a long ragged family may come from
#' a single member that reaches further than the others, and a repeat found there
#' says nothing about the family.
#'
#' So the arms are counted separately: how many members carry the 5'-arm, how many
#' carry the 3'-arm, and how many carry both. The family is accepted when one of
#' the arms is present in at least `min.arm` members and at least `min.both`
#' members carry the pair.
#'
#' @param seqs Named vector of the sequences of the members, on one strand.
#' @param cons Consensus of the family; the candidate arms are taken from it, and
#' if it carries no repeat, from the member with the strongest one.
#' @param w Size of the terminal window to look for an arm in.
#' @param min.id Minimum identity of an arm to the candidate.
#' @param min.cov.arm Minimum share of the candidate arm that has to be covered.
#' @param min.arm Minimum number of members carrying one of the arms.
#' @param min.both Minimum number of members carrying both arms.
#' @param tir.par Parameters of \code{\link{findTIR}}.
#'
#' @return A list with the arms, their source, the counts `n5`, `n3`, `n.both`
#' and the verdict `tir`; or `NULL` if there is no candidate repeat at all.
#'
#' @author Anna A. Igolkina
#' @export
svFamilyArms <- function(seqs, cons, w = 150, min.id = 0.8, min.cov.arm = 0.7,
                         min.arm = 3, min.both = 1, tir.par = list()){

  arms <- function(s, r) list(a5 = substr(s, r$beg, r$beg + r$len - 1),
                              a3 = substr(revCompl(s), r$off.3, r$off.3 + r$len - 1))

  cand = NULL
  r = do.call(findTIR, c(list(s = cons), tir.par))
  if(!is.null(r)){
    cand = arms(cons, r); cand$src = 'cons'; cand$score = r$score
  } else {
    for(nm in names(seqs)){
      rr = do.call(findTIR, c(list(s = seqs[[nm]]), tir.par))
      if(is.null(rr)) next
      if(is.null(cand) || (rr$score > cand$score)){
        cand = arms(seqs[[nm]], rr); cand$src = nm; cand$score = rr$score
      }
    }
  }
  if(is.null(cand)) return(NULL)

  hasArm <- function(x, arm, side){
    n = nchar(x)
    w.cur = min(w, n)
    win = if(side == 5) substr(x, 1, w.cur) else revCompl(substr(x, n - w.cur + 1, n))
    h = svAlnLocal(seq2nt(win), seq2nt(arm), mism = -4)
    !is.null(h) && (h$len >= min.cov.arm * nchar(arm)) && (h$ident >= min.id)
  }

  v5 = vapply(seqs, hasArm, logical(1), arm = cand$a5, side = 5)
  v3 = vapply(seqs, hasArm, logical(1), arm = cand$a3, side = 3)

  res = list(a5 = cand$a5, a3 = cand$a3, src = cand$src,
             n5 = sum(v5), n3 = sum(v3), n.both = sum(v5 & v3), n = length(seqs))
  res$tir = (max(res$n5, res$n3) >= min.arm) && (res$n.both >= min.both)

  return(res)
}


#' Are the two sequences a circular permutation of each other
#'
#' The boundaries of an SV do not have to fall at the boundaries of an element:
#' when the element sits in a cluster of repeats, the caller may cut it at
#' different points, and the copies come out rotated - the beginning of one is the
#' end of the other. Such a pair covers each other completely in the similarity
#' search (the hits simply come as two blocks), but no aligner can put them into
#' one collinear alignment, and a family of them falls apart into a mosaic.
#'
#' The shift between the shared k-mers is taken modulo the length: for a rotated
#' pair it is one and the same value, for an ordinary pair it is also one but
#' without the wrap, and for unrelated sequences the shifts are scattered.
#'
#' Taken modulo the length, the shift is circular: `d` and `nb - d` are the same
#' rotation seen from the two sides, so a shift of `nb - 4` is a shift of 4 and no
#' rotation at all. Two collinear copies that differ by a few nucleotides of indel
#' drift give exactly that, and because the drift makes `j - i` change sign along
#' the sequence, a good share of the k-mers ends up on either side of zero - which
#' `side` alone reads as a rotation point. The rotation point therefore also has to
#' be far from the termini: the circular shift has to reach `min.shift` of the length.
#'
#' @param a,b Two sequences (strings).
#' @param k Length of the k-mer.
#' @param tol Tolerance of the shift, in nucleotides.
#' @param min.frac Minimum share of the shared k-mers agreeing on the shift.
#' @param min.side Minimum share of the k-mers on the smaller side of the rotation
#' point; below it the pair is collinear rather than rotated.
#' @param min.shift Minimum circular shift, as a share of the length of `b`. It
#' keeps the indel drift of a collinear pair from being read as a rotation.
#'
#' @return A list with `perm`, the shift `shift`, its circular value `shift.circ`,
#' the share `frac` of the agreeing k-mers and the share `side` of the smaller
#' block; or `NULL` if there are no shared k-mers at all.
#'
#' @author Anna A. Igolkina
#' @export
svPermuted <- function(a, b, k = 12, tol = 10, min.frac = 0.5, min.side = 0.15,
                      min.shift = 0.05){

  a = toupper(a); b = toupper(b)
  na = nchar(a); nb = nchar(b)
  if((na < 2 * k) || (nb < 2 * k)) return(NULL)

  ka = substring(a, 1:(na - k + 1), k:na)
  kb = substring(b, 1:(nb - k + 1), k:nb)
  pos.b = split(seq_along(kb), kb)
  pos.b = pos.b[lengths(pos.b) == 1]              # unique k-mers only: no tandem noise
  hit = pos.b[ka]
  keep = lengths(hit) == 1
  if(sum(keep) < 10) return(NULL)

  i = which(keep)
  j = unlist(hit[keep], use.names = FALSE)
  d = (j - i) %% nb                               # the shift modulo the length

  # the dominant shift, with the tolerance
  br = round(d / tol)
  tb = table(br)
  d.best = as.numeric(names(tb)[which.max(tb)]) * tol
  ok = abs(((d - d.best + nb/2) %% nb) - nb/2) <= tol

  # where the rotation point is: the shift wraps around the end of b
  wrap = (j[ok] < i[ok])
  side = min(mean(wrap), 1 - mean(wrap))

  # the shift is circular, so a shift of nb - 4 is a shift of 4
  shift.circ = min(d.best, nb - d.best)

  res = list(shift = d.best, shift.circ = shift.circ, frac = mean(ok), side = side,
             perm = (mean(ok) >= min.frac) && (side >= min.side) &&
                    (shift.circ >= min.shift * nb))

  return(res)
}


#' Share of rotated pairs among the members of a family
#'
#' @param seqs Named vector of the sequences of the members.
#' @param n.max Maximum number of members to compare pairwise.
#' @param ... Parameters of \code{\link{svPermuted}}.
#'
#' @return A list with the number of compared pairs, the number of rotated ones
#' and their share.
#'
#' @author Anna A. Igolkina
#' @export
svFamilyPermuted <- function(seqs, n.max = 12, ...){

  s = seqs[1:min(length(seqs), n.max)]
  n = length(s)
  if(n < 2) return(list(n.pairs = 0, n.perm = 0, frac = 0))

  n.pairs = 0; n.perm = 0
  for(i in 1:(n - 1)) for(j in (i + 1):n){
    r = svPermuted(s[[i]], s[[j]], ...)
    if(is.null(r)) next
    n.pairs = n.pairs + 1
    if(r$perm) n.perm = n.perm + 1
  }

  frac = if(n.pairs > 0) n.perm / n.pairs else 0
  list(n.pairs = n.pairs, n.perm = n.perm, frac = frac, n.cmp = n)
}


#' Structural order of a family, with the rotated families set apart
#'
#' Before looking for the termini it is worth asking whether the family has
#' termini at all. When the members are rotated relative to each other (see
#' \code{\link{svPermuted}}), no consensus of them means anything: the alignment
#' falls apart into a mosaic, the consensus comes out longer than any member, and
#' its ends belong to no one. Such families are given their own class instead of
#' being classified by the termini of an artefact.
#'
#' @param seqs Named vector of the sequences of the members, on one strand.
#' @param cons Consensus of the family.
#' @param max.perm Share of rotated pairs above which the family is called
#' `Rearranged` -- but only when the consensus has shown nothing more specific,
#' see the hierarchy LTR > TIR > Helitron > Rearranged > LINE_like. A rotated
#' family always gets `rotated members` in its `note`, whatever the class.
#' @param perm.par Parameters of \code{\link{svFamilyPermuted}}.
#' @param ... Parameters of \code{\link{svTeOrder}}.
#'
#' @return A one-row data frame of \code{\link{svTeOrder}} with the columns
#' `perm.frac` and `perm.pairs` added; the class is `Rearranged` when the members
#' are rotated.
#'
#' @author Anna A. Igolkina
#' @export
svFamilyOrder <- function(seqs, cons, max.perm = 0.2, perm.par = list(), ...){

  pm = do.call(svFamilyPermuted, c(list(seqs = seqs), perm.par))

  res = svTeOrder(cons, name = 'cons', ...)
  res$perm.frac = round(pm$frac, 3)
  res$perm.pairs = pm$n.pairs

  # Rearranged has its place in the same hierarchy as the rest:
  # LTR > TIR > Helitron > Rearranged > LINE_like. A family whose consensus carries
  # a terminal repeat or a Helitron hairpin is named by it even when some of its
  # members are rotated -- the signature is the specific evidence and the rotation
  # is not. Below them the rotation is worth more than a poly-A tail, so it takes
  # over LINE_like and unknown.
  if((pm$frac > max.perm) && (res$te.class %in% c('LINE_like', 'unknown'))){
    res$te.class = 'Rearranged'
  }
  if(pm$frac > max.perm){
    res$note = paste(c(res$note[nzchar(res$note)], 'rotated members'), collapse = ',')
  }

  return(res)
}
