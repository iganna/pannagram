# Structural order of mobile elements hidden in SVs:
#  - LTR   : terminal direct repeats
#  - TIR   : terminal inverted repeats
#  - LINE  : one-sided poly-A tail (poly-T on the opposite strand)
#
# Everything is done on the SV sequence itself, in base R:
# k-mer seeds -> the best band of diagonals -> chaining of ungapped segments.
# No external aligners and no Biostrings::pairwiseAlignment are used.


# ***********************************************************************
# ---- Low-level helpers ----

#' Exact k-mer seeds between two nucleotide vectors
#'
#' @param a Character vector of nucleotides (the first sequence).
#' @param b Character vector of nucleotides (the second sequence).
#' @param k Length of the seed.
#' @param max.rep Maximum number of occurrences of a k-mer in `b`;
#' more frequent k-mers are masked out as low-complexity.
#'
#' @return A data frame with columns `i` (position in `a`), `j` (position in `b`)
#' and `d` (diagonal, `j - i`), or `NULL` if there are no seeds.
#'
#' @author Anna A. Igolkina
#' @export
svTeSeeds <- function(a, b, k = 8, max.rep = 20){

  n.a = length(a)
  n.b = length(b)
  if((n.a < k) || (n.b < k)) return(NULL)

  km.a = substring(paste0(a, collapse = ''), 1:(n.a - k + 1), k:n.a)
  km.b = substring(paste0(b, collapse = ''), 1:(n.b - k + 1), k:n.b)

  # Ambiguous nucleotides do not produce seeds
  km.a[grepl('[^ACGT]', km.a)] = NA
  km.b[grepl('[^ACGT]', km.b)] = NA

  pos.b = split(seq_along(km.b), km.b)              # NAs are dropped by split
  pos.b = pos.b[lengths(pos.b) <= max.rep]          # mask low-complexity k-mers
  if(length(pos.b) == 0) return(NULL)

  hits = pos.b[km.a]                                # missing k-mers give NULL
  n.hits = lengths(hits)
  if(sum(n.hits) == 0) return(NULL)

  seeds = data.frame(i = rep(seq_along(km.a), n.hits),
                     j = unlist(hits, use.names = FALSE))
  seeds$d = seeds$j - seeds$i

  return(seeds)
}


#' Chain k-mer seeds into one terminal hit
#'
#' Finds the densest band of diagonals, merges seeds into ungapped segments
#' and chains them, allowing small shifts of the diagonal (indels).
#'
#' @param a Character vector of nucleotides (the first sequence).
#' @param b Character vector of nucleotides (the second sequence).
#' @param k Length of the seed.
#' @param band Half-width of the band of diagonals.
#' @param max.gap Maximum unaligned stretch between two chained segments.
#' @param max.rep Maximum number of occurrences of a k-mer, see \code{\link{svTeSeeds}}.
#' @param anchor Pair `c(off.a, off.b)`: only the diagonals having a seed within
#' `off.a` of the beginning of `a` and a seed within `off.b` of the end of `b`
#' are considered. Terminal repeats are what is looked for, and a strong internal
#' repeat should not outvote them.
#' @param d.range Range `c(min, max)` of the diagonals to look at, or `NULL` for all of them.
#' When the position of the hit is known in advance (an inverted repeat is expected right
#' at the termini, so its diagonal is close to zero), this keeps the terminal hit from
#' being outvoted by a stronger internal repeat somewhere else in the window.
#'
#' @return A list with `beg.a`, `end.a`, `beg.b`, `end.b`, `n.match`, `aln.len`,
#' `ident` and `n.seg`, or `NULL` if nothing was chained.
#'
#' @author Anna A. Igolkina
#' @export
svTeChain <- function(a, b, k = 8, band = 25, max.gap = 200, max.rep = 20,
                      d.range = NULL, anchor = NULL){

  seeds = svTeSeeds(a, b, k = k, max.rep = max.rep)
  if(is.null(seeds)) return(NULL)

  if(!is.null(d.range)){
    seeds = seeds[(seeds$d >= d.range[1]) & (seeds$d <= d.range[2]), , drop = FALSE]
    if(nrow(seeds) == 0) return(NULL)
  }

  # Only the diagonals that touch both termini are of interest. Without this the
  # densest band wins, and a strong internal repeat hides the terminal one -
  # which is exactly what happens in the long elements.
  if(!is.null(anchor)){
    d.ok = unique(seeds$d[seeds$i <= anchor[1]])
    d.ok = intersect(d.ok, unique(seeds$d[seeds$j >= (length(b) - anchor[2] + 1)]))
    if(length(d.ok) == 0) return(NULL)
    seeds = seeds[seeds$d %in% d.ok, , drop = FALSE]
  }

  n.b = length(b)

  # ---- The densest band of diagonals ----
  d.min = min(seeds$d)
  cnt = tabulate(seeds$d - d.min + 1)
  n.d = length(cnt)
  cs = cumsum(c(0, cnt))
  idx.d = seq_len(n.d)
  cnt.band = cs[pmin(n.d, idx.d + band) + 1] - cs[pmax(1, idx.d - band)]
  d.best = which.max(cnt.band) + d.min - 1

  seeds = seeds[abs(seeds$d - d.best) <= band, , drop = FALSE]
  if(nrow(seeds) == 0) return(NULL)

  # ---- Merge seeds of every diagonal into ungapped segments ----
  segs = c()
  for(d in sort(unique(seeds$d))){
    beg = sort(seeds$i[seeds$d == d])
    end = beg + k - 1
    br = c(TRUE, beg[-1] > (end[-length(end)] + 1))
    g = cumsum(br)
    segs = rbind(segs, data.frame(d   = d,
                                  beg = as.numeric(tapply(beg, g, min)),
                                  end = as.numeric(tapply(end, g, max))))
  }
  segs$len = segs$end - segs$beg + 1
  segs = segs[order(segs$beg, segs$end), , drop = FALSE]
  rownames(segs) = NULL

  # ---- Chaining of segments ----
  n.s = nrow(segs)
  score = segs$len
  prev = rep(0, n.s)
  if(n.s > 1){
    for(t in 2:n.s){
      for(u in 1:(t - 1)){
        # both coordinates should increase
        if(segs$beg[u] >= segs$beg[t]) next
        if((segs$beg[u] + segs$d[u]) >= (segs$beg[t] + segs$d[t])) next

        gap.a = segs$beg[t] - segs$end[u] - 1
        gap.b = (segs$beg[t] + segs$d[t]) - (segs$end[u] + segs$d[u]) - 1
        gap = max(gap.a, gap.b)
        if(gap > max.gap) next

        sc = score[u] + segs$len[t] - max(0, gap)
        if(sc > score[t]){
          score[t] = sc
          prev[t] = u
        }
      }
    }
  }

  # ---- Restore the best chain ----
  t = which.max(score)
  chain = c()
  while(t > 0){
    chain = c(t, chain)
    t = prev[t]
  }
  segs = segs[chain, , drop = FALSE]

  # Trim overlaps between consecutive segments of different diagonals
  if(nrow(segs) > 1){
    for(t in 2:nrow(segs)){
      shift = max(0,
                  segs$end[t-1] + 1 - segs$beg[t],
                  (segs$end[t-1] + segs$d[t-1] + 1) - (segs$beg[t] + segs$d[t]))
      segs$beg[t] = segs$beg[t] + shift
    }
    segs = segs[segs$beg <= segs$end, , drop = FALSE]
  }
  if(nrow(segs) == 0) return(NULL)

  # ---- Count matches along the chain ----
  n.match = 0
  aln.len = 0
  for(t in 1:nrow(segs)){
    i1 = segs$beg[t]; i2 = segs$end[t]; d = segs$d[t]
    j1 = i1 + d;      j2 = i2 + d
    if(j1 < 1)   { i1 = i1 + (1 - j1);   j1 = 1 }
    if(j2 > n.b) { i2 = i2 - (j2 - n.b); j2 = n.b }
    if(i2 < i1) next
    nt.a = a[i1:i2]
    n.match = n.match + sum((nt.a == b[j1:j2]) & (nt.a %in% c('A','C','G','T')))
    aln.len = aln.len + (i2 - i1 + 1)
  }
  if(nrow(segs) > 1){
    n.seg = nrow(segs)
    gap.a = segs$beg[-1] - segs$end[-n.seg] - 1
    gap.b = (segs$beg[-1] + segs$d[-1]) - (segs$end[-n.seg] + segs$d[-n.seg]) - 1
    aln.len = aln.len + sum(pmax(gap.a, gap.b, 0))
  }
  if(aln.len == 0) return(NULL)

  res = list(beg.a   = segs$beg[1],
             end.a   = segs$end[nrow(segs)],
             beg.b   = segs$beg[1] + segs$d[1],
             end.b   = segs$end[nrow(segs)] + segs$d[nrow(segs)],
             n.match = n.match,
             aln.len = aln.len,
             ident   = n.match / aln.len,
             n.seg   = nrow(segs))

  return(res)
}


#' Ungapped extension of a homopolymer run from the beginning of a vector
#'
#' @param v Character vector of nucleotides, the run is searched from `v[1]`.
#' @param nt Nucleotide of the run, `A` or `T`.
#' @param mism.penalty Penalty for a non-matching position.
#'
#' @return A list with `len` and `ident` of the run (`len = 0` if there is no run).
#'
#' @author Anna A. Igolkina
#' @export
svTeRun <- function(v, nt = 'A', mism.penalty = 2){

  if(length(v) == 0) return(list(len = 0, ident = 0))

  m = (v == nt)
  sc = cumsum(ifelse(m, 1, -mism.penalty))
  k.best = which.max(sc)
  if(sc[k.best] <= 0) return(list(len = 0, ident = 0))

  return(list(len = k.best, ident = mean(m[1:k.best])))
}


#' Check that a sequence is not low-complexity
#'
#' Homopolymer and AT-rich termini easily produce formally perfect terminal
#' repeats, which have nothing to do with mobile elements.
#'
#' @param s A single sequence (string).
#' @param max.nt.frac Maximum fraction of one nucleotide.
#' @param min.nt Minimum number of distinct nucleotides.
#'
#' @return `TRUE` if the sequence is complex enough.
#'
#' @author Anna A. Igolkina
#' @export
svTeComplex <- function(s, max.nt.frac = 0.7, min.nt = 3, min.k2 = 0){

  v = seq2nt(toupper(s))
  v = v[v %in% c('A','C','G','T')]
  n = length(v)
  if(n == 0) return(FALSE)

  tbl = table(v)
  if(length(tbl) < min.nt) return(FALSE)
  if((max(tbl) / n) > max.nt.frac) return(FALSE)

  # The share of the commonest nucleotide alone throws away real signatures:
  # GGGGTGGGC is 78% G and is still a terminal repeat, while ATATATAT is only
  # 50% A and is nothing. The number of distinct dinucleotides separates them.
  if(min.k2 > 0){
    if(n < 2) return(FALSE)
    if(length(unique(paste0(v[-n], v[-1]))) < min.k2) return(FALSE)
  }

  return(TRUE)
}


# ***********************************************************************
# ---- Detectors ----

#' Optimal local alignment of two short sequences
#'
#' A plain Smith-Waterman with linear gaps (Smith & Waterman, 1981), the same
#' way the terminal repeats are looked for by `einverted` (EMBOSS) and by the
#' structural part of LTRharvest (Ellinghaus et al., 2008): first the optimal
#' repeat is found, and only then the constraints on its position and length are
#' applied. Seed-and-extend does the opposite and loses the degenerate repeats,
#' where no exact seed of the required length survives.
#'
#' Only the terminal windows are aligned, so the quadratic cost stays small.
#'
#' @param a,b Character vectors of nucleotides.
#' @param match,mism,gap Scores of a match, of a mismatch and of a gap.
#'
#' @return A list with the coordinates of the best local alignment in both
#' sequences, its length, the number of matches, the identity and the score;
#' or `NULL` if there is no positive-scoring alignment.
#'
#' @author Anna A. Igolkina
#' @export
svAlnLocal <- function(a, b, match = 2, mism = -3, gap = -5){

  n = length(a)
  m = length(b)
  if((n == 0) || (m == 0)) return(NULL)

  H = matrix(0, n + 1, m + 1)     # scores
  P = matrix(0L, n + 1, m + 1)    # 1 diagonal, 2 up, 3 left

  for(i in 1:n){
    s.row = ifelse(a[i] == b, match, mism)
    s.row[!(a[i] %in% c('A','C','G','T'))] = mism
    for(j in 1:m){
      v = c(H[i, j] + s.row[j], H[i, j+1] + gap, H[i+1, j] + gap, 0)
      k = which.max(v)
      H[i+1, j+1] = v[k]
      P[i+1, j+1] = k
    }
  }

  if(max(H) <= 0) return(NULL)
  ij = which(H == max(H), arr.ind = TRUE)[1, ]
  i = ij[[1]]; j = ij[[2]]
  end.a = i - 1; end.b = j - 1
  n.match = 0; aln.len = 0

  while(H[i, j] > 0){
    k = P[i, j]
    if(k == 1){
      if(a[i-1] == b[j-1]) n.match = n.match + 1
      i = i - 1; j = j - 1
    } else if(k == 2){
      i = i - 1
    } else {
      j = j - 1
    }
    aln.len = aln.len + 1
  }

  res = list(beg.a = i, end.a = end.a, beg.b = j, end.b = end.b,
             len = aln.len, n.match = n.match, ident = n.match / aln.len,
             score = max(H))

  return(res)
}


#' Find terminal direct repeats (LTR) in a sequence
#'
#' The first and the last window of the sequence are compared on the same strand.
#' The hit is accepted only if it is anchored at both termini of the sequence.
#'
#' @param s A single sequence (string).
#' @param w.max,w.min Sizes of the nested terminal windows.
#' @param min.len Minimum length of the LTR.
#' @param min.ident Minimum identity between the two LTR copies.
#' @param max.offset Maximum distance between the hit and the terminus of the sequence.
#' @param k,band,max.gap Parameters of \code{\link{svTeChain}}.
#'
#' @return A list describing the LTR pair, or `NULL`.
#'
#' @author Anna A. Igolkina
#' @export
findLTR <- function(s, w.max = 2000, min.len = 40, min.ident = 0.7,
                    max.offset = 300, off.free = 20, ident.per.off = 0.0005,
                    max.nt.frac = 0.7, k = 8, band = 25, max.gap = 200){

  if(length(s) != 1) stop('There should be only one sequence')
  s = toupper(s)
  len = nchar(s)
  if(len < 2 * min.len) return(NULL)

  w = min(w.max, floor(len / 2))
  a = seq2nt(substr(s, 1, w))
  b = seq2nt(substr(s, len - w + 1, len))

  h = svTeChain(a, b, k = k, band = band, max.gap = max.gap,
                anchor = c(max.offset, max.offset))
  if(is.null(h)) return(NULL)

  # Anchoring: the left copy starts at the beginning, the right one ends at the end.
  # A copy that does not touch the terminus is a weaker claim, so it has to be
  # a better repeat - the required identity grows with the offset.
  off = max(h$beg.a - 1, w - h$end.b)
  if(off > max.offset) return(NULL)

  ltr.len = min(h$end.a - h$beg.a + 1, h$end.b - h$beg.b + 1)
  if(ltr.len < min.len) return(NULL)
  if(h$ident < (min.ident + ident.per.off * max(0, off - off.free))) return(NULL)

  # Low-complexity termini are not LTRs
  if(!svTeComplex(substr(s, h$beg.a, h$end.a), max.nt.frac = max.nt.frac)) return(NULL)

  # The two copies should not overlap in the sequence
  beg.right = (len - w) + h$beg.b
  end.right = (len - w) + h$end.b
  if(h$end.a >= beg.right) return(NULL)

  # The same measure as in findTIR, so that the two kinds of evidence can be
  # compared at all: how unlikely is a repeat this good, this close to the ends.
  off5 = h$beg.a - 1
  off3 = w - h$end.b
  sc = 2 * round(ltr.len * h$ident) - ltr.len
  e = 0.33 * (off5 + 1) * (off3 + 1) * exp(-log(3) * sc)

  res = list(len       = ltr.len,
             ident     = h$ident,
             score     = sc,
             e         = e,
             beg.left  = h$beg.a,
             end.left  = h$end.a,
             beg.right = beg.right,
             end.right = end.right,
             off5      = off5,
             off3      = off3,
             off       = off,
             tg.ca     = (substr(s, 1, 2) == 'TG') && (substr(s, len - 1, len) == 'CA'))

  return(res)
}


#' Find terminal inverted repeats (TIR) in a sequence
#'
#' The beginning of the sequence is aligned with the reverse complement of its
#' end by \code{\link{svAlnLocal}}, and the found repeat is then checked against
#' the constraints. This order matters: the boundaries of an SV do not have to
#' coincide with the boundaries of the element, so the arms can start tens of
#' nucleotides away from the termini, and a repeat that is looked for only at
#' the exact termini is lost.
#'
#' @param s A single sequence (string).
#' @param w.max,w.min Sizes of the nested terminal windows.
#' @param min.len Minimum length of the TIR.
#' @param min.ident Minimum identity between the two copies.
#' @param min.score Minimum score of the local alignment; together with the
#' identity it keeps random hits out, see the calibration in the tests.
#' @param max.offset Maximum distance between an arm and the terminus.
#' @param off.free Distance from the terminus that costs nothing; further away
#' the required score grows by `score.per.off` for every nucleotide. An arm that
#' does not touch the terminus is a weaker claim and has to be a better repeat.
#' @param score.per.off Growth of the required score per nucleotide of the offset.
#' @param max.nt.frac Maximum share of one nucleotide, against low-complexity termini.
#'
#' @return A list describing the TIR pair, or `NULL`.
#'
#' @author Anna A. Igolkina
#' @export
findTIR <- function(s, w.max = 300, w.min = 60, mism = c(-4, -3), min.len = 5,
                    max.offset = 100, alpha = 0.01, n.cand = 4,
                    max.nt.frac = 0.85, min.k2 = 3){

  if(length(s) != 1) stop('There should be only one sequence')
  s = toupper(s)
  len = nchar(s)
  if(len < 2 * min.len) return(NULL)

  # What makes a terminal repeat believable is not its length but the fact that
  # it sits at the termini. A perfect 5-mer at offset zero and a 9-mer at offset
  # ten are both far less likely by chance than the same match somewhere within
  # a hundred nucleotides of the end, and a fixed length cutoff cannot say that.
  # So a candidate is kept when the expected number of chance hits at least as
  # good and at least as close to the ends is below alpha:
  #
  #   E = K * (off5 + 1) * (off3 + 1) * exp(-lambda * S),   S = matches - mismatches
  #
  # lambda = log(3) is exact for the +1/-1 score over four equiprobable letters,
  # K is the usual ungapped-alignment constant. This is the Karlin-Altschul form
  # used by BLAST, with the search space taken as the offsets actually reached
  # rather than the whole window (Karlin & Altschul 1990, PNAS 87:2264).
  lambda = log(3)
  k.const = 0.33

  # Nested windows: in a narrow one a strong internal repeat is simply absent
  # and cannot outscore the terminal signature.
  w.all = unique(pmin(c(w.min, 150, w.max), floor(len / 2)))
  w.all = w.all[w.all >= min.len]

  # Every window, every mismatch penalty and every declumped candidate is one
  # more look at the same sequence, and without a correction a random sequence
  # yields a hit about once per run. The budget alpha is for the whole search.
  n.tests = length(mism) * length(w.all) * n.cand

  best = NULL
  for(ms in mism) for(w in w.all){
    a = seq2nt(substr(s, 1, w))
    b = seq2nt(revCompl(substr(s, len - w + 1, len)))

    # The optimal local alignment is not necessarily the terminal one, and
    # discarding the window because of that loses the repeat. The found segment
    # is masked out and the search repeated, so the next candidates are seen too
    # (declumping, as in Waterman & Eggert 1987, J Mol Biol 197:723).
    for(it in seq_len(n.cand)){
      cur = svAlnLocal(a, b, mism = ms)
      if(is.null(cur)) break
      a[cur$beg.a:cur$end.a] = 'N'
      b[cur$beg.b:cur$end.b] = 'N'

      off5 = cur$beg.a - 1
      off3 = cur$beg.b - 1
      if(max(off5, off3) > max.offset) next
      if(cur$len < min.len) next

      sc = 2 * cur$n.match - cur$len            # +1 per match, -1 per mismatch or gap
      if(sc < min.len) next
      e = k.const * (off5 + 1) * (off3 + 1) * exp(-lambda * sc) * n.tests
      if(e > alpha) next
      if(!svTeComplex(substr(s, cur$beg.a, cur$end.a),
                      max.nt.frac = max.nt.frac, min.k2 = min.k2)) next

      cur$off5 = off5; cur$off3 = off3; cur$sc = sc; cur$e = e
      if(is.null(best) || (e < best$e)) best = cur
    }
  }
  if(is.null(best)) return(NULL)

  res = list(len    = best$len,
             ident  = best$ident,
             score  = best$sc,
             e      = best$e,
             beg    = best$beg.a,
             off.3  = best$beg.b,
             off5   = best$off5,
             off3   = best$off3,
             off    = max(best$off5, best$off3),
             method = 'local',
             motif  = substr(s, best$beg.a, min(best$beg.a + 5, len)),
             arm    = substr(s, best$beg.a, best$end.a))

  return(res)
}


#' Find a hairpin (a palindrome with a loop) near a terminus of a sequence
#'
#' The search is vectorised over anti-diagonals: for a fixed sum of the positions
#' of the two arms, the complementarity of all pairs is checked at once, and the
#' length of the stem is the run of complementary pairs adjacent to the loop.
#'
#' @param s A single sequence (string).
#' @param w Size of the terminal window to look at.
#' @param side Terminus to look at, `right` (3') or `left` (5').
#' @param min.stem Minimum length of the stem.
#' @param loop.min,loop.max Range of the length of the loop.
#'
#' @return A list with `stem`, `loop`, `beg`, `end` and `dist` to the terminus, or `NULL`.
#'
#' @author Anna A. Igolkina
#' @export
findHairpin <- function(s, w = 80, side = 'right', min.stem = 7, loop.min = 3, loop.max = 10){

  if(length(s) != 1) stop('There should be only one sequence')
  s = toupper(s)
  len = nchar(s)
  w = min(w, len)
  if(w < 2 * min.stem + loop.min) return(NULL)

  shift = if(side == 'right') len - w else 0
  v = seq2nt(substr(s, shift + 1, shift + w))
  compl = c(A = 'T', T = 'A', G = 'C', C = 'G')
  v.c = compl[v]

  best = NULL
  for(d in (loop.min + 3):(2 * w - loop.min - 1)){

    j.lo = max(1, d - w)
    j.hi = floor((d - loop.min - 1) / 2)
    if((j.hi - j.lo + 1) < min.stem) next
    j = j.lo:j.hi

    b = (v.c[j] == v[d - j])
    b[is.na(b)] = FALSE
    z = ifelse(b, 0L, seq_along(b))
    r = seq_along(b) - cummax(z)

    for(loop in loop.min:loop.max){
      if(((d - loop - 1) %% 2) != 0) next
      j0 = (d - loop - 1) / 2
      i0 = j0 - j.lo + 1
      if((i0 < 1) || (i0 > length(j))) next
      stem = r[i0]
      if(stem < min.stem) next
      if(!is.null(best) && (stem <= best$stem)) next
      best = list(stem = unname(stem), loop = loop,
                  beg = unname(shift + j0 - stem + 1),
                  end = unname(shift + (d - (j0 - stem + 1))))
    }
  }

  if(is.null(best)) return(NULL)
  best$dist = if(side == 'right') len - best$end else best$beg - 1

  return(best)
}


#' Find the 3'-terminal signature of a Helitron
#'
#' Helitrons roll out by the rolling-circle mechanism and have neither terminal
#' repeats nor a poly-A tail. Their conserved end is the 3' one: a hairpin
#' followed by the `CTRR` terminus (Kapitonov & Jurka, 2001; the same pair of
#' features is what HelitronScanner scores). The `TC` at the 5'-terminus is only
#' reported, because 5'-truncated copies are common. Both strands are checked.
#'
#' @param s A single sequence (string).
#' @param w Size of the terminal window to look for the hairpin in.
#' @param min.stem Minimum length of the stem of the hairpin.
#' @param loop.min,loop.max Range of the length of the loop.
#' @param max.dist Maximum distance between the hairpin and the 3'-terminus.
#' @param ctrr Allowed 3'-termini.
#' @param max.offset Maximum shift of the termini.
#'
#' @return A list with `strand`, `stem`, `loop`, `dist`, `ctrr`, `ctrr.offset` and `tc`,
#' or `NULL`.
#'
#' @author Anna A. Igolkina
#' @export
findHelitron <- function(s, w = 80, min.stem = 7, loop.min = 3, loop.max = 10,
                         max.dist = 80, ctrr = c('CTAA', 'CTAG', 'CTGA', 'CTGG'),
                         max.offset = 5){

  if(length(s) != 1) stop('There should be only one sequence')
  s = toupper(s)
  if(nchar(s) < 2 * min.stem + loop.min + 4) return(NULL)

  res = NULL
  for(strand in c('+', '-')){

    x = if(strand == '+') s else revCompl(s)
    len = nchar(x)

    i.ctrr = NA
    for(off in 0:max.offset){
      m = substr(x, len - 3 - off, len - off)
      if(m %in% ctrr){ i.ctrr = off; s.ctrr = m; break }
    }
    if(is.na(i.ctrr)) next

    hp = findHairpin(x, w = w, side = 'right', min.stem = min.stem,
                     loop.min = loop.min, loop.max = loop.max)
    if(is.null(hp)) next
    if(hp$dist > max.dist) next

    tc = FALSE
    for(off in 0:max.offset) if(substr(x, 1 + off, 2 + off) == 'TC'){ tc = TRUE; break }

    cur = list(strand = strand, stem = hp$stem, loop = hp$loop, dist = hp$dist,
               ctrr = s.ctrr, ctrr.offset = i.ctrr, tc = tc)
    if(is.null(res) || (cur$stem + cur$tc) > (res$stem + res$tc)) res = cur
  }

  return(res)
}


#' Find a known short terminal signature
#'
#' Some superfamilies of DNA transposons have very short but strictly conserved
#' termini, which are too short for \code{\link{findTIR}}: the classical example
#' is `CACTA`, where the element starts with `CACTA` and ends with `TAGTG`.
#' The requirement to match a fixed motif at *both* termini makes the test
#' specific: a random sequence gives such a pair with probability ~4^(-2L).
#'
#' The check does not depend on the orientation of the SV, because the pair
#' (motif, reverse complement of the motif) is the same on both strands.
#'
#' @param s A single sequence (string).
#' @param motifs Vector of motifs to look for at the 5'-terminus; at the 3'-terminus
#' the reverse complement of the motif is expected.
#' @param max.offset Maximum shift of the motif from the terminus,
#' the boundaries of an SV can be inaccurate by a couple of nucleotides.
#' @param max.mism Maximum total number of mismatches in both termini.
#'
#' @return A list with `motif`, `offset` and `mism`, or `NULL`.
#'
#' @author Anna A. Igolkina
#' @export
findTermMotif <- function(s, motifs = c('CACTA', 'CACTG'), max.offset = 2, max.mism = 0){

  if(length(s) != 1) stop('There should be only one sequence')
  s = toupper(s)
  len = nchar(s)

  for(m in toupper(motifs)){
    n.m = nchar(m)
    m.nt = seq2nt(m)
    m.rc = seq2nt(revCompl(m))

    for(off in 0:max.offset){
      if(len < 2 * (n.m + off)) next

      n.mism = sum(seq2nt(substr(s, 1 + off, off + n.m)) != m.nt) +
               sum(seq2nt(substr(s, len - off - n.m + 1, len - off)) != m.rc)

      if(n.mism <= max.mism) return(list(motif = m, offset = off, mism = n.mism))
    }
  }

  return(NULL)
}


#' Find a one-sided poly-A tail in a sequence
#'
#' A poly-A tail at the right terminus and a poly-T head at the left one
#' correspond to the two orientations of a LINE-like element.
#'
#' @param s A single sequence (string).
#' @param w Size of the terminal window to look at.
#' @param min.len Minimum length of the run.
#' @param min.ident Minimum fraction of A (T) inside the run.
#' @param mism.penalty Penalty for a non-matching position.
#' @param w.tail,min.frac A degraded tail is accepted if the share of A in the
#' terminal window of `w.tail` reaches `min.frac` on one side only.
#'
#' @return A list with `len`, `ident`, `side` (`right`, `left`, `both`, `none`)
#' and the lengths of both runs.
#'
#' @author Anna A. Igolkina
#' @export
findPolyA <- function(s, w = 60, min.len = 7, min.ident = 0.8, mism.penalty = 2,
                      w.tail = 40, min.frac = 0.70){

  if(length(s) != 1) stop('There should be only one sequence')
  s = toupper(s)
  len = nchar(s)
  w = min(w, len)

  # poly-A at the 3'-end, read from the right terminus inside
  run.r = svTeRun(rev(seq2nt(substr(s, len - w + 1, len))), 'A', mism.penalty)
  # poly-T at the 5'-end, read from the left terminus inside
  run.l = svTeRun(seq2nt(substr(s, 1, w)), 'T', mism.penalty)

  ok.r = (run.r$len >= min.len) && (run.r$ident >= min.ident)
  ok.l = (run.l$len >= min.len) && (run.l$ident >= min.ident)

  # A tail of an old copy is broken by substitutions and no longer reads as a run,
  # but the terminus stays A-rich. The share of A in the terminal window catches
  # such tails; the requirement of one-sidedness below keeps it specific.
  if(!ok.r && !ok.l && (len >= 2 * w.tail)){
    f.r = mean(seq2nt(substr(s, len - w.tail + 1, len)) == 'A')
    f.l = mean(seq2nt(substr(s, 1, w.tail)) == 'T')
    if((f.r >= min.frac) && (f.l < min.frac)){
      run.r = list(len = round(f.r * w.tail), ident = f.r); ok.r = TRUE
    } else if((f.l >= min.frac) && (f.r < min.frac)){
      run.l = list(len = round(f.l * w.tail), ident = f.l); ok.l = TRUE
    }
  }

  if(ok.r && !ok.l){
    side = 'right'; run = run.r
  } else if(ok.l && !ok.r){
    side = 'left';  run = run.l
  } else if(ok.l && ok.r){
    side = 'both';  run = if(run.r$len >= run.l$len) run.r else run.l
  } else {
    side = 'none';  run = list(len = 0, ident = 0)
  }

  res = list(len       = run$len,
             ident     = run$ident,
             side      = side,
             len.right = run.r$len,
             len.left  = run.l$len)

  return(res)
}


# ***********************************************************************
# ---- Classification ----

#' Structural order of one SV sequence
#'
#' Combines \code{\link{findLTR}}, \code{\link{findTIR}}, \code{\link{findTermMotif}},
#' \code{\link{findHelitron}} and \code{\link{findPolyA}} into one class:
#' `LTR`, `TIR`, `Helitron`, `LINE_like` or `unknown`.
#' A known short terminal signature (`CACTA`) is reported as `TIR` with `tir.method = 'motif'`.
#'
#' @param s A single sequence (string).
#' @param name Name of the sequence, used in the output.
#' @param ltr.par,tir.par,polya.par,motif.par,hel.par Named lists with parameters of the detectors.
#'
#' @return A one-row data frame with the class and all supporting values.
#'
#' @author Anna A. Igolkina
#' @export
svTeOrder <- function(s, name = NULL, ltr.par = list(), tir.par = list(),
                      polya.par = list(), motif.par = list(), hel.par = list()){

  if(length(s) != 1) stop('There should be only one sequence')
  if(is.null(name)) name = if(!is.null(names(s))) names(s) else 'seq'
  s = toupper(unname(s))

  ltr = do.call(findLTR,   c(list(s = s), ltr.par))
  tir = do.call(findTIR,   c(list(s = s), tir.par))
  pa  = do.call(findPolyA, c(list(s = s), polya.par))
  mtf = do.call(findTermMotif, c(list(s = s), motif.par))
  hel = do.call(findHelitron,  c(list(s = s), hel.par))

  # CACTA gets no separate rule: findTIR judges a terminal repeat by how unlikely
  # it is at the position where it sits, and a five-letter signature at the very
  # ends passes that test on its own. The motif is kept as an annotation only.

  note = c()
  if(!is.null(ltr) && !is.null(tir)) note = c(note, 'ltr_and_tir')
  if(!is.null(mtf))                  note = c(note, paste0('motif_', mtf$motif))
  if(!is.null(hel) && hel$tc)        note = c(note, 'helitron_tc')
  if(!is.null(ltr) && ltr$tg.ca)     note = c(note, 'tg_ca')
  if(pa$side == 'both')              note = c(note, 'polyA_both_ends')

  # The detectors are not equally specific, so they are read in a fixed order and
  # not by the strength of the evidence. A terminal repeat, direct or inverted, is
  # a statement about both termini at once and is hard to come by at random; the
  # Helitron hairpin with a CTRR terminus is next; the poly-A tail is the weakest
  # of them all -- an A-rich stretch turns up at the end of anything in an AT-rich
  # genome -- so LINE_like is only what is left when nothing else has fired.
  # Both terminal repeats are still weighed on one scale -- minus the decimal
  # logarithm of the expected number of chance hits -- so that the reported
  # evidence stays comparable between them.
  ev = c(LTR       = if(is.null(ltr)) 0 else -log10(ltr$e) / 3,
         TIR       = if(is.null(tir)) 0 else -log10(tir$e) / 3,
         Helitron  = if(is.null(hel)) 0 else hel$stem / 8 + 0.5 * hel$tc,
         LINE_like = if(pa$side %in% c('left', 'right')) pa$len / 12 else 0)

  hierarchy = c('LTR', 'TIR', 'Helitron', 'LINE_like')
  fired = hierarchy[ev[hierarchy] > 0]
  te.class = if(length(fired) > 0) fired[1] else 'unknown'
  if((te.class != 'LINE_like') && (pa$side %in% c('left', 'right'))) note = c(note, 'polyA')
  if((te.class != 'Helitron') && !is.null(hel))                      note = c(note, 'helitron')
  if((te.class != 'TIR') && !is.null(tir))                           note = c(note, 'tir')
  if((te.class != 'LTR') && !is.null(ltr))                           note = c(note, 'ltr')

  res = data.frame(
    sv          = name,
    len         = nchar(s),
    te.class    = te.class,
    ltr.len     = if(is.null(ltr)) 0 else ltr$len,
    ltr.ident   = if(is.null(ltr)) 0 else round(ltr$ident, 4),
    ltr.off     = if(is.null(ltr)) NA else ltr$off,
    ltr.ev      = round(unname(ev[['LTR']]), 3),
    tir.len     = if(is.null(tir)) 0 else tir$len,
    tir.ident   = if(is.null(tir)) 0 else round(tir$ident, 4),
    tir.off     = if(is.null(tir)) NA else tir$off,
    tir.ev      = round(unname(ev[['TIR']]), 3),
    tir.motif   = if(is.null(tir)) '' else tir$motif,
    tir.method  = if(is.null(tir)) '' else tir$method,
    hel.stem    = if(is.null(hel)) 0  else hel$stem,
    hel.dist    = if(is.null(hel)) 0  else hel$dist,
    hel.ctrr    = if(is.null(hel)) '' else hel$ctrr,
    hel.ev      = round(unname(ev[['Helitron']]), 3),
    line.ev     = round(unname(ev[['LINE_like']]), 3),
    polyA.len   = pa$len,
    polyA.ident = round(pa$ident, 4),
    polyA.side  = pa$side,
    note        = paste0(note, collapse = ','),
    stringsAsFactors = FALSE)

  return(res)
}


#' Structural order of SV families
#'
#' Every family is described by its longest members (representatives),
#' because truncated copies simply have no termini to look at.
#' The class of the family is the most frequent class among the representatives.
#'
#' @param seqs Named vector of SV sequences.
#' @param families Named vector `SV -> family`, or a two-column data frame.
#' @param n.repr Maximum number of representatives per family.
#' @param len.frac Minimum length of a representative relative to the longest member.
#' @param min.support Minimum fraction of representatives supporting the class.
#' @param echo Show the progress.
#' @param ... Parameters passed to \code{\link{svTeOrder}}.
#'
#' @return A list with two data frames: `families` and `repr`.
#'
#' @author Anna A. Igolkina
#' @export
svTeOrderFamilies <- function(seqs, families, n.repr = 10, len.frac = 0.9,
                              min.support = 1/3, echo = TRUE, ...){

  if(is.data.frame(families)){
    fam = as.character(families[[2]])
    names(fam) = as.character(families[[1]])
  } else {
    fam = as.character(families)
    names(fam) = names(families)
  }
  if(is.null(names(fam))) stop('Families should be named by SV names')

  idx.absent = !(names(fam) %in% names(seqs))
  if(sum(idx.absent) > 0){
    pokazAttention('Sequences are not found for', sum(idx.absent), 'SVs, they are skipped')
    fam = fam[!idx.absent]
  }
  if(length(fam) == 0) stop('No SV sequences correspond to the families')

  fam.names = unique(fam)
  res.repr = c()
  res.fam = c()

  for(i.f in seq_along(fam.names)){
    f = fam.names[i.f]
    if(echo && (i.f %% 100 == 0)) pokaz('Family', i.f, 'of', length(fam.names))

    sv.names = names(fam)[fam == f]
    sv.len = nchar(seqs[sv.names])
    sv.names = sv.names[order(sv.len, decreasing = TRUE)]
    sv.len = sort(sv.len, decreasing = TRUE)

    # Representatives: the (nearly) full-length copies
    idx.repr = which(sv.len >= len.frac * sv.len[1])
    idx.repr = idx.repr[1:min(length(idx.repr), n.repr)]

    df = c()
    for(i.r in idx.repr){
      df = rbind(df, svTeOrder(seqs[[sv.names[i.r]]], name = sv.names[i.r], ...))
    }
    df = cbind(family = f, df, stringsAsFactors = FALSE)
    res.repr = rbind(res.repr, df)

    # ---- Vote ----
    cl = df$te.class[df$te.class != 'unknown']
    if(length(cl) == 0){
      te.class = 'unknown'
      support = 0
      i.best = 1
    } else {
      tbl = sort(table(cl), decreasing = TRUE)
      te.class = names(tbl)[1]
      support = tbl[[1]] / nrow(df)
      if(support < min.support){
        te.class = 'unknown'
        i.best = 1
      } else {
        i.best = which(df$te.class == te.class)[1]   # the longest one, they are sorted
      }
    }

    row.fam = data.frame(
      family      = f,
      n.members   = length(sv.names),
      n.repr      = nrow(df),
      len.max     = sv.len[1],
      te.class    = te.class,
      support     = round(support, 4),
      sv.repr     = df$sv[i.best],
      ltr.len     = df$ltr.len[i.best],
      ltr.ident   = df$ltr.ident[i.best],
      tir.len     = df$tir.len[i.best],
      tir.ident   = df$tir.ident[i.best],
      tir.motif   = df$tir.motif[i.best],
      tir.method  = df$tir.method[i.best],
      hel.stem    = df$hel.stem[i.best],
      hel.dist    = df$hel.dist[i.best],
      hel.ctrr    = df$hel.ctrr[i.best],
      polyA.len   = df$polyA.len[i.best],
      polyA.side  = df$polyA.side[i.best],
      note        = df$note[i.best],
      stringsAsFactors = FALSE)

    res.fam = rbind(res.fam, row.fam)
  }

  rownames(res.fam) = NULL
  rownames(res.repr) = NULL

  return(list(families = res.fam, repr = res.repr))
}
