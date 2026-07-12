#' Build a coordinate-compressed score matrix from alignment intervals
#'
#' Instead of a per-base matrix of width `chr.len` (up to gigabases), the
#' chromosome is cut into segments whose boundaries are the endpoints of all
#' intervals. Within a segment every query row has a constant score, so the
#' column-wise arg-max used by \code{findBestChromosome} is identical to the
#' dense per-base version, while memory drops from O(chr.len) to O(#intervals).
#'
#' @param intervals List (one element per query chromosome) of data.frames with
#'   columns \code{beg}, \code{end} (1-based, inclusive, beg<=end) and \code{score}.
#'   Intervals of the same row are applied in order (last one wins on overlap),
#'   matching the dense assignment \code{pos[i, beg:end] <- score}.
#' @param chr.len Integer, length of the chromosome (number of bases).
#' @return list(pos = k x S score matrix, seg.pos = real start of each segment,
#'   seg.len = length in bases of each segment).
buildSegScore <- function(intervals, chr.len) {
  k <- length(intervals)

  # Breakpoints: interval starts and ends+1, plus the chromosome bounds.
  bps <- c(1, chr.len + 1)
  for (iv in intervals) {
    if (!is.null(iv) && nrow(iv) > 0) bps <- c(bps, iv$beg, iv$end + 1)
  }
  bps <- sort(unique(bps))
  bps <- bps[bps >= 1 & bps <= chr.len + 1]
  if (bps[1] != 1) bps <- c(1, bps)
  if (bps[length(bps)] != chr.len + 1) bps <- c(bps, chr.len + 1)

  S <- length(bps) - 1
  seg.pos <- bps[-length(bps)]          # real start coordinate of each segment
  seg.len <- diff(bps)                  # length in bases of each segment

  pos <- matrix(0, nrow = k, ncol = S)
  for (i in seq_len(k)) {
    iv <- intervals[[i]]
    if (is.null(iv) || nrow(iv) == 0) next
    s.beg <- findInterval(iv$beg, bps)  # segment index containing beg
    s.end <- findInterval(iv$end, bps)  # segment index containing end
    for (r in seq_len(nrow(iv))) {
      pos[i, s.beg[r]:s.end[r]] <- iv$score[r]   # overwrite -> last wins
    }
  }

  list(pos = pos, seg.pos = seg.pos, seg.len = seg.len)
}


#' Select Best Correspondence Regions Among Chromosomes
#'
#' Finds the best matching regions between reference and query chromosomes based
#' on a coordinate-compressed BLAST-scoring matrix (see \code{buildSegScore}).
#'
#' @param pos Numeric matrix with scores, rows = queries, columns = segments.
#' @param seg.pos Integer vector, real start coordinate of each segment column.
#' @param seg.len Integer vector, length in bases of each segment column.
#' @param i.chr.ref.len Integer, length of the reference chromosome.
#' @param i.chr.ref Integer or string, ID of the reference chromosome.
#' @param i.chr.corresp Vector of identifiers for query chromosomes.
#' @param min.len Numeric, minimal length (in bases) for a valid region.
# @export
findBestChromosome <- function(pos,
                       seg.pos,
                       seg.len,
                       i.chr.ref.len,
                       i.chr.ref,
                       i.chr.corresp,
                       min.len, to.extend = T) {

  # Find maximum by chromosomes (column-wise arg-max over segments)
  pos.max <- pos[1, ]
  pos.i <- as.integer(pos.max != 0)
  if (nrow(pos) >= 2) {
    for (i in 2:nrow(pos)) {
      idx.i <- which(pos[i, ] > pos.max)
      pos.max[idx.i] <- pos[i, idx.i]
      pos.i[idx.i] <- i
    }
  }
  rm('pos')

  # Keep only segments with a winner; carry their real length and coordinates.
  idx.remain = pos.i != 0
  seg.i = pos.i[idx.remain]                       # winner per retained segment
  seg.w = seg.len[idx.remain]                     # bases per retained segment
  seg.b = seg.pos[idx.remain]                     # real start of retained segment
  seg.e = (seg.pos + seg.len - 1)[idx.remain]     # real end of retained segment
  if (length(seg.i) == 0) return(NULL)

  # Get regions (runs of a single winner over the retained-segment array).
  df.all = c()
  for(i in 1:length(i.chr.corresp)){
    df = findOnes((seg.i == i) * 1)               # beg/end are retained-seg indices
    if(nrow(df) == 0) next
    df$len = mapply(function(b, e) sum(seg.w[b:e]), df$beg, df$end)  # bases
    df$rbeg = seg.b[df$beg]
    df$rend = seg.e[df$end]
    df$i.ref = i.chr.ref
    df$i.acc =  i.chr.corresp[i]
    df.all = rbind(df.all, df)
  }
  df.all = df.all[order(df.all$beg),]             # spatial order (retained-seg index)

  # Gradually remove the smallest, merging same-query neighbours across the gap.
  while(min(df.all$len) < min.len){
    idx = which.min(df.all$len)[1]
    df.all = df.all[-idx,,drop=F]
    if(nrow(df.all) == 0) break

    idx.merge = which(diff(df.all$i.acc) == 0)
    if(length(idx.merge) == 0) next
    idx.merge = idx.merge[1]
    df.all$end[idx.merge]  = df.all$end[idx.merge + 1]
    df.all$rend[idx.merge] = df.all$rend[idx.merge + 1]
    df.all$len[idx.merge]  = sum(seg.w[df.all$beg[idx.merge]:df.all$end[idx.merge]])
    df.all = df.all[-(idx.merge + 1),,drop=F]
    if(nrow(df.all) == 0) stop("no correspondence is left")
  }
  if(nrow(df.all) == 0) return(NULL)

  # Map retained-segment indices back to real chromosome coordinates.
  df.all$beg = df.all$rbeg
  df.all$end = df.all$rend
  df.all$len = df.all$end - df.all$beg + 1
  df.all$rbeg = NULL
  df.all$rend = NULL

  if(to.extend){

    df.all = df.all[order(df.all$beg),]
    if(nrow(df.all) > 1){
      for (i in 1:(nrow(df.all)-1)) {
        gap <- df.all$beg[i+1] - df.all$end[i] - 1
        shift <- floor(gap / 2)

        df.all$end[i]   <- df.all$end[i] + shift
        df.all$beg[i+1] <- df.all$end[i] + 1
      }
    }

    df.all$beg[1] = 1
    df.all$end[nrow(df.all)] = i.chr.ref.len
    df.all$len = df.all$end - df.all$beg + 1
  }

  return(df.all)
}
