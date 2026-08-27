# ============================================================================
#  Interval representation of the pangenome <-> accession correspondence.
#
#  LEGACY format (chunk_hdf5.R):
#    /accs/<acc> -- a vector of length(pangenome); the value is the signed
#    position in the accession chromosome, 0 means a gap.
#
#  INTERVAL format (this file):
#    /accs/<acc> -- an integer matrix [n x 5] with the columns
#        pan.beg | acc.beg | len | delta | block
#    acc(pan.beg + k) = acc.beg + delta * k,   k = 0 .. len-1,  acc.beg > 0
#    `delta` is +1/-1 and carries the strand, so no sign trick on `acc.beg`.
#    `block` is the id of the synteny block the interval belongs to.
#    /len -- the length of the pangenome (no longer implied by the dataset dim).
#
#  Relation to the legacy vector:   v[pan.beg + k] = delta * acc.beg + k
#  An interval is broken by ANY discontinuity, a 1 nt indel included.
# ============================================================================

ivl.cols <- c("pan.beg", "acc.beg", "len", "delta", "block")

#' Empty interval table
#' @export
ivlEmpty <- function(){
  data.frame(pan.beg = integer(0), acc.beg = integer(0), len = integer(0),
             delta = integer(0), block = integer(0))
}

#' Synteny blocks of an interval table
#'
#' Reproduces the block definition used by the legacy code (`findRuns` over the
#' sign-flipped ranks of abs(v)): two intervals adjacent along the pangenome stay
#' in one block if they have the same strand and their touching ends are
#' neighbours in the order of the accession coordinate.
#' @export
ivlBlocks <- function(ivl){
  m <- nrow(ivl)
  if(m == 0) return(integer(0))
  if(m == 1) return(1L)
  acc.end <- ivl$acc.beg + ivl$delta * (ivl$len - 1L)
  lo <- pmin(ivl$acc.beg, acc.end)
  o <- order(lo)
  cum <- integer(m); cum[o] <- cumsum(c(0L, ivl$len[o][-m]))   # elements before interval j
  rk <- function(x, j) cum[j] + (x - lo[j]) + 1L               # global rank of an element
  j <- 1:(m-1); d <- ivl$delta
  same <- (d[j] == d[j+1]) &
          (abs(d[j] * rk(acc.end[j], j) - d[j+1] * rk(ivl$acc.beg[j+1], j+1)) == 1L)
  cumsum(c(TRUE, !same))
}

#' Legacy vector -> interval table
#' @export
ivlFromVec <- function(v, blocks = TRUE){
  v[is.na(v)] <- 0
  n <- length(v)
  if(n == 0) return(ivlEmpty())
  nz <- v != 0
  if(!any(nz)) return(ivlEmpty())
  cont <- c(FALSE, nz[-1] & nz[-n] & (sign(v[-1]) == sign(v[-n])) & (diff(v) == 1))
  beg <- which(nz & !cont)
  end <- which(nz & !c(cont[-1], FALSE))
  # The table is stored as int32: guard the (currently theoretical) case of a
  # chromosome or a pangenome longer than 2^31-1 instead of silently writing NAs.
  if(max(n, max(abs(v))) > .Machine$integer.max)
    stop("interval_func: coordinates exceed int32; the interval table needs a wider type")
  ivl <- data.frame(pan.beg = as.integer(beg),
                    acc.beg = as.integer(abs(v[beg])),
                    len     = as.integer(end - beg + 1L),
                    delta   = as.integer(sign(v[beg])))
  ivl$block <- if(blocks) as.integer(ivlBlocks(ivl)) else 0L
  ivl
}

#' Interval table -> legacy vector
#' @export
ivlToVec <- function(ivl, len.pan){
  v <- numeric(len.pan)
  if(nrow(ivl) == 0) return(v)
  # expand without a per-row loop: sequence() gives the concatenated offsets
  pos <- sequence(ivl$len, from = ivl$pan.beg)
  off <- sequence(ivl$len, from = 0L)
  v[pos] <- rep(ivl$delta * ivl$acc.beg, ivl$len) + off
  v
}

#' Expand only a window [from, to] of the pangenome axis
#'
#' The point of the interval format: a region can be materialised without touching
#' the rest of the alignment. Returns a vector of length (to - from + 1).
#' @export
ivlToVecRange <- function(ivl, from, to){
  n <- to - from + 1L
  if(n <= 0) return(numeric(0))
  v <- numeric(n)
  if(nrow(ivl) == 0) return(v)
  pan.end <- ivl$pan.beg + ivl$len - 1L
  keep <- (pan.end >= from) & (ivl$pan.beg <= to)
  if(!any(keep)) return(v)
  d <- ivl[keep, , drop = FALSE]
  pe <- pan.end[keep]
  b <- pmax(d$pan.beg, from)          # clipped begin on the pangenome axis
  e <- pmin(pe, to)                   # clipped end
  acc.b <- d$delta * d$acc.beg + (b - d$pan.beg)   # legacy value at the clipped begin
  ln <- as.integer(e - b + 1L)
  pos <- sequence(ln, from = as.integer(b - from + 1L))
  v[pos] <- rep(acc.b, ln) + sequence(ln, from = 0L)
  v
}

#' Legacy value at given pangenome positions (0 where the accession has a gap)
#'
#' The inverse of ivlPanAt(). Point lookup by binary search: what `v[p]` was in the
#' legacy vector, without expanding anything.
#' @export
ivlAccAt <- function(ivl, p){
  res <- numeric(length(p))
  if(nrow(ivl) == 0 || length(p) == 0) return(res)
  o  <- order(ivl$pan.beg)
  pb <- ivl$pan.beg[o]; ln <- ivl$len[o]; ab <- ivl$acc.beg[o]; dl <- ivl$delta[o]
  i  <- findInterval(p, pb)
  ok <- (i > 0)
  ok[ok] <- p[ok] <= (pb[i[ok]] + ln[i[ok]] - 1L)
  j <- i[ok]
  res[ok] <- dl[j] * ab[j] + (p[ok] - pb[j])
  res
}

#' Accession-coordinate spans of the synteny blocks
#'
#' Mirrors the legacy `/blocks/<acc>` vector, which filled the WHOLE span of each
#' block in accession coordinates (gaps inside the block included).
#' @export
ivlBlockSpans <- function(ivl){
  if(nrow(ivl) == 0) return(data.frame(block = integer(0), lo = integer(0), hi = integer(0)))
  acc.end <- ivl$acc.beg + ivl$delta * (ivl$len - 1L)
  lo <- pmin(ivl$acc.beg, acc.end); hi <- pmax(ivl$acc.beg, acc.end)
  s <- data.frame(block = ivl$block, lo = lo, hi = hi)
  s <- do.call(rbind, lapply(split(s, s$block),
                             function(d) data.frame(block = d$block[1],
                                                    lo = min(d$lo), hi = max(d$hi))))
  s[order(s$lo), , drop = FALSE]
}

#' Block id at given accession coordinates (0 where not covered)
#'
#' Replacement for the legacy `b.acc[abs(x)]` lookup.
#' @export
ivlBlockAt <- function(ivl, q){
  res <- integer(length(q))
  if(nrow(ivl) == 0 || length(q) == 0) return(res)
  s <- ivlBlockSpans(ivl)
  i <- findInterval(q, s$lo)
  ok <- (i > 0)
  ok[ok] <- q[ok] <= s$hi[i[ok]]
  res[ok] <- s$block[i[ok]]
  res
}

#' Pangenome position of given accession coordinates (0 where not aligned)
#'
#' Replacement for the legacy `which(v == pos)` scan.
#' @export
ivlPanAt <- function(ivl, q){
  res <- integer(length(q))
  if(nrow(ivl) == 0 || length(q) == 0) return(res)
  acc.end <- ivl$acc.beg + ivl$delta * (ivl$len - 1L)
  lo <- pmin(ivl$acc.beg, acc.end); hi <- pmax(ivl$acc.beg, acc.end)
  o <- order(lo); lo <- lo[o]; hi <- hi[o]
  pb <- ivl$pan.beg[o]; ab <- ivl$acc.beg[o]; dl <- ivl$delta[o]
  i <- findInterval(q, lo)
  ok <- (i > 0)
  ok[ok] <- q[ok] <= hi[i[ok]]
  j <- i[ok]
  res[ok] <- pb[j] + dl[j] * (q[ok] - ab[j])
  res
}

# ---------------------------------------------------------------------------
# ---- HDF5 layer -----------------------------------------------------------
# ---------------------------------------------------------------------------

#' Names of the accessions stored in a file
#' @export
h5AccNames <- function(file, gr.b = "/accs"){
  g <- rhdf5::h5ls(file)
  g$name[g$group == gr.b]
}

#' Length of the pangenome stored in a file
#' @export
h5PanLen <- function(file){
  g <- rhdf5::h5ls(file, recursive = FALSE)
  if(!("len" %in% g$name)) stop(paste("No '/len' in", file))
  as.numeric(rhdf5::h5read(file, "len"))[1]
}

#' @export
h5PanLenSet <- function(file, len.pan){
  suppressMessages({
    try(rhdf5::h5delete(file, "len"), silent = TRUE)
    rhdf5::h5write(as.numeric(len.pan), file, "len")
  })
}

#' Write an interval table for one accession
#'
#' The length of the pangenome axis this table lives on is stored as the dataset
#' attribute `len.pan`. In the legacy format that length was carried by the dataset
#' dimension; keeping it per dataset (and not only in the file-level /len) is what
#' makes resume and the two-stage rewrites of comb_03 behave exactly as before.
#' @export
h5IvlWrite <- function(file, acc, ivl, len.pan, gr = "accs/"){
  m <- as.matrix(ivl[, ivl.cols, drop = FALSE])
  storage.mode(m) <- "integer"
  if(nrow(m) == 0) m <- matrix(0L, nrow = 1, ncol = length(ivl.cols))  # sentinel: len == 0
  s <- paste0(gr, acc)
  suppressMessages({
    try(rhdf5::h5delete(file, s), silent = TRUE)
    rhdf5::h5createDataset(file, s, dim(m), storage.mode = "integer",
                           chunk = c(min(65536L, nrow(m)), ncol(m)),
                           level = 6, filter = "gzip", shuffle = TRUE)
    rhdf5::h5write(m, file, s)
    fid <- rhdf5::H5Fopen(file)
    did <- rhdf5::H5Dopen(fid, s)
    rhdf5::h5writeAttribute(as.numeric(len.pan), did, "len.pan")
    rhdf5::H5Dclose(did); rhdf5::H5Fclose(fid)
  })
  invisible(NULL)
}

#' Read an interval table for one accession
#' @export
h5IvlRead <- function(file, acc, gr = "accs/"){
  m <- rhdf5::h5read(file, paste0(gr, acc))
  ivl <- data.frame(pan.beg = as.integer(m[, 1]), acc.beg = as.integer(m[, 2]),
                    len = as.integer(m[, 3]), delta = as.integer(m[, 4]),
                    block = as.integer(m[, 5]))
  ivl[ivl$len > 0, , drop = FALSE]
}

#' Length of the pangenome axis of one accession dataset
#' @export
h5AccLen <- function(file, acc, gr = "accs/"){
  a <- rhdf5::h5readAttributes(file, paste0(gr, acc))
  if(is.null(a$len.pan)) stop(paste("No 'len.pan' attribute for", acc, "in", file))
  as.numeric(a$len.pan)[1]
}

#' Write a legacy vector, storing it as an interval table
#' @export
h5VecWrite <- function(file, acc, v, gr = "accs/", blocks = TRUE){
  h5IvlWrite(file, acc, ivlFromVec(v, blocks = blocks), length(v), gr = gr)
}

#' Read an interval table and expand it into a legacy vector
#' @export
h5VecRead <- function(file, acc, len.pan = NULL, gr = "accs/"){
  if(is.null(len.pan)) len.pan <- h5AccLen(file, acc, gr = gr)
  ivlToVec(h5IvlRead(file, acc, gr = gr), len.pan)
}

#' Number of intervals of an accession without reading the payload
#' @export
h5IvlNrow <- function(file, acc, gr.b = "/accs"){
  g <- rhdf5::h5ls(file)
  d <- g$dim[(g$group == gr.b) & (g$name == acc)]
  if(length(d) == 0) return(NA_integer_)
  as.integer(strsplit(d[1], " x ")[[1]][1])
}
