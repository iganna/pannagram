library(testthat)

# The shared runner (tests/testthat/test.R) only sources utils/utils.R, so pull in
# the interval codec explicitly when it is not already on the search path.
if(!exists("findRuns"))    source(system.file("utils/utils.R", package = "pannagram"))
if(!exists("ivlFromVec"))  source(system.file("utils/interval_func.R", package = "pannagram"))

# Reference implementations copied from the legacy (vector) code path, so the
# interval codec is checked against what pannagram actually did before.
.refBlocksAcc <- function(v){                       # comb_03 / defineBlocks
  v.idx <- seq_along(v); v.idx <- v.idx[v != 0]; w <- v[v != 0]
  if(!length(w)) return(integer(0))
  r <- rank(abs(w)); r[w < 0] <- -r[w < 0]
  b <- findRuns(r)
  b$v.beg <- w[b$beg]; b$v.end <- w[b$end]
  out <- rep(0L, max(abs(w)))
  lo <- pmin(abs(b$v.beg), abs(b$v.end)); hi <- pmax(abs(b$v.beg), abs(b$v.end))
  out[sequence(hi - lo + 1L, from = lo)] <- rep(seq_len(nrow(b)), hi - lo + 1L)
  out
}

.randomVec <- function(L){
  v <- rep(0, L); p <- 1; a <- sample(1:40, 1)
  while(p <= L - 3){
    run <- sample(1:12, 1); if(p + run - 1 > L) break
    if(runif(1) < 0.35) v[p:(p+run-1)] <- -((a+run-1):a) else v[p:(p+run-1)] <- a:(a+run-1)
    a <- a + run + sample(0:2, 1)      # 0 = joint, 1 = a 1 nt insertion in the accession
    p <- p + run + sample(0:2, 1)      # 0 = joint, 1 = a 1 nt deletion in the accession
  }
  v
}

test_that("ivlFromVec: a toy pangenome with inversions and 1 nt indels", {
  v <- rep(0, 46)
  v[1:10]  <- 101:110
  v[12:16] <- 111:115          # 1 nt deletion in the accession at pan 11
  v[17:20] <- 117:120          # 1 nt insertion in the accession (116 is skipped)
  v[21:26] <- -(140:135)       # inversion
  v[31]    <- 121
  v[33:38] <- 122:127

  ivl <- ivlFromVec(v)
  expect_equal(ivl$pan.beg, c(1L, 12L, 17L, 21L, 31L, 33L))
  expect_equal(ivl$acc.beg, c(101L, 111L, 117L, 140L, 121L, 122L))
  expect_equal(ivl$len,     c(10L, 5L, 4L, 6L, 1L, 6L))
  expect_equal(ivl$delta,   c(1L, 1L, 1L, -1L, 1L, 1L))
  # every SV breaks the interval, but the block keeps the pieces together
  expect_equal(ivl$block,   c(1L, 1L, 1L, 2L, 3L, 3L))
})

test_that("ivlFromVec: an inversion stays a single row", {
  v <- rep(0, 20); v[5:12] <- -(200:193)
  ivl <- ivlFromVec(v)
  expect_equal(nrow(ivl), 1L)
  expect_equal(ivl$delta, -1L)
  expect_equal(ivl$acc.beg, 200L)
  expect_equal(ivl$len, 8L)
})

test_that("ivlFromVec: empty and all-zero input", {
  expect_equal(nrow(ivlFromVec(numeric(0))), 0L)
  expect_equal(nrow(ivlFromVec(rep(0, 100))), 0L)
  expect_equal(ivlToVec(ivlFromVec(rep(0, 100)), 100), rep(0, 100))
})

test_that("ivlToVec: round-trip is exact on random alignments", {
  set.seed(42)
  for(i in 1:200){
    L <- sample(50:250, 1); v <- .randomVec(L)
    expect_identical(as.numeric(ivlToVec(ivlFromVec(v), L)), as.numeric(v))
  }
})

test_that("ivlBlockAt: matches the legacy /blocks vector element by element", {
  set.seed(43)
  for(i in 1:200){
    L <- sample(50:250, 1); v <- .randomVec(L)
    if(anyDuplicated(abs(v[v != 0]))) next      # the pipeline de-duplicates first
    ref <- .refBlocksAcc(v)
    if(!length(ref)) next
    expect_identical(as.integer(ivlBlockAt(ivlFromVec(v), seq_along(ref))), as.integer(ref))
  }
})

test_that("ivlPanAt: matches the legacy which(v == pos) scan", {
  set.seed(44)
  for(i in 1:200){
    L <- sample(50:250, 1); v <- .randomVec(L)
    if(anyDuplicated(abs(v[v != 0]))) next
    ivl <- ivlFromVec(v); w <- v[v != 0]
    ref <- vapply(w, function(x) which(v == x)[1], integer(1))
    expect_identical(as.integer(ivlPanAt(ivl, abs(w))), as.integer(ref))
  }
})

test_that("ivlToVecRange: a window equals the corresponding slice of the full vector", {
  set.seed(45)
  for(i in 1:100){
    L <- sample(60:300, 1); v <- .randomVec(L); ivl <- ivlFromVec(v)
    for(k in 1:5){
      b <- sample(L, 1); e <- b + sample(0:(L - b), 1)
      expect_identical(as.numeric(ivlToVecRange(ivl, b, e)), as.numeric(v[b:e]))
    }
  }
})

test_that("ivlAccAt: point lookup equals the legacy v[p]", {
  set.seed(46)
  for(i in 1:200){
    L <- sample(50:250, 1); v <- .randomVec(L); ivl <- ivlFromVec(v)
    p <- sample(L, min(L, 40))
    expect_identical(as.numeric(ivlAccAt(ivl, p)), as.numeric(v[p]))
    expect_identical(as.numeric(ivlAccAt(ivl, seq_len(L))), as.numeric(v))
  }
  expect_equal(length(ivlAccAt(ivlEmpty(), 1:5)), 5L)
  expect_true(all(ivlAccAt(ivlEmpty(), 1:5) == 0))
})

test_that("ivlAccAt / ivlPanAt are inverse on aligned positions", {
  set.seed(47)
  for(i in 1:100){
    L <- sample(50:250, 1); v <- .randomVec(L)
    if(anyDuplicated(abs(v[v != 0]))) next
    ivl <- ivlFromVec(v)
    p <- which(v != 0)
    expect_identical(as.integer(ivlPanAt(ivl, abs(ivlAccAt(ivl, p)))), as.integer(p))
  }
})

test_that("ivlBlockAt: 0 between blocks is not the same as past the end", {
  # Regression: filterBlocks() keeps a feature whose both ends fall between two
  # blocks (legacy value 0 on both sides) and drops it only when a position runs
  # past max(abs(v)) (legacy NA). ivlBlockAt() returns 0 for both, so callers must
  # tell them apart by the accession length -- this pins the two cases down.
  v <- rep(0, 30)
  v[1:5]   <- 101:105        # block A: acc 101..105, plus strand
  v[11:15] <- -(310:306)     # block B: acc 306..310, minus strand
  ivl <- ivlFromVec(v)
  ref <- defineBlocks(v)     # the legacy vector, length max(abs(v)) == 310

  expect_equal(length(ref), 310L)
  expect_equal(ivlBlockAt(ivl, 103), ref[103])     # inside block A
  expect_equal(ivlBlockAt(ivl, 308), ref[308])     # inside block B
  expect_equal(ivlBlockAt(ivl, 200), 0L)           # between the blocks
  expect_equal(ref[200], 0)                        # ... and the legacy agrees
  expect_equal(ivlBlockAt(ivl, 50), 0L)            # before the first block
  expect_equal(ref[50], 0)
  expect_true(is.na(ref[400]))                     # past the end: legacy NA
  expect_equal(ivlBlockAt(ivl, 400), 0L)           # ... but ivlBlockAt says 0

  # The whole vector, position by position, with the two cases distinguished
  q <- seq_along(ref)
  got <- ivlBlockAt(ivl, q)
  expect_identical(as.integer(got), as.integer(ref))
})

test_that("ivlFromVec: refuses coordinates that do not fit into int32", {
  v <- c(1, 2, .Machine$integer.max + 10)
  expect_error(ivlFromVec(v), "int32")
})
