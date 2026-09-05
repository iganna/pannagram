library(testthat)
source(system.file("utils/utils.R", package = "pannagram"))

# the copy in the repo wins over the installed one: it may not be rebuilt yet
path.func <- Filter(file.exists, c("../../analys/sv_te_func.R", "inst/analys/sv_te_func.R"))
path.func <- if (length(path.func)) path.func[1] else system.file("analys/sv_te_func.R", package = "pannagram")
source(path.func)

# Deterministic sequences, independent of the RNG of R
lcg <- function(n, seed = 1) {
  x <- seed
  v <- numeric(n)
  for (i in seq_len(n)) {
    x <- (69069 * x + 1) %% 2147483648
    v[i] <- x / 2147483648
  }
  v
}
rndSeq <- function(n, seed = 1) paste0(c("A", "C", "G", "T")[floor(lcg(n, seed) * 4) + 1], collapse = "")
mutSeq <- function(s, p, seed = 1) {
  v <- seq2nt(s)
  u <- lcg(2 * length(v), seed)
  idx <- which(u[seq_along(v)] < p)
  v[idx] <- c("A", "C", "G", "T")[floor(u[length(v) + idx] * 4) + 1]
  nt2seq(v)
}

test_that("findLTR: finds terminal direct repeats", {
  ltr <- rndSeq(300, seed = 11)
  s <- paste0(ltr, rndSeq(1500, seed = 12), mutSeq(ltr, 0.03, seed = 13))

  res <- findLTR(s)
  expect_type(res, "list")
  expect_gt(res$len, 250)
  expect_lte(res$len, 300)
  expect_gt(res$ident, 0.85)
  expect_lte(res$beg.left, 20)  # the hit is anchored at the terminus
})

test_that("findLTR: tolerates an indel inside the LTR", {
  ltr <- rndSeq(400, seed = 21)
  ltr.2 <- paste0(substr(ltr, 1, 200), rndSeq(7, seed = 22), substr(ltr, 201, 400))
  s <- paste0(ltr, rndSeq(1500, seed = 23), mutSeq(ltr.2, 0.03, seed = 24))

  res <- findLTR(s)
  expect_type(res, "list")
  expect_gt(res$len, 300)
})

test_that("findLTR: no repeat in a random sequence", {
  expect_null(findLTR(rndSeq(3000, seed = 31)))
})

test_that("findLTR: internal repeats are not terminal", {
  rpt <- rndSeq(300, seed = 41)

  # within max.offset the pair still counts: the boundaries of an SV are inaccurate
  s.near <- paste0(rndSeq(200, seed = 42), rpt, rndSeq(1500, seed = 43), rpt, rndSeq(200, seed = 44))
  expect_type(findLTR(s.near), "list")
  expect_equal(findLTR(s.near)$off, 200)

  # further away it is an internal duplication, not a terminal repeat
  s.far <- paste0(rndSeq(600, seed = 42), rpt, rndSeq(1500, seed = 43), rpt, rndSeq(600, seed = 44))
  expect_null(findLTR(s.far))
})

test_that("findLTR: low-complexity termini are not LTRs", {
  a.rich <- paste0(rep("A", 200), collapse = "")
  s <- paste0(a.rich, rndSeq(1500, seed = 51), a.rich)

  expect_null(findLTR(s))
})

test_that("findLTR: too short sequence", {
  expect_null(findLTR(rndSeq(150, seed = 61)))
})
