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

test_that("findTIR: finds a short perfect inverted repeat", {
  tir <- rndSeq(20, seed = 11)
  s <- paste0(tir, rndSeq(1000, seed = 12), revCompl(tir))

  res <- findTIR(s)
  expect_type(res, "list")
  expect_equal(res$len, 20)
  expect_equal(res$ident, 1)
})

test_that("findTIR: finds a long diverged inverted repeat", {
  tir <- rndSeq(250, seed = 21)
  s <- paste0(tir, rndSeq(2000, seed = 22), revCompl(mutSeq(tir, 0.05, seed = 23)))

  res <- findTIR(s)
  expect_type(res, "list")
  expect_gt(res$len, 150)
  expect_gt(res$ident, 0.85)
})

test_that("findTIR: no repeat in a random sequence", {
  expect_null(findTIR(rndSeq(2000, seed = 31)))
})

test_that("findTIR: poly-T / poly-A termini are not TIRs", {
  s <- paste0(paste0(rep("T", 20), collapse = ""),
              rndSeq(2000, seed = 41),
              paste0(rep("A", 20), collapse = ""))

  expect_null(findTIR(s))
})

test_that("findTIR: direct terminal repeats are not inverted ones", {
  ltr <- rndSeq(300, seed = 51)
  s <- paste0(ltr, rndSeq(1500, seed = 52), ltr)

  expect_null(findTIR(s))
})

test_that("findTIR: the two termini of an SV may be shifted differently", {
  tir <- rndSeq(25, seed = 61)
  s <- paste0(rndSeq(4, seed = 62), tir, rndSeq(1000, seed = 63), revCompl(tir), rndSeq(7, seed = 64))

  # the extension is anchored, but the chaining tolerates the shift by itself
  res <- findTIR(s)
  expect_type(res, "list")
  expect_equal(res$method, "local")
  expect_gte(res$len, 20)
  expect_equal(res$beg, 5)                 # the hit starts after the 4 extra nucleotides

  expect_null(findTIR(s, max.offset = 2))  # the shift is larger than allowed

})

test_that("findTIR: a strong internal inverted repeat does not hide the terminal one", {
  tir <- rndSeq(20, seed = 71)
  big <- rndSeq(200, seed = 72)
  s <- paste0(tir, rndSeq(300, seed = 73), big, rndSeq(300, seed = 74),
              revCompl(big), rndSeq(300, seed = 75), revCompl(tir))

  res <- findTIR(s)
  expect_type(res, "list")
  expect_gte(res$len, 20)
})
