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

test_that("findPolyA: poly-A at the right terminus", {
  s <- paste0(rndSeq(1000, seed = 11), paste0(rep("A", 25), collapse = ""))

  res <- findPolyA(s)
  expect_equal(res$side, "right")
  expect_gte(res$len, 25)
})

test_that("findPolyA: poly-T at the left terminus", {
  s <- paste0(paste0(rep("T", 25), collapse = ""), rndSeq(1000, seed = 21))

  res <- findPolyA(s)
  expect_equal(res$side, "left")
  expect_gte(res$len, 25)
})

test_that("findPolyA: an interrupted tail is still a tail", {
  s <- paste0(rndSeq(1000, seed = 31), "AAAAAAAAAAGAAAAAAAAAA")

  res <- findPolyA(s)
  expect_equal(res$side, "right")
  expect_gte(res$len, 20)
  expect_gt(res$ident, 0.9)
})

test_that("findPolyA: tails at both termini are ambiguous", {
  s <- paste0(paste0(rep("T", 20), collapse = ""),
              rndSeq(1000, seed = 41),
              paste0(rep("A", 20), collapse = ""))

  expect_equal(findPolyA(s)$side, "both")
})

test_that("findPolyA: an internal run is not a tail", {
  s <- paste0(rndSeq(500, seed = 51), paste0(rep("A", 30), collapse = ""), rndSeq(500, seed = 52))

  res <- findPolyA(s)
  expect_equal(res$side, "none")
  expect_equal(res$len, 0)
})

test_that("findPolyA: no tail in a random sequence", {
  expect_equal(findPolyA(rndSeq(1000, seed = 61))$side, "none")
})
