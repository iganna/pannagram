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

test_that("findTermMotif: finds the CACTA signature", {
  s <- paste0("CACTA", rndSeq(1000, seed = 11), revCompl("CACTA"))

  res <- findTermMotif(s)
  expect_type(res, "list")
  expect_equal(res$motif, "CACTA")
  expect_equal(res$offset, 0)
  expect_equal(res$mism, 0)
})

test_that("findTermMotif: does not depend on the orientation of the SV", {
  s <- paste0("CACTA", rndSeq(1000, seed = 21), revCompl("CACTA"))

  expect_equal(findTermMotif(revCompl(s))$motif, "CACTA")
})

test_that("findTermMotif: tolerates inaccurate SV boundaries", {
  s <- paste0("GT", "CACTA", rndSeq(1000, seed = 31), revCompl("CACTA"), "AC")

  res <- findTermMotif(s)
  expect_equal(res$motif, "CACTA")
  expect_equal(res$offset, 2)

  expect_null(findTermMotif(s, max.offset = 1))
})

test_that("findTermMotif: both termini are required", {
  expect_null(findTermMotif(paste0("CACTA", rndSeq(1000, seed = 41))))
  expect_null(findTermMotif(paste0(rndSeq(1000, seed = 42), revCompl("CACTA"))))
})

test_that("findTermMotif: nothing in a random sequence", {
  expect_null(findTermMotif(rndSeq(1000, seed = 51)))
})

test_that("findTermMotif: an empty list of motifs switches the check off", {
  s <- paste0("CACTA", rndSeq(1000, seed = 61), revCompl("CACTA"))

  expect_null(findTermMotif(s, motifs = character(0)))
})

test_that("findTermMotif: mismatches are allowed on demand", {
  s <- paste0("CACTT", rndSeq(1000, seed = 71), revCompl("CACTA"))

  expect_null(findTermMotif(s))
  expect_equal(findTermMotif(s, max.mism = 1)$mism, 1)
})
