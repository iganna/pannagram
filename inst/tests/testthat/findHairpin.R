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

# "A" on both sides of the stem is an anchor: A is not complementary to A,
# so the stem cannot be accidentally extended by the random flanks
test_that("findHairpin: finds a hairpin near the 3'-terminus", {
  stem <- rndSeq(10, seed = 11)
  s <- paste0(rndSeq(500, seed = 12), "A", stem, "AAAA", revCompl(stem), "A", rndSeq(15, seed = 13))

  res <- findHairpin(s)
  expect_type(res, "list")
  expect_equal(res$stem, 10)
  expect_equal(res$loop, 4)
  expect_equal(res$dist, 16)
  expect_equal(substr(s, res$beg, res$beg + res$stem - 1), stem)
  expect_equal(substr(s, res$end - res$stem + 1, res$end), revCompl(stem))
})

test_that("findHairpin: works at the 5'-terminus too", {
  stem <- rndSeq(9, seed = 21)
  s <- paste0(rndSeq(11, seed = 22), "A", stem, "ACGTA", revCompl(stem), "A", rndSeq(500, seed = 23))

  res <- findHairpin(s, side = "left")
  expect_equal(res$stem, 9)
  expect_equal(res$loop, 5)
  expect_equal(res$dist, 12)
})

test_that("findHairpin: a hairpin outside the window is not found", {
  stem <- rndSeq(10, seed = 31)
  s <- paste0(rndSeq(500, seed = 32), "A", stem, "AAAA", revCompl(stem), "A", rndSeq(300, seed = 33))

  expect_null(findHairpin(s))              # the default window is 60 nt
  expect_type(findHairpin(s, w = 400), "list")
})

test_that("findHairpin: a mismatch breaks the stem", {
  stem <- rndSeq(12, seed = 41)
  arm <- seq2nt(revCompl(stem))
  arm[3] <- c(A = "C", C = "A", G = "T", T = "G")[arm[3]]   # a mismatch 3 nt from the loop

  s.ok <- paste0(rndSeq(200, seed = 42), "A", stem, "AAA", revCompl(stem), "A", rndSeq(10, seed = 43))
  s.mm <- paste0(rndSeq(200, seed = 42), "A", stem, "AAA", nt2seq(arm), "A", rndSeq(10, seed = 43))

  expect_equal(findHairpin(s.ok, min.stem = 12, loop.min = 3, loop.max = 3)$stem, 12)
  expect_null(findHairpin(s.mm, min.stem = 12, loop.min = 3, loop.max = 3))

  # a free loop is allowed to absorb the mismatch and restores most of the stem
  expect_gte(findHairpin(s.mm, min.stem = 2)$stem, 9)
})

test_that("findHairpin: no long stem in a random sequence", {
  # a stem of 7 occurs by chance, a long one does not
  expect_null(findHairpin(rndSeq(500, seed = 51), min.stem = 10))
})

test_that("findHairpin: too short sequence", {
  expect_null(findHairpin(rndSeq(15, seed = 61)))
})
