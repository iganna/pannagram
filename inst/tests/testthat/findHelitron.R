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

# 5'-TC ... hairpin ... CTRR-3'; "A" around the stem is an anchor against
# an accidental extension of the stem by the random flanks
makeHelitron <- function(stem.len = 10, loop = "AAAA", tail = 11, ctrr = "CTAG", seed = 1) {
  stem <- rndSeq(stem.len, seed = seed)
  paste0("TC", rndSeq(1000, seed = seed + 1), "A", stem, loop, revCompl(stem), "A",
         rndSeq(tail, seed = seed + 2), ctrr)
}

test_that("findHelitron: finds the 3'-signature", {
  res <- findHelitron(makeHelitron(seed = 11))

  expect_type(res, "list")
  expect_equal(res$strand, "+")
  expect_equal(res$stem, 10)
  expect_equal(res$loop, 4)
  expect_equal(res$dist, 16)
  expect_equal(res$ctrr, "CTAG")
  expect_true(res$tc)
})

test_that("findHelitron: does not depend on the orientation of the SV", {
  s <- makeHelitron(seed = 21)

  res <- findHelitron(revCompl(s))
  expect_equal(res$strand, "-")
  expect_equal(res$stem, 10)
  expect_equal(res$ctrr, "CTAG")
})

test_that("findHelitron: all the variants of CTRR are accepted", {
  for (m in c("CTAA", "CTAG", "CTGA", "CTGG")) {
    expect_equal(findHelitron(makeHelitron(ctrr = m, seed = 31))$ctrr, m)
  }
  expect_null(findHelitron(makeHelitron(ctrr = "CTTT", seed = 31)))
})

test_that("findHelitron: the hairpin should be close to the terminus", {
  expect_type(findHelitron(makeHelitron(tail = 30, seed = 41)), "list")
  expect_null(findHelitron(makeHelitron(tail = 200, seed = 41)))
})

test_that("findHelitron: the stem should be long enough", {
  expect_null(findHelitron(makeHelitron(stem.len = 6, seed = 51)))
  expect_equal(findHelitron(makeHelitron(stem.len = 6, seed = 51), min.stem = 6)$stem, 6)
})

test_that("findHelitron: TC at the 5'-terminus is reported, but not required", {
  s <- makeHelitron(seed = 61)
  s.no.tc <- paste0("GG", substr(s, 3, nchar(s)))

  res <- findHelitron(s.no.tc)
  expect_type(res, "list")
  expect_false(res$tc)
})

test_that("findHelitron: nothing in a random sequence", {
  expect_null(findHelitron(rndSeq(2000, seed = 71)))
})
