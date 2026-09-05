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

test_that("svTeOrder: LTR element", {
  ltr <- rndSeq(300, seed = 11)
  s <- paste0(ltr, rndSeq(1500, seed = 12), mutSeq(ltr, 0.03, seed = 13))

  res <- svTeOrder(s, name = "sv1")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_equal(res$sv, "sv1")
  expect_equal(res$te.class, "LTR")
  expect_gt(res$ltr.len, 250)
  expect_equal(res$tir.len, 0)
})

test_that("svTeOrder: TIR element", {
  tir <- rndSeq(24, seed = 21)
  s <- paste0(tir, rndSeq(1200, seed = 22), revCompl(tir))

  res <- svTeOrder(s)
  expect_equal(res$te.class, "TIR")
  expect_equal(res$tir.len, 24)
  expect_equal(res$tir.motif, substr(tir, 1, 6))
})

test_that("svTeOrder: LINE-like element, both orientations", {
  s <- paste0(rndSeq(2000, seed = 31), paste0(rep("A", 25), collapse = ""))

  res <- svTeOrder(s)
  expect_equal(res$te.class, "LINE_like")
  expect_equal(res$polyA.side, "right")

  res.rc <- svTeOrder(revCompl(s))
  expect_equal(res.rc$te.class, "LINE_like")
  expect_equal(res.rc$polyA.side, "left")
})

test_that("svTeOrder: CACTA gets no rule of its own", {
  s <- paste0("CACTA", rndSeq(2000, seed = 91), revCompl("CACTA"))

  # Five letters in a random context are not significant once the search has
  # looked in many places, and no motif is privileged to say otherwise.
  expect_null(findTIR(s))
  expect_equal(svTeOrder(s)$tir.method, "")
  expect_equal(svTeOrder(s)$note, "motif_CACTA")   # the motif is still annotated

  # With the terminal repeat that a real CACTA element carries, it is found
  tir <- rndSeq(14, seed = 92)
  s2 <- paste0("CACTA", tir, rndSeq(2000, seed = 93), revCompl(tir), revCompl("CACTA"))
  res <- svTeOrder(s2)
  expect_equal(res$te.class, "TIR")
  expect_gte(res$tir.len, 19)
  expect_equal(res$tir.method, "local")
})

test_that("svTeOrder: a long terminal repeat is found by the alignment", {
  tir <- rndSeq(24, seed = 95)
  s <- paste0(tir, rndSeq(1200, seed = 96), revCompl(tir))

  res <- svTeOrder(s)
  expect_equal(res$tir.len, 24)
  expect_equal(res$tir.method, "local")
})

test_that("svTeOrder: Helitron has no terminal repeats at all", {
  stem <- rndSeq(10, seed = 101)
  s <- paste0("TC", rndSeq(1500, seed = 102), "A", stem, "AAAA", revCompl(stem), "A",
              rndSeq(11, seed = 103), "CTAG")

  expect_null(findLTR(s))
  expect_null(findTIR(s))

  res <- svTeOrder(s)
  expect_equal(res$te.class, "Helitron")
  expect_equal(res$hel.stem, 10)
  expect_equal(res$hel.ctrr, "CTAG")
  expect_equal(res$polyA.side, "none")
  expect_true(grepl("helitron_tc", res$note))
})

test_that("svTeOrder: nothing is found in a random sequence", {
  res <- svTeOrder(rndSeq(2000, seed = 41))
  expect_equal(res$te.class, "unknown")
  expect_equal(res$polyA.side, "none")
  expect_equal(res$note, "")
})

test_that("svTeOrder: an A-rich terminus does not hide the LTR", {
  ltr <- rndSeq(300, seed = 51)
  s <- paste0(ltr, rndSeq(1500, seed = 52), mutSeq(ltr, 0.03, seed = 53),
              paste0(rep("A", 12), collapse = ""))

  res <- svTeOrder(s)
  expect_equal(res$te.class, "LTR")
  expect_true(grepl("polyA", res$note))
})

test_that("svTeOrder: thresholds are passed to the detectors", {
  ltr <- rndSeq(300, seed = 61)
  s <- paste0(ltr, rndSeq(1500, seed = 62), mutSeq(ltr, 0.03, seed = 63))

  res <- svTeOrder(s, ltr.par = list(min.len = 500))
  expect_equal(res$te.class, "unknown")
})
