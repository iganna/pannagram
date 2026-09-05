library(testthat)
source(system.file("utils/utils.R", package = "pannagram"))

# the copy in the repo wins over the installed one: it may not be rebuilt yet
path.func <- Filter(file.exists, c("../../analys/sv_aln_func.R", "inst/analys/sv_aln_func.R"))
path.func <- if (length(path.func)) path.func[1] else system.file("analys/sv_aln_func.R", package = "pannagram")
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

test_that("svFamilyRepr: takes the longest ones", {
  len <- c(500, 1000, 800, 300, 900)

  expect_equal(len[svFamilyRepr(len, n.repr = 3)], c(1000, 900, 800))
  expect_length(svFamilyRepr(len, n.repr = 10), 5)
})

test_that("svFamilyRepr: only the members with a mutual edge", {
  len <- c(1000, 900, 800, 700, 600)
  mut <- c(FALSE, FALSE, TRUE, TRUE, TRUE)

  expect_equal(len[svFamilyRepr(len, n.repr = 3, mutual = mut)], c(800, 700, 600))
})

test_that("svFamilyRepr: without mutual edges the family is left as it is", {
  len <- c(1000, 900, 800)

  expect_equal(len[svFamilyRepr(len, n.repr = 2, mutual = c(FALSE, FALSE, FALSE))], c(1000, 900))
  expect_equal(len[svFamilyRepr(len, n.repr = 2, mutual = c(FALSE, TRUE, FALSE))], c(1000, 900))
})

test_that("svFamilyRepr: one member", {
  expect_equal(svFamilyRepr(c(500), n.repr = 20), 1)
})
