library(testthat)
source(system.file("utils", "utils.R", package = "pannagram"))


test_that("aln2mx: Couldn't correctly convert aligned sequences to nucleotide matrix", {
  s.aln <- c("ACGTACTTACGGACGT", "ACGCACGTACGTACAT", "ACGTGCGTACGTTCGT")
  result <- aln2mx(s.aln)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 16))
  expect_equal(rownames(result), c("s1", "s2", "s3"))
  expect_equal(result[1, ], seq2nt("ACGTACTTACGGACGT"))
  expect_equal(result[2, ], seq2nt("ACGCACGTACGTACAT"))
  expect_equal(result[3, ], seq2nt("ACGTGCGTACGTTCGT"))
})

test_that("aln2mx: Couldn't handle sequences with different lengths", {
  s.aln <- c("ACGTACTTACGGACGT", "ACGCACGTACGTACAT", "ACGTGCGTACGTT")
  expect_error(aln2mx(s.aln), "Aligned sequences have different lengths")
})

test_that("aln2mx: Couldn't handle sequences without names", {
  s.aln <- c("ACGTACTTACGGACGT", "ACGCACGTACGTACAT", "ACGTGCGTACGTTCGT")
  names(s.aln) <- NULL
  result <- aln2mx(s.aln)
  
  expect_equal(rownames(result), c("s1", "s2", "s3"))
})

test_that("aln2mx: Couldn't handle sequences with names", {
  s.aln <- c("ACGTACTTACGGACGT", "ACGCACGTACGTACAT", "ACGTGCGTACGTTCGT")
  names(s.aln) <- c("seq1", "seq2", "seq3")
  result <- aln2mx(s.aln)
  
  expect_equal(rownames(result), c("seq1", "seq2", "seq3"))
})

test_that("aln2mx: Couldn't handle empty input", {
  s.aln <- character(0)
  result <- aln2mx(s.aln)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(0, 0))
})
