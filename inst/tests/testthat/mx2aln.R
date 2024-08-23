library(testthat)
source(system.file("utils", "utils.R", package = "pannagram"))


test_that("mx2aln: Couldn't correctly convert nucleotide matrix to sequences", {
  mx <- matrix(c('A', 'C', 'G', 'T', 'A', 'C', 'T', 'T', 'A', 'C', 'G', 'G', 'A', 'C', 'G', 'T',
                 'A', 'C', 'G', 'C', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'A', 'T',
                 'A', 'C', 'G', 'T', 'G', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'T', 'C', 'G', 'T'),
               nrow = 3, byrow = TRUE, dimnames = list(c('s1', 's2', 's3'), NULL))
  result <- mx2aln(mx)
  
  expect_type(result, "character")
  expect_length(result, 3)
  expect_named(result, c("s1", "s2", "s3"))
  expect_equal(result[["s1"]], "ACGTACTTACGGACGT")
  expect_equal(result[["s2"]], "ACGCACGTACGTACAT")
  expect_equal(result[["s3"]], "ACGTGCGTACGTTCGT")
})

test_that("mx2aln: Couldn't handle matrix without row names", {
  mx <- matrix(c('A', 'C', 'G', 'T', 'A', 'C', 'T', 'T', 'A', 'C', 'G', 'G', 'A', 'C', 'G', 'T',
                 'A', 'C', 'G', 'C', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'A', 'T',
                 'A', 'C', 'G', 'T', 'G', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'T', 'C', 'G', 'T'),
               nrow = 3, byrow = TRUE)
  result <- mx2aln(mx)
  
  expect_named(result, c("s1", "s2", "s3"))
  expect_equal(result[["s1"]], "ACGTACTTACGGACGT")
  expect_equal(result[["s2"]], "ACGCACGTACGTACAT")
  expect_equal(result[["s3"]], "ACGTGCGTACGTTCGT")
})

test_that("mx2aln: Couldn't handle empty matrix", {
  mx <- matrix(character(0), nrow = 0, ncol = 0)
  result <- mx2aln(mx)
  
  expect_type(result, "character")
  expect_length(result, 0)
})

test_that("mx2aln: Couldn't handle matrix with single row", {
  mx <- matrix(c('A', 'C', 'G', 'T', 'A', 'C', 'T', 'T', 'A', 'C', 'G', 'G', 'A', 'C', 'G', 'T'),
               nrow = 1, byrow = TRUE, dimnames = list(c('seq1'), NULL))
  result <- mx2aln(mx)
  
  expect_type(result, "character")
  expect_length(result, 1)
  expect_named(result, "seq1")
  expect_equal(result[["seq1"]], "ACGTACTTACGGACGT")
})

test_that("mx2aln: Couldn't handle matrix with single column", {
  mx <- matrix(c('A', 'C', 'G'), nrow = 3, byrow = TRUE, dimnames = list(c('s1', 's2', 's3'), NULL))
  result <- mx2aln(mx)
  
  expect_type(result, "character")
  expect_length(result, 3)
  expect_named(result, c("s1", "s2", "s3"))
  expect_equal(result[["s1"]], "A")
  expect_equal(result[["s2"]], "C")
  expect_equal(result[["s3"]], "G")
})