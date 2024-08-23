library(testthat)
source(system.file("utils", "utils.R", package = "pannagram"))

test_that("mx2seq: Couldn't correctly convert nucleotide matrix to sequences without gaps", {
  mx <- matrix(c('A', 'C', 'G', 'T', '-', 'C', 'T', 'T', 'A', '-', 'G', 'G', 'A', 'C', '-', 'T',
                 'A', 'C', 'G', 'C', 'A', 'C', 'G', 'T', 'A', '-', 'G', 'T', 'A', 'C', '-', 'T',
                 'A', 'C', 'G', 'T', 'G', 'C', 'G', 'T', '-', 'C', 'G', 'T', 'T', 'C', '-', 'T'),
               nrow = 3, byrow = TRUE, dimnames = list(c('s1', 's2', 's3'), NULL))
  result <- mx2seq(mx)
  
  expect_type(result, "character")
  expect_length(result, 3)
  expect_named(result, c("s1", "s2", "s3"))
  expect_equal(result[["s1"]], "ACGTCTTAGGACT")
  expect_equal(result[["s2"]], "ACGCACGTAGTACT")
  expect_equal(result[["s3"]], "ACGTGCGTCGTTCT")
})

test_that("mx2seq: Couldn't handle matrix without row names", {
  mx <- matrix(c('A', 'C', 'G', 'T', '-', 'C', 'T', 'T', 'A', '-', 'G', 'G', 'A', 'C', '-', 'T',
                 'A', 'C', 'G', 'C', 'A', 'C', 'G', 'T', 'A', '-', 'G', 'T', 'A', 'C', '-', 'T',
                 'A', 'C', 'G', 'T', 'G', 'C', 'G', 'T', '-', 'C', 'G', 'T', 'T', 'C', '-', 'T'),
               nrow = 3, byrow = TRUE)
  result <- mx2seq(mx)
  
  expect_named(result, c("s1", "s2", "s3"))
  expect_equal(result[["s1"]], "ACGTCTTAGGACT")
  expect_equal(result[["s2"]], "ACGCACGTAGTACT")
  expect_equal(result[["s3"]], "ACGTGCGTCGTTCT")
})

test_that("mx2seq: Couldn't handle empty matrix", {
  mx <- matrix(character(0), nrow = 0, ncol = 0)
  result <- mx2seq(mx)
  
  expect_type(result, "character")
  expect_length(result, 0)
})

test_that("mx2seq: Couldn't handle matrix with single row", {
  mx <- matrix(c('A', 'C', 'G', 'T', '-', 'C', 'T', 'T', 'A', '-', 'G', 'G', 'A', 'C', '-', 'T'),
               nrow = 1, byrow = TRUE, dimnames = list(c('seq1'), NULL))
  result <- mx2seq(mx)
  
  expect_type(result, "character")
  expect_length(result, 1)
  expect_named(result, "seq1")
  expect_equal(result[["seq1"]], "ACGTCTTAGGACT")
})

test_that("mx2seq: Couldn't handle matrix with single column", {
  mx <- matrix(c('A', '-', 'G'), nrow = 3, byrow = TRUE, dimnames = list(c('s1', 's2', 's3'), NULL))
  result <- mx2seq(mx)
  
  expect_type(result, "character")
  expect_length(result, 3)
  expect_named(result, c("s1", "s2", "s3"))
  expect_equal(result[["s1"]], "A")
  expect_equal(result[["s2"]], "")
  expect_equal(result[["s3"]], "G")
})