library(testthat)
source(system.file("utils", "utils.R", package = "pannagram"))


test_that("seq2nt: Couldn't convert string to vector of nucleotides", {
  result <- seq2nt("ATGC")
  expect_type(result, "character")
  expect_length(result, 4)
  expect_equal(result, c("A", "T", "G", "C"))
})

test_that("seq2nt: Couldn't handle empty string", {
  result <- seq2nt("")
  expect_type(result, "character")
  expect_length(result, 0)
  expect_equal(result, character(0))
})

test_that("seq2nt: Didn't throw error for multiple sequences", {
  expect_error(seq2nt(c("ATGC", "CGTA")), "There should be only one sequence")
})