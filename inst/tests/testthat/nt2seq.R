library(testthat)
source(system.file("utils", "utils.R", package = "pannagram"))

test_that("nt2seq: Converts vector of nucleotides to string", {
  result <- nt2seq(c("A", "T", "G", "C"))
  expect_type(result, "character")
  expect_length(result, 1)
  expect_equal(result, "ATGC")
})

test_that("nt2seq: Handles empty vector", {
  result <- nt2seq(character(0))
  expect_type(result, "character")
  expect_length(result, 1)
  expect_equal(result, "")
})

test_that("nt2seq: Handles single nucleotide", {
  result <- nt2seq(c("A"))
  expect_type(result, "character")
  expect_length(result, 1)
  expect_equal(result, "A")
})