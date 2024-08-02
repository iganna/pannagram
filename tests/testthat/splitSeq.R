library(testthat)
source("../../utils/utils.R")

test_that("splitSeq: Could't split sequence into chunks of specified length", {
  seq <- "1111222233334444"
  result <- splitSeq(seq, n = 4)
  
  expect_type(result, "character")
  expect_length(result, 4)
  expect_equal(result, c("1111", "2222", "3333", "4444"))
})

test_that("splitSeq: Couldn't handle sequences shorter than chunk length", {
  seq <- "ACGT"
  result <- splitSeq(seq, n = 5)
  
  expect_type(result, "character")
  expect_length(result, 1)
  expect_equal(result, "ACGT")
})

test_that("splitSeq: Didn't pad the last chunk with empty characters", {
  seq <- "ACGTAnna0123456789"
  result <- splitSeq(seq, n = 5)
  
  expect_type(result, "character")
  expect_length(result, 4)
  expect_equal(result, c("ACGTA", "nna01", "23456", "789"))
})

test_that("splitSeq: Couldn't handle empty sequence", {
  seq <- ""
  result <- splitSeq(seq, n = 5)
  
  expect_type(result, "character")
  expect_length(result, 0)
  expect_equal(result, character(0))
})

test_that("splitSeq: Couldn't handle step parameter", {
  seq <- "ACGTAnna0123456789"
  result <- splitSeq(seq, n = 4, step = 2)
  
  expect_type(result, "character")
  expect_length(result, 4)
  expect_equal(result, c("GTAn", "na01", "2345", "6789"))
})