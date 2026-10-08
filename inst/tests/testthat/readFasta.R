library(testthat)
source(system.file("utils", "utils.R", package = "pannagram"))

test_that("readFasta: Couldn't correctly read FASTA file", {
  test_file <- tempfile(fileext = ".fasta")
  writeLines(c(
    ">seq1",
    "ATGCATGC",
    ">seq2",
    "GCATGCAT",
    ">seq3",
    "TGCATGCA"
  ), test_file)
  
  result <- readFasta(test_file)
  
  expect_type(result, "character")
  expect_length(result, 3)
  expect_named(result, c("seq1", "seq2", "seq3"))
  expect_equal(result[[1]], "ATGCATGC")
  expect_equal(result[[2]], "GCATGCAT")
  expect_equal(result[[3]], "TGCATGCA")
  
  unlink(test_file)
})

test_that("readFasta: Couldn't handles multiline sequences", {
  multiline_file <- tempfile(fileext = ".fasta")
  writeLines(c(
    ">seq1",
    "ATGC",
    "ATGC",
    ">seq2",
    "GCAT",
    "GCAT"
  ), multiline_file)
  
  result <- readFasta(multiline_file)
  
  expect_type(result, "character")
  expect_length(result, 2)
  expect_named(result, c("seq1", "seq2"))
  expect_equal(result[["seq1"]], "ATGCATGC")
  expect_equal(result[["seq2"]], "GCATGCAT")
  
  unlink(multiline_file)
})

test_that("readFasta: Couldn't handle empty file", {
  empty_file <- tempfile(fileext = ".fasta")
  writeLines(character(0), empty_file)

  result <- readFasta(empty_file, stop.on.error=FALSE)
  expect_type(result, "character")
  expect_length(result, 0)
  
  expect_error(readFasta(empty_file, stop.on.error=TRUE))
  unlink(empty_file)
})

test_that("readFasta: Couldn't handle file without headers", {
  no_headers_file <- tempfile(fileext = ".fasta")
  writeLines(c("ATGC", "GCTA"), no_headers_file)

  result <- readFasta(no_headers_file, stop.on.error=FALSE)
  expect_type(result, "character")
  expect_length(result, 0)

  expect_error(readFasta(no_headers_file, stop.on.error=TRUE))
  
  unlink(no_headers_file)
})

test_that("readFasta: Couldn't handle file with only headers", {
  headers_file <- tempfile(fileext = ".fasta")
  writeLines(c(">seq1", ">seq2", ">seq3"), headers_file)
  result <- readFasta(headers_file, stop.on.error=FALSE)
  
  expect_type(result, "character")
  expect_length(result, 3)
  expect_named(result, c("seq1", "seq2", "seq3"))
  expect_equal(result[["seq1"]], "")
  expect_equal(result[["seq2"]], "")
  expect_equal(result[["seq3"]], "")

  expect_error(readFasta(headers_file, stop.on.error = TRUE))
  
  unlink(headers_file)
})