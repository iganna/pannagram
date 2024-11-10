library(testthat)
library(crayon)

# Capture output for testing (reuse from previous example)
capture_output <- function(expr) {
  capture.output(suppressWarnings(expr))
}

test_that("pokazStage: Basic message display", {
  output <- capture_output(pokazStage("Stage 1"))
  expect_equal(output, "* Stage 1")
})

test_that("pokazStage: Concatenates multiple arguments", {
  output <- capture_output(pokazStage("Processing", "step", "2"))
  expect_equal(output, "* Processing step 2")
})

test_that("pokazStage: Handles vector arguments", {
  output <- capture_output(pokazStage(c("Vector", "argument"), "stage"))
  expect_equal(output, "* Vector argument stage")
})

test_that("pokazStage: Writes to file correctly", {
  temp_file <- tempfile()
  pokazStage("File output stage.", file = temp_file, echo = FALSE)
  file_content <- readLines(temp_file)
  expect_equal(file_content, "* File output stage.")
  unlink(temp_file)
})

test_that("pokazStage: Handles both console and file output", {
  temp_file <- tempfile()
  output <- capture_output(pokazStage("Both output stage.", file = temp_file))
  expect_equal(output, "* Both output stage.")
  file_content <- readLines(temp_file)
  expect_equal(file_content, "* Both output stage.")
  unlink(temp_file)
})

test_that("pokazStage: Correctly handles empty input", {
  output <- capture_output(pokazStage())
  expect_equal(output, "* ")
})

test_that("pokazStage: Handles non-character arguments", {
  output <- capture_output(pokazStage(1, 2, 3))
  expect_equal(output, "* 1 2 3")
})