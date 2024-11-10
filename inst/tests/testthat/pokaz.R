library(testthat)
library(crayon)

# Capture output for testing
capture_output <- function(expr) {
  capture.output(suppressWarnings(expr)) # Suppress crayon warnings for cleaner tests
}

test_that("pokaz: Basic message display", {
  output <- capture_output(pokaz("Hello, world!"))
  expect_equal(output, "  Hello, world!") # Check console output
})

test_that("pokaz: Concatenates multiple arguments", {
  output <- capture_output(pokaz("This", "is", "a", "test."))
  expect_equal(output, "  This is a test.")
})

test_that("pokaz: Handles vector arguments", {
  output <- capture_output(pokaz(c("Vector", "argument"), "test."))
  expect_equal(output, "  Vector argument test.")
})


test_that("pokaz: Writes to file correctly", {
  temp_file <- tempfile()
  pokaz("File output test.", file = temp_file, echo = FALSE) # echo = FALSE to avoid console output in test
  file_content <- readLines(temp_file)
  expect_equal(file_content, "  File output test.") # Check file content
  unlink(temp_file) # Remove the temporary file
})

test_that("pokaz: Handles both console and file output", {
  temp_file <- tempfile()
  output <- capture_output(pokaz("Both output test.", file = temp_file))
  expect_equal(output, "  Both output test.") # Check console output
  file_content <- readLines(temp_file)
  expect_equal(file_content, "  Both output test.") # Check file content
  unlink(temp_file)
})

test_that("pokaz: Correctly handles empty input", {
  output <- capture_output(pokaz())
  expect_equal(output, "  ")
})


test_that("pokaz: Handles non-character arguments", {
  output <- capture_output(pokaz(123, 4.56, TRUE))
  expect_equal(output, "  123 4.56 TRUE")
})