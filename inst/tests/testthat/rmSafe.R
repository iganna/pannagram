library(testthat)

test_that("rmSafe: Removes an existing variable", {
  x <<- 42  # Assign to the global environment
  expect_true(exists("x", envir = globalenv()))
  rmSafe(x)
  expect_false(exists("x", envir = globalenv()))
})

test_that("rmSafe: Doesn't throw an error for a non-existent variable", {
  expect_false(exists("y", envir = globalenv()))
  expect_error(rmSafe(y), NA) # Expect no error
})

test_that("rmSafe: Removes an existing global variable", {
  x <<- 42
  expect_true(exists("x", envir = globalenv()))
  rmSafe(x)
  expect_false(exists("x", envir = globalenv()))
})

test_that("rmSafe: Handles non-existent global variables without error", {
  expect_false(exists("y", envir = globalenv())) # Ensure y doesn't exist
  expect_error(rmSafe(y), NA) # Expect no error
})

# Clean up in case previous tests failed
if (exists("x", envir = globalenv())) rm(x, envir = globalenv())

