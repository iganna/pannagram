library(testthat)

test_that("aln2pos: Converts aligned sequences to position matrix correctly", {
  seqs <- c(seq1 = "AT-GC", seq2 = "A--GC")
  expected_matrix <- matrix(c(1, 2, 0, 3, 4, 1, 0, 0, 3, 4), nrow = 2, byrow = TRUE)
  rownames(expected_matrix) <- c("seq1", "seq2")
  result <- aln2pos(seqs)
  expect_equal(result, expected_matrix)
})

test_that("aln2pos: Handles empty sequences", {
  seqs <- c(seq1 = "", seq2 = "")
  expected_matrix <- matrix(numeric(0), nrow = 2, ncol = 0)
  rownames(expected_matrix) <- c("seq1", "seq2")
  result <- aln2pos(seqs)
  expect_equal(result, expected_matrix)
})


test_that("aln2pos: Handles sequences with only gaps", {
  seqs <- c(seq1 = "---", seq2 = "---")
  expected_matrix <- matrix(0, nrow = 2, ncol = 3)
  rownames(expected_matrix) <- c("seq1", "seq2")
  result <- aln2pos(seqs)
  expect_equal(result, expected_matrix)
})

test_that("aln2pos: Throws error for unequal sequence lengths", {
  seqs <- c(seq1 = "AT-GC", seq2 = "A--G")
  expect_error(aln2pos(seqs), "All sequences must be of the same length.")
})

test_that("aln2pos: Works with unnamed sequences", {
  seqs <- c("AT-GC", "A--GC")
  expected_matrix <- matrix(c(1, 2, 0, 4, 3, 1, 0, 0, 4, 3), nrow = 2, byrow = TRUE)
  result <- aln2pos(seqs)
  expect_equal(result, expected_matrix)
  expect_null(rownames(result)) # Expect no row names when input is unnamed
})

test_that("aln2pos: Works with a single sequence", {
    seqs <- c(seq1 = "AT-GC")
    expected_matrix <- matrix(c(1, 2, 0, 4, 3), nrow = 1, byrow = TRUE)
    rownames(expected_matrix) <- c("seq1")
    result <- aln2pos(seqs)
    expect_equal(result, expected_matrix)
})

test_that("aln2pos: Handles sequences with internal gaps followed by characters", {
    seqs <- c(seq1 = "A-T-G-", seq2 = "AT--G-")
    expected_matrix <- matrix(c(1, 0, 2, 0, 3, 0, 1, 2, 0, 0, 3, 0), nrow = 2, byrow = TRUE)
    rownames(expected_matrix) <- c("seq1", "seq2")
    result <- aln2pos(seqs)
    expect_equal(result, expected_matrix)
})