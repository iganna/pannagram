library(testthat)
source("../../utils/utils.R")

create_temp_fasta <- function() {
  tempfile(fileext = ".fasta")
}

# Вспомогательная функция для чтения содержимого файла
read_fasta_content <- function(file) {
  readLines(file)
}

test_that("writeFastaMy: Couldn't write sequences correctly", {
  temp_file <- create_temp_fasta()
  sequences <- c("ATGC", "GCTA")
  names(sequences) <- c("seq1", "seq2")
  
  writeFastaMy(sequences, temp_file)
  
  content <- read_fasta_content(temp_file)
  expect_equal(content, c(">seq1", "ATGC", ">seq2", "GCTA"))
  
  unlink(temp_file)
})

test_that("writeFastaMy: Couldn't handle sequences without names", {
  temp_file <- create_temp_fasta()
  sequences <- c("ATGC", "GCTA")
  
  writeFastaMy(sequences, temp_file)
  
  content <- read_fasta_content(temp_file)
  expect_equal(content, c(">seq_1", "ATGC", ">seq_2", "GCTA"))
  
  unlink(temp_file)
})

test_that("writeFastaMy: Couldn't use custom prefix", {
  temp_file <- create_temp_fasta()
  sequences <- c("ATGC", "GCTA")
  
  writeFastaMy(sequences, temp_file, pref = "custom_")
  
  content <- read_fasta_content(temp_file)
  expect_equal(content, c(">custom_1", "ATGC", ">custom_2", "GCTA"))
  
  unlink(temp_file)
})

test_that("writeFastaMy: Couldn't append sequences to existing file", {
  temp_file <- create_temp_fasta()
  sequences1 <- c("ATGC")
  sequences2 <- c("GCTA")
  names(sequences1) <- "seq1"
  names(sequences2) <- "seq2"
  
  writeFastaMy(sequences1, temp_file)
  writeFastaMy(sequences2, temp_file, append = TRUE)
  
  content <- read_fasta_content(temp_file)
  expect_equal(content, c(">seq1", "ATGC", ">seq2", "GCTA"))
  
  unlink(temp_file)
})

test_that("writeFastaMy: Couldn't overwrite existing file when append is FALSE", {
  temp_file <- create_temp_fasta()
  sequences1 <- c("ATGC")
  sequences2 <- c("GCTA")
  names(sequences1) <- "seq1"
  names(sequences2) <- "seq2"
  
  writeFastaMy(sequences1, temp_file)
  writeFastaMy(sequences2, temp_file, append = FALSE)
  
  content <- read_fasta_content(temp_file)
  expect_equal(content, c(">seq2", "GCTA"))
  
  unlink(temp_file)
})

test_that("writeFastaMy: Couldn't handle empty sequences vector", {
  temp_file <- create_temp_fasta()
  sequences <- character(0)
  
  writeFastaMy(sequences, temp_file)
  
  expect_true(file.exists(temp_file))
  content <- read_fasta_content(temp_file)
  expect_equal(content, character(0))
  
  expect_equal(file.size(temp_file), 0)
  
  unlink(temp_file)
})

test_that("writeFastaMy: Couldn't handle empty sequences vector with append = TRUE", {
  temp_file <- create_temp_fasta()
  
  initial_sequences <- c("ATGC")
  names(initial_sequences) <- "seq1"
  writeFastaMy(initial_sequences, temp_file)
  
  empty_sequences <- character(0)
  writeFastaMy(empty_sequences, temp_file, append = TRUE)
  
  content <- read_fasta_content(temp_file)
  expect_equal(content, c(">seq1", "ATGC"))
  
  unlink(temp_file)
})

test_that("writeFastaMy: Couldn't use provided seq.names", {
  temp_file <- create_temp_fasta()
  sequences <- c("ATGC", "GCTA")
  seq.names <- c("custom1", "custom2")
  
  writeFastaMy(sequences, temp_file, seq.names = seq.names)
  
  content <- read_fasta_content(temp_file)
  expect_equal(content, c(">custom1", "ATGC", ">custom2", "GCTA"))
  
  unlink(temp_file)
})

test_that("writeFastaMy: Couldn't handle long sequences", {
  temp_file <- create_temp_fasta()
  long_seq <- paste(rep("ATGC", 100), collapse = "")
  sequences <- c(long_seq)
  names(sequences) <- "long_seq"
  
  writeFastaMy(sequences, temp_file)
  
  content <- read_fasta_content(temp_file)
  expect_equal(content, c(">long_seq", long_seq))
  
  unlink(temp_file)
})

test_that("writeFastaMy: Couldn't handle special characters in sequence names", {
  temp_file <- create_temp_fasta()
  sequences <- c("ATGC", "GCTA", "TCGC")
  names(sequences) <- c("seq with spaces", "seq_with_underscore", "seq:with|(other/ special%chars^#-=_)")
  
  writeFastaMy(sequences, temp_file)
  
  content <- read_fasta_content(temp_file)
  expect_equal(content, c(">seq with spaces", "ATGC", ">seq_with_underscore", "GCTA", ">seq:with|(other/ special%chars^#-=_)", "TCGC"))
  
  unlink(temp_file)
})

test_that("writeFastaMy: Couldn't handle non-standard nucleotides", {
  temp_file <- create_temp_fasta()
  sequences <- c("ATG-CN", "RYM*KWS")
  names(sequences) <- c("seq1", "seq2")
  
  writeFastaMy(sequences, temp_file)
  
  content <- read_fasta_content(temp_file)
  expect_equal(content, c(">seq1", "ATG-CN", ">seq2", "RYM*KWS"))
  
  unlink(temp_file)
})

test_that("writeFastaMy: Couldn't handle sequences of different lengths", {
  temp_file <- create_temp_fasta()
  sequences <- c("ATGC", "GCTAGCT", "A")
  names(sequences) <- c("seq1", "seq2", "seq3")
  
  writeFastaMy(sequences, temp_file)
  
  content <- read_fasta_content(temp_file)
  expect_equal(content, c(">seq1", "ATGC", ">seq2", "GCTAGCT", ">seq3", "A"))
  
  unlink(temp_file)
})

test_that("writeFastaMy: Couldn't handle file paths with spaces", {
  temp_dir <- tempdir()
  temp_file <- file.path(temp_dir, "test file.fasta")
  sequences <- c("ATGC", "GCTA")
  names(sequences) <- c("seq1", "seq2")
  
  writeFastaMy(sequences, temp_file)
  
  expect_true(file.exists(temp_file))
  content <- read_fasta_content(temp_file)
  expect_equal(content, c(">seq1", "ATGC", ">seq2", "GCTA"))
  
  unlink(temp_file)
})