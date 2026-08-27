library(testthat)
source(system.file("utils/utils.R", package = "pannagram"))

# the copy in the repo wins over the installed one: it may not be rebuilt yet
path.func <- Filter(file.exists, c("../../analys/graph_func.R", "inst/analys/graph_func.R"))
path.func <- if (length(path.func)) path.func[1] else system.file("analys/graph_func.R", package = "pannagram")
source(path.func)

# Deterministic sequences, independent of the RNG of R
lcg <- function(n, seed = 1) {
  x <- seed
  v <- numeric(n)
  for (i in seq_len(n)) {
    x <- (69069 * x + 1) %% 2147483648
    v[i] <- x / 2147483648
  }
  v
}
rndSeq <- function(n, seed = 1) paste0(c("A", "C", "G", "T")[floor(lcg(n, seed) * 4) + 1], collapse = "")
mutSeq <- function(s, p, seed = 1) {
  v <- seq2nt(s)
  u <- lcg(2 * length(v), seed)
  idx <- which(u[seq_along(v)] < p)
  v[idx] <- c("A", "C", "G", "T")[floor(u[length(v) + idx] * 4) + 1]
  nt2seq(v)
}

suppressMessages(library(pannagram))   # orfFinder, revComplSeq, parseStrings

# names of SVs always carry the length after "|"
nm <- function(x, s) paste0(x, "|", nchar(s))

test_that("getComponentSequences: brings the members to one strand", {
  s <- rndSeq(600, seed = 11)
  seqs <- c(s, revComplSeq(s))
  names(seqs) <- c(nm("a", s), nm("b", s))

  nst <- data.frame(name.query = nm("a", s), name.target = nm("b", s), strand = "-",
                    coverage.query = 1, coverage.target = 1, stringsAsFactors = FALSE)

  res <- getComponentSequences(names(seqs), seqs, nst)
  expect_length(res, 2)
  expect_equal(toupper(res[[1]]), toupper(res[[2]]))
})

test_that("getComponentSequences: a family without ORFs keeps the seed as is", {
  # stop codons in every frame on both strands: no ORF at all
  s <- paste0(rep("TAATAGTGA", 8), collapse = "")
  expect_null(orfFinder(s)$pos)

  seqs <- c(s, revComplSeq(s))
  names(seqs) <- c(nm("a", s), nm("b", s))
  nst <- data.frame(name.query = nm("a", s), name.target = nm("b", s), strand = "-",
                    coverage.query = 1, coverage.target = 1, stringsAsFactors = FALSE)

  res <- getComponentSequences(names(seqs), seqs, nst)
  expect_length(res, 2)
  expect_equal(toupper(res[[nm("a", s)]]), toupper(s))   # the seed is taken as "+"
  expect_equal(toupper(res[[1]]), toupper(res[[2]]))
})

test_that("getComponentSequences: the orientation is propagated along the edges", {
  s <- rndSeq(600, seed = 21)
  seqs <- c(s, revComplSeq(s), s)
  names(seqs) <- c(nm("a", s), nm("b", s), nm("c", s))

  # "c" is connected to "b" only, not to the seed "a"
  nst <- data.frame(name.query = c(nm("a", s), nm("b", s)),
                    name.target = c(nm("b", s), nm("c", s)),
                    strand = c("-", "-"), coverage.query = 1, coverage.target = 1,
                    stringsAsFactors = FALSE)

  res <- getComponentSequences(names(seqs), seqs, nst)
  expect_length(res, 3)
  expect_equal(toupper(res[[1]]), toupper(res[[2]]))
  expect_equal(toupper(res[[1]]), toupper(res[[3]]))
})
