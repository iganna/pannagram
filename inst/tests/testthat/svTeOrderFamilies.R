library(testthat)
source(system.file("utils/utils.R", package = "pannagram"))

# the copy in the repo wins over the installed one: it may not be rebuilt yet
path.func <- Filter(file.exists, c("../../analys/sv_te_func.R", "inst/analys/sv_te_func.R"))
path.func <- if (length(path.func)) path.func[1] else system.file("analys/sv_te_func.R", package = "pannagram")
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

makeFamilies <- function() {
  ltr <- rndSeq(300, seed = 11)
  el.ltr <- function(seed) paste0(mutSeq(ltr, 0.02, seed = seed),
                                  rndSeq(1500, seed = 12),
                                  mutSeq(ltr, 0.03, seed = seed + 1))
  tir <- rndSeq(24, seed = 21)
  el.tir <- function(seed) paste0(tir, mutSeq(rndSeq(1200, seed = 22), 0.03, seed = seed), revCompl(tir))

  s.ltr.1 <- el.ltr(101)
  s.ltr.2 <- el.ltr(201)

  seqs <- c(sv_ltr_1  = s.ltr.1,
            sv_ltr_2  = s.ltr.2,
            sv_ltr_cut = substr(s.ltr.1, 1, 900),   # a truncated copy, has no termini
            sv_tir_1  = el.tir(301),
            sv_tir_2  = el.tir(401))

  families <- c(sv_ltr_1 = "fam_1", sv_ltr_2 = "fam_1", sv_ltr_cut = "fam_1",
                sv_tir_1 = "fam_2", sv_tir_2 = "fam_2")

  list(seqs = seqs, families = families)
}

test_that("svTeOrderFamilies: classifies families and skips truncated copies", {
  d <- makeFamilies()
  res <- svTeOrderFamilies(d$seqs, d$families, echo = FALSE)

  expect_named(res, c("families", "repr"))
  expect_equal(nrow(res$families), 2)

  fam <- res$families[order(res$families$family), ]
  expect_equal(fam$te.class, c("LTR", "TIR"))
  expect_equal(fam$n.members, c(3, 2))
  expect_equal(fam$n.repr, c(2, 2))          # the truncated copy is not a representative
  expect_equal(fam$support, c(1, 1))
  expect_equal(nrow(res$repr), 4)
})

test_that("svTeOrderFamilies: accepts a data frame of families", {
  d <- makeFamilies()
  df <- data.frame(sv = names(d$families), family = as.character(d$families),
                   stringsAsFactors = FALSE)

  res <- svTeOrderFamilies(d$seqs, df, echo = FALSE)
  expect_equal(nrow(res$families), 2)
  expect_equal(sort(res$families$te.class), c("LTR", "TIR"))
})

test_that("svTeOrderFamilies: a family without structure is unknown", {
  seqs <- c(sv_1 = rndSeq(2000, seed = 71), sv_2 = rndSeq(2000, seed = 81))
  families <- c(sv_1 = "fam_x", sv_2 = "fam_x")

  res <- svTeOrderFamilies(seqs, families, echo = FALSE)
  expect_equal(res$families$te.class, "unknown")
  expect_equal(res$families$support, 0)
})

test_that("svTeOrderFamilies: SVs without sequences are skipped", {
  d <- makeFamilies()
  families <- c(d$families, sv_absent = "fam_3")

  capture.output(res <- svTeOrderFamilies(d$seqs, families, echo = FALSE))
  expect_equal(nrow(res$families), 2)
  expect_false("fam_3" %in% res$families$family)
})
