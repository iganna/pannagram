library(testthat)

test_file(system.file("tests/testthat/readFastaMy.R", package = "pannagram"))
test_file(system.file("tests/testthat/writeFastaMy.R", package = "pannagram"))
test_file(system.file("tests/testthat/seq2int.R", package = "pannagram"))
test_file(system.file("tests/testthat/nt2seq.R", package = "pannagram"))
test_file(system.file("tests/testthat/splitSeq.R", package = "pannagram"))
test_file(system.file("tests/testthat/aln2mx.R", package = "pannagram"))
test_file(system.file("tests/testthat/mx2aln.R", package = "pannagram"))
test_file(system.file("tests/testthat/mx2seq.R", package = "pannagram"))