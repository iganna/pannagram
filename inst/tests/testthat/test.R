library(testthat)

utils_file <- system.file("utils/utils.R", package = "pannagram")
source(utils_file)
all_functions <- ls()

test_files <- list.files(system.file("tests/testthat", package = "pannagram"), pattern = "\\.R$")
test_functions <- gsub("\\.R$", "", test_files)
test_functions <- setdiff(test_functions, "test")

for (test_file in test_files) {
  if (test_file != "test.R") {
    test_file(system.file(paste0("tests/testthat/", test_file), package = "pannagram"))
  }
}

untested_functions <- setdiff(all_functions, test_functions)
n_untested <- length(untested_functions)

if (n_untested > 0) {
  cat(paste0("\nWarning: The following ", n_untested, " functions are not covered by tests:\n"))
  cat(paste("- ", untested_functions, collapse = "\n"))
  cat("\n")
} else {
  cat("\nAll functions are covered by tests.\n")
}