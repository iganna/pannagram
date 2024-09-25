# Testing Pannagram

Currently only some functions from `utils/utils.R` are covered with tests. Here's how to get names of uncovered ones:
```bash
UTILS_FILE="inst/utils/utils.R"
TESTS_DIR="inst/tests/testthat/"

function_names=$(grep -oP '^\s*([a-zA-Z0-9_]+)\s*<-\s*function' "$UTILS_FILE" | awk '{print $1}')
test_files=$(ls "$TESTS_DIR" | grep -oP '^\s*([a-zA-Z0-9_]+)\.R' | sed 's/\.R$//')

# prints uncovered functions:
comm -23 <(grep -oP '^\s*([a-zA-Z0-9_]+)\s*<-\s*function' "$UTILS_FILE" | awk '{print $1}' | sort) <(ls "$TESTS_DIR" | grep -oP '^\s*([a-zA-Z0-9_]+)\.R' | sed 's/\.R$//' | sort)
```

1. Activate the pannagram environment from `pannagram_dev.yaml`
2. Run tests (from the root of the repo): `Rscript inst/tests/testthat/test.R`
