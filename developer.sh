#!/bin/bash

START_TIME=$(date +%s.%N)

CONDAENV_NAME="pannagram"
PACKAGE_NAME="pannagram"

source pannagram_checks.sh $CONDAENV_NAME $PACKAGE_NAME

# Quick installation: no documentation
Rscript -e "
devtools::install(
  quick=TRUE,
  force = FALSE,
  upgrade = 'never',
  build_vignettes = FALSE,
  dependencies=FALSE
)
"

END_TIME=$(date +%s.%N)
ELAPSED_TIME=$(awk "BEGIN {printf \"%.1f\", $END_TIME - $START_TIME}")

echo -e "\n\033[32mDeveloper mode: Package $PACKAGE_NAME check and (re)installation process completed in $ELAPSED_TIME seconds.\033[0m\n"