#!/bin/bash

START_TIME=$(date +%s.%N)

CONDAENV_NAME="pannagram"
PACKAGE_NAME="pannagram"

source pannagram_checks.sh $CONDAENV_NAME $PACKAGE_NAME


# Full installation: documentation + R package
Rscript -e "
devtools::document()
devtools::install()
"

END_TIME=$(date +%s.%N)
ELAPSED_TIME=$(awk "BEGIN {printf \"%.1f\", $END_TIME - $START_TIME}")

# Format elapsed time in seconds
echo -e "\n\033[32mUser mode: Package $PACKAGE_NAME check and (re)installation process completed in $ELAPSED_TIME seconds.\033[0m\n"
