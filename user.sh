#!/bin/bash

START_TIME=$(date +%s.%N)

CONDAENV_NAME="pannagram"
PACKAGE_NAME="pannagram"

source pannagram_checks.sh $CONDAENV_NAME $PACKAGE_NAME


# Full installation: documentation + R package
echo -e "[5/6] \033[34mR package + documentation installation\033[0m"

Rscript -e "
suppressMessages(devtools::document(quiet=TRUE))
suppressMessages(devtools::install(quiet=TRUE))
"

END_TIME=$(date +%s.%N)
ELAPSED_TIME=$(awk "BEGIN {printf \"%.1f\", $END_TIME - $START_TIME}")

# Format elapsed time in seconds
echo -e "[6/6] \033[32mUser mode: Package $PACKAGE_NAME check and (re)installation process completed in $ELAPSED_TIME seconds.\033[0m"
