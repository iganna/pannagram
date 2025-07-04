#!/bin/bash

START_TIME=$(date +%s.%N)

CONDAENV_NAME="pannagram"
PACKAGE_NAME="pannagram"

if [ -f "pannagram_checks.sh" ]; then
  source pannagram_checks.sh "$CONDAENV_NAME" "$PACKAGE_NAME"
else
  echo -e "\n\033[31mError: 'pannagram_checks.sh' not found! Run '${0}' from the root of the repo please!\033[0m\n"
  exit 1
fi

echo -e "[4] \033[34mPannagram documentation installation\033[0m"
Rscript -e "
suppressMessages(devtools::document(quiet=TRUE))
"
echo -e "[5] \033[34mPannagram R package installation\033[0m"
Rscript -e "
suppressMessages(devtools::install(quiet=TRUE))
" && {
  END_TIME=$(date +%s.%N)
  ELAPSED_TIME=$(awk "BEGIN {printf \"%.1f\", $END_TIME - $START_TIME}")
  echo -e "\033[32mUser mode: Package $PACKAGE_NAME check and (re)installation process completed in $ELAPSED_TIME seconds.\033[0m"
}