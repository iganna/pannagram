#!/bin/bash

START_TIME=$(date +%s)

CONDAENV_NAME="pannagram"
PACKAGE_NAME="pannagram"

if [ -f "pannagram_checks.sh" ]; then
  source pannagram_checks.sh "$CONDAENV_NAME" "$PACKAGE_NAME"
else
  echo -e "\n\033[31mError: 'pannagram_checks.sh' not found! Run '${0}' from the root of the repo please!\033[0m\n"
  exit 1
fi

# Run a command quietly; print its full output only if it fails.
run_quiet() {
  local log; log=$(mktemp)
  if ! "$@" >"$log" 2>&1; then
    echo -e "\033[31mFailed — full output:\033[0m"
    cat "$log"; rm -f "$log"; exit 1
  fi
  rm -f "$log"
}

echo -e "[4] \033[34mPannagram documentation installation\033[0m"
rm -rf man/
run_quiet Rscript -e "suppressMessages(devtools::document(quiet=TRUE))"

echo -e "[5] \033[34mPannagram R package installation\033[0m"
# Use plain 'R CMD INSTALL' instead of devtools::install(): dependencies are supplied
# by the conda env, and devtools' pak backend is broken on some macOS/arm builds
# (platform string 'darwinXX.0.0' is misread as unsupported).
run_quiet R CMD INSTALL --no-multiarch --with-keep.source .

END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
echo -e "\033[32mUser mode: Package $PACKAGE_NAME check and (re)installation process completed in $ELAPSED_TIME seconds.\033[0m"
