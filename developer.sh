#!/bin/bash

CONDAENV_NAME="pannagram"
PACKAGE_NAME="pannagram"

# Conda env chack
if [[ "$CONDA_DEFAULT_ENV" != "$CONDAENV_NAME" ]]; then
  echo -e "\n\033[31mThis script must be run inside the '$CONDAENV_NAME' conda environment.\033[0m\n"
  exit 1
fi

# Messing with R env
Rscript -e "
if (requireNamespace('$PACKAGE_NAME', quietly = TRUE)){
    cat('\n\033[34mPackage $PACKAGE_NAME has been found in R environment.\nReinstalling...\n\033[0m\n')
    remove.packages('$PACKAGE_NAME')
} else {
    cat('\n\033[34mPackage $PACKAGE_NAME has not been found in R environment.\nInstalling...\033[0m\n')
}
"

# Force updating for symlinks
echo -e "\n\033[34mForce symlinks update\033[0m\n"
rm -fr R/
mkdir R

find inst -type f -name "*.R" -exec sh -c '
  for file; do
    if grep -q "@export" "$file"; then
      ln -sf "$(realpath "$file")" "R/$(basename "$file")"
    fi
  done
' sh {} +
echo -e "\n\033[34mSymlinks created\033[0m\n"

rm -rf man/

Rscript -e "
devtools::document()
devtools::install()
"

echo -e "\n\033[32mPackage $PACKAGE_NAME check and installation process completed.\033[0m\n"