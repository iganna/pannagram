#!/bin/bash

# Exit if script is run directly (not sourced)
[[ "${BASH_SOURCE[0]}" == "${0}" ]] && {
    echo -e "\n\033[31mDo not run '${0}' directly!\033[0m"
    echo -e "Please use: 'developer.sh' or 'user.sh' instead.\n"
    exit 1
}

CONDAENV_NAME=$1 # "pannagram"
PACKAGE_NAME=$2  # "pannagram"

# Conda environment check
if [[ "$CONDA_DEFAULT_ENV" != "$CONDAENV_NAME" ]]; then
  echo -e "\n\033[31mThis script must be run inside the '$CONDAENV_NAME' conda environment. Use '$CONDAENV_NAME.yaml' to recreate it.\033[0m\n"
  exit 1
fi

# Check if R is within the Conda environment
R_PATH=$(which R)
if [[ "$R_PATH" != "$CONDA_PREFIX/bin/R" ]]; then
  echo -e "\n\033[31mPath to your R interpreter leads to '$R_PATH' which is not inside current Conda env (has to be '$CONDA_PREFIX/bin/R'). Use original '$CONDAENV_NAME.yaml' to recreate the proper env. \033[0m\n"
  exit 1
fi

# Remove old installation of pannagram
echo -e "[1] \033[34mRemoving any potential old installations of $PACKAGE_NAME\033[0m"
Rscript -e "
invisible(suppressMessages(tryCatch(remove.packages('$PACKAGE_NAME'), error = function(e) NULL)))
"
rm -fr $CONDA_PREFIX/lib/R/library/$PACKAGE_NAME # double take


# Updating symlinks for R package exports
echo -e "[2] \033[34mUpdating symlinks for exported R functions\033[0m"
rm -fr R/
mkdir R

find inst -type f -name "*.R" -exec sh -c '
  for file; do
    if grep -q "@export" "$file"; then
      ln -sf "$(realpath "$file")" "R/$(basename "$file")"
    fi
  done
' sh {} +

# Remove old manuals directory
rm -rf man/


# Scripts force linking to Conda env
echo -e "[3] \033[34mUpdating symlinks for Bash scripts\033[0m"
ln -sf "$(realpath ./inst/analys.sh)" "$CONDA_PREFIX/bin/analys"
ln -sf "$(realpath ./inst/pannagram.sh)" "$CONDA_PREFIX/bin/pannagram"
ln -sf "$(realpath ./inst/simsearch.sh)" "$CONDA_PREFIX/bin/simsearch"
ln -sf "$(realpath ./inst/chromotools.sh)" "$CONDA_PREFIX/bin/chromotools"

