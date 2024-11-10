#!/bin/bash

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
Rscript -e "
if (requireNamespace('$PACKAGE_NAME', quietly = TRUE)){
    cat('[1/6] \033[34mOld installation of $PACKAGE_NAME has been found in R environment. Reinstalling...\n\033[0m')
    suppressMessages(remove.packages('$PACKAGE_NAME'))
} else {
    cat('[1/6] \033[34mPackage $PACKAGE_NAME has not been found in R environment. Installing...\n\033[0m')
}
"

# Updating symlinks for R package exports
echo -e "[2/6] \033[34mForce R symlinks update\033[0m"
rm -fr R/
mkdir R

find inst -type f -name "*.R" -exec sh -c '
  for file; do
    if grep -q "@export" "$file"; then
      ln -sf "$(realpath "$file")" "R/$(basename "$file")"
    fi
  done
' sh {} +
echo -e "[3/6] \033[34mSymlinks created\033[0m"

# Remove old manuals directory
rm -rf man/


# Scripts force linking to Conda env
echo -e "[4/6] \033[34mForce scripts symlinks update\033[0m"
ln -sf "$(realpath ./inst/analys.sh)" "$CONDA_PREFIX/bin/analys"
ln -sf "$(realpath ./inst/pannagram.sh)" "$CONDA_PREFIX/bin/pannagram"
ln -sf "$(realpath ./inst/simsearch.sh)" "$CONDA_PREFIX/bin/simsearch"