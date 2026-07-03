#!/bin/bash

# Dynamically creating R directory upon building
# And populating it with links to R scripts with `@export`-s
# Now documentation will be correctly generated
mkdir R
find inst -type f -name "*.R" -exec sh -c '
  for file; do
    if grep -q "@export" "$file"; then
      ln -sf "$(realpath "$file")" "R/$(basename "$file")"
    fi
  done
' sh {} +

Rscript -e 'devtools::document()'
Rscript -e 'devtools::install()'

# Linking executable scripts to PATH
ln -sf $PREFIX/lib/R/library/pannagram/features.sh $PREFIX/bin/features
ln -sf $PREFIX/lib/R/library/pannagram/pannagram.sh $PREFIX/bin/pannagram
ln -sf $PREFIX/lib/R/library/pannagram/simsearch.sh $PREFIX/bin/simsearch
ln -sf $PREFIX/lib/R/library/pannagram/chromotools.sh $PREFIX/bin/chromotools
