#!/usr/bin/env bash
# ===============================================
# Мalidation script for Pannagram installation
# ===============================================

set -e  # Exit immediately if any command fails

# Check required CLI tools
for tool in pannagram simsearch features chromotools Rscript; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: module $tool not found."
        exit 1
    fi
done

# Check Rscript
if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript not found."
    exit 1
fi

# Check pannagram R package
if ! Rscript -e "if (!requireNamespace('pannagram', quietly = TRUE)) quit(status = 1)" >/dev/null 2>&1; then
    echo "ERROR: pannagram R package not installed."
    exit 1
fi

# Test pannagram function
test_command="library(pannagram); pannagram::pokaz('Pannagram has been installed successfully')"
if ! Rscript -e "${test_command}" >/dev/null 2>&1; then
    echo "ERROR: pannagram test command failed."
    exit 1
fi

echo "Pannagram has been installed successfully!"
