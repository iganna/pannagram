#!/bin/bash

# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

# Exit immediately if any command returns a non-zero status
set -e

# Keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# Define a trap for the EXIT signal
trap 'catch $?' EXIT

# Function to handle the exit signal
catch() {
    # Check if the exit code is non-zero
    if [ $1 -ne 0 ]; then
        echo "\"${last_command}\" command failed with exit code $1."
    fi
}

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source utils_bash.sh


print_usage() {
    echo "Usage: $0 [-blocks] [-seq] [-sv]"
    echo "  -blocks    Run analys_01_blocks.R script"
    echo "  -seq       Run analys_02_seq.R script"
    echo "  -sv        Run analys_03_sv.R script"
    echo "  -h, --help Display this help message"
}

# ----------------------------------------------------------------------------
#             PARAMETERS
# ----------------------------------------------------------------------------


# Initialize variables to determine which scripts to run
run_blocks=false
run_seq=false
run_sv=false

# Parse command line arguments
for arg in "$@"; do
    case $arg in
        -blocks) run_blocks=true ;;
        -seq)    run_seq=true ;;
        -sv)     run_sv=true ;;
        -h|--help) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done


# ----------------------------------------------------------------------------
#             MAIN
# ----------------------------------------------------------------------------

# Execute scripts based on the provided keys
if [ "$run_blocks" = true ]; then
    Rscript analys_01_blocks.R
fi

if [ "$run_seq" = true ]; then
    Rscript analys_02_seq.R
fi

if [ "$run_sv" = true ]; then
    Rscript analys_03_sv.R
fi