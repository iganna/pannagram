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
    echo "  -h, --help Display this help message"
}

# ----------------------------------------------------------------------------
#             PARAMETERS
# ----------------------------------------------------------------------------


# Initialize variables to determine which scripts to run
run_blocks=false
run_seq=false

# Parse command line arguments
while [ $# -gt 0 ]; do
    case $1 in
        -pref_global) pref_global=$2; shift 2;;
        -ref_pref) ref_pref=$2; shift 2;;
        -path_consensus) path_consensus=$2; shift 2;;
        -blocks) run_blocks=true; shift;;
        -seq)    run_seq=true; shift;;
        -h|--help) print_usage; exit 0;;
        *) print_usage; exit 1;;
    esac
done

path_consensus="${path_consensus:-${pref_global}consensus/}"
path_consensus=$(add_symbol_if_missing "$path_consensus" "/")

# ----------------------------------------------------------------------------
#             MAIN
# ----------------------------------------------------------------------------

# Execute scripts based on the provided keys
if [ "$run_blocks" = true ]; then
    Rscript analys/analys_01_blocks.R --path.cons ${path_consensus} --ref.pref  ${ref_pref}
fi

if [ "$run_seq" = true ]; then
    Rscript analys/analys_02_seq.R
fi
