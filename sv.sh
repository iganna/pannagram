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

source utils/utils_bash.sh


print_usage() {
    pokaz_help
    cat << EOF
Usage: ${0##*/} [-pref_global PREFIX] [-ref_pref REF_PREFIX] [-path_consensus PATH_CONSENSUS]
                [-gff] [-te] [-graph] [-h|--help]

This script provides various functionalities depending on the options provided.

Options:
    -pref_global PREFIX        Set a global preference or configuration.
    -ref_pref REF_PREFIX       Specify a reference preference or parameter.
    -path_consensus PATH_CONSENSUS
                               Define the path to the consensus data.
    -gff                       Run the script in GFF mode, processing data in Generic Feature Format.
    -te                        Enable processing in Transposable Element mode.
    -graph                     Activate graph mode for graphical data processing or visualization.
    -h, --help                 Display this help message and exit.

Examples:
    ${0##*/} -pref_global 'global_config' -ref_pref 'ref1' -path_consensus '/path/to/data'
    ${0##*/} -gff -path_consensus '/path/to/gff_data'

EOF
}


# ----------------------------------------------------------------------------
#             PARAMETERS
# ----------------------------------------------------------------------------


# Initialize variables to determine which scripts to run
run_gff=false
run_te=false
run_graph=false

# Parse command line arguments
while [ $# -gt 0 ]; do
    case $1 in
        -pref_global) pref_global=$2; shift 2;;
        -ref_pref) ref_pref=$2; shift 2;;
        -path_consensus) path_consensus=$2; shift 2;;
        -gff) run_gff=true; shift;;
        -te)    run_te=true; shift;;
        -graph)    run_graph=true; shift;;
        -h|--help) print_usage; exit 0;;
        *) print_usage; exit 1;;
    esac
done

check_missing_variable "ref_pref"

check_missing_variable "pref_global"
pref_global=$(add_symbol_if_missing "$pref_global" "/")

path_consensus="${path_consensus:-${pref_global}consensus/}"
path_consensus=$(add_symbol_if_missing "$path_consensus" "/")

# ----------------------------------------------------------------------------
#             MAIN
# ----------------------------------------------------------------------------

# Execute scripts based on the provided keys
if [ "$run_gff" = true ]; then

    # Philosophy: GFF does not make any sense without a pangenome consensus fasta. 
    # So, sonsensus should be run before GFF
    # Therefore, sequences of seSVs could also be produced together with GFFs.


    Rscript sv/sv_01_calling.R --path.cons ${path_consensus} --ref.pref  ${ref_pref}
fi

if [ "$run_te" = true ]; then
    Rscript sv/sv_02_te.R --path.cons ${path_consensus} --ref.pref  ${ref_pref}
fi

if [ "$run_graph" = true ]; then
    Rscript sv/sv_02_graph.R --path.cons ${path_consensus} --ref.pref  ${ref_pref}
fi
