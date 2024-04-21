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
    echo "Usage: $0 [-blocks] [-seq] [-sv]"
    echo "  -blocks    Run analys_01_blocks.R script"
    echo "  -seq       Run analys_02_seq.R script"
    # echo "  -sv        Run analys_03_sv.R script"
    echo "  -h, --help Display this help message"
}

# ----------------------------------------------------------------------------
#             PARAMETERS
# ----------------------------------------------------------------------------


# Initialize variables to determine which scripts to run
run_blocks=false
run_seq=false
run_aln=false
run_snp=false
aln_type='msa_'
# run_sv=false

# Parse command line arguments
while [ $# -gt 0 ]; do
    case $1 in
        -path_msa) pref_global=$2; shift 2;;
        -ref) ref_pref=$2; shift 2;;
        -aln_type) aln_type=$2; shift 2;;
        -path_consensus) path_consensus=$2; shift 2;;
        -path_chromosomes) path_chromosomes=$2; shift 2 ;;
        -cores) cores=$2; shift 2 ;;
        -blocks) run_blocks=true; shift;;
        -seq)    run_seq=true; shift;;
        -aln)    run_aln=true; shift;;
        -snp)    run_snp=true; shift;;
        # -sv)     run_sv=true; shift;;
        -h|--help) print_usage; exit 0;;
        *) print_usage; exit 1;;
    esac
done

cores="${cores:-1}"  # Number of cores
pokaz_message "Number of cores: ${cores}"

check_missing_variable "ref_pref"

check_missing_variable "pref_global"
pref_global=$(add_symbol_if_missing "$pref_global" "/")

path_consensus="${path_consensus:-${pref_global}consensus/}"
path_consensus=$(add_symbol_if_missing "$path_consensus" "/")

path_chromosomes="${path_chromosomes:-${pref_global}chromosomes/}"
path_chromosomes=$(add_symbol_if_missing "$path_chromosomes" "/")

# ----------------------------------------------------------------------------
#             MAIN
# ----------------------------------------------------------------------------

# Execute scripts based on the provided keys
if [ "$run_blocks" = true ]; then

    Rscript analys/analys_01_blocks.R --path.cons ${path_consensus} --ref.pref  ${ref_pref} --cores ${cores} --aln.type ${aln_type}
fi

if [ "$run_seq" = true ]; then

    Rscript analys/analys_02_seq_cons.R --path.cons ${path_consensus} --ref.pref  ${ref_pref} --path.chromosomes ${path_chromosomes}  --aln.type ${aln_type} --cores ${cores}
fi

if [ "$run_aln" = true ]; then

    Rscript analys/analys_03_seq_aln.R --path.cons ${path_consensus} --ref.pref  ${ref_pref} --path.chromosomes ${path_chromosomes} --aln.type ${aln_type} --cores ${cores}
fi


if [ "$run_snp" = true ]; then

    Rscript analys/analys_04_snp.R --path.cons ${path_consensus} --ref.pref  ${ref_pref} --path.chromosomes ${path_chromosomes}  --aln.type ${aln_type} --cores ${cores}
fi


# if [ "$run_sv" = true ]; then
#     Rscript analys_03_sv.R
# fi