#!/bin/bash

# This script contains a pipeline to rearrange one genome in the order of another, 
# solely by rearranging parts between chromosomes, 
# not by rearranging fragments within chromosomes

# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

source utils/error_block.sh

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source utils/utils_bash.sh


print_usage() {
    pokaz_help
#     cat << EOF
# Usage: ${0##*/} [-pref_global PREFIX] [-ref_pref REF_PREFIX] [-path_consensus PATH_CONSENSUS]
                

# Options:
   
# Examples:


# EOF
}


# ----------------------------------------------------------------------------
#             PARAMETERS
# ----------------------------------------------------------------------------


# Initialize variables to determine which scripts to run


# Parse command line arguments
while [ $# -gt 0 ]; do
    case $1 in
        -pref_global) pref_global=$2; shift 2;;
        -ref_pref) ref_pref=$2; shift 2;;
        -path_consensus) path_consensus=$2; shift 2;;

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


# search major
# search boundaries
# compare with the guide
# rearrange