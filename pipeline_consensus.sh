# This escript performs all stages of the alignment of genomes, 
# when the reference genome is already identified

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
#             USAGE
# ----------------------------------------------------------------------------

# ./pipeline.sh -pref_global 'rhiz' -ref_pref 'ref_1021' -n_chr_ref 1 -path_in '../rhizobia/' -n_chr_query 1 -sort_chr_len T 

# ./pipeline.sh -pref_global 'ly' -ref_pref '0' -path_chr_ref "../pb_chromosomes/" -n_chr_ref 5 -path_in '../lyrata/' -n_chr_query 8


# ./pipeline.sh -pref_global 'ly2' -ref_pref '0' -path_chr_ref "../pb_chromosomes/" -n_chr_ref 5 -path_in '../lyrata/' -n_chr_query 8

# ./pipeline.sh -pref_global 'toy' -ref_pref '0'  -n_chr_ref 5 -path_in '../pb_genomes/' -n_chr_query 5 -all_cmp F -acc_anal 'acc_analysis.txt'


# ----------------------------------------------------------------------------
#             MAIN
# ----------------------------------------------------------------------------


ref_set=""
additional_params=""

# Iterate through all arguments
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -ref_set)
        # Get the value of the ref_set parameter
        ref_set="$2"
        shift # Skip the value of the parameter
        shift # Skip the parameter itself
        ;;
    *)
        # Add the remaining arguments as additional parameters
        additional_params+=" $1"
        shift
        ;;
  esac
done

# Split the value of ref_set into separate words
IFS=',' read -ra words <<< "$ref_set"


# Iterate over each word in ref_set
for word in "${words[@]}"; do
    # Call the second script pipeline.sh, passing the current word as an argument
    # ./pipeline.sh "$word" $additional_params

    # command="./pipeline.sh -ref_pref ${word} ${additional_params}"
    command="./work.sh -ref_pref ${word} ${additional_params}"
    echo "Executing command: ${command}"
    eval "${command}"

done

# ----------------------------------------------------------------------------
#             COMMON CONSENSUS
# ----------------------------------------------------------------------------

path_consensus="${path_consensus:-${pref_global}consensus/}"
path_consensus=$(add_symbol_if_missing "$path_consensus" "/")