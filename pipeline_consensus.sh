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
#             FUNCTIONS
# ----------------------------------------------------------------------------

source utils_bash.sh

# ----------------------------------------------------------------------------
#             MAIN
# ----------------------------------------------------------------------------


# ref_set=""
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
    -pref_global)
        # Get the value of the pref_global parameter
        pref_global="$2"
        additional_params+=" $1 $2"
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

if [[ -z "$ref_set" ]]; then
  echo "Error: -ref_set parameter is not set" >&2
  exit 1
fi



# Split the value of ref_set into separate words
IFS=',' read -ra refs_all <<< "$ref_set"

# Iterate over each word in ref_set
for ref0 in "${refs_all[@]}"; do
    # Call the second script pipeline.sh, passing the current word as an argument
    # ./pipeline.sh "$ref0" $additional_params

    # command="./pipeline.sh -ref_pref ${ref0} ${additional_params}"
    command="./work.sh -ref_pref ${ref0} ${additional_params}"
    echo "Executing command: ${command}"
    # eval "${command}"

done

# ----------------------------------------------------------------------------
#             COMMON CONSENSUS
# ----------------------------------------------------------------------------


check_missing_variable "pref_global"
pref_global=$(add_symbol_if_missing "$pref_global" "/")


path_consensus="${path_consensus:-${pref_global}consensus/}"
path_consensus=$(add_symbol_if_missing "$path_consensus" "/")

echo ${path_consensus}



cores=10

ref0=${refs_all[0]}
for ((i = 1; i < ${#refs_all[@]}; i++)); do
    ref1=${refs_all[i]}
    
    Rscript comb_ref_free.R --path.cons ${path_consensus} --ref0 ${ref0} --ref1 ${ref1} --cores ${cores}

done