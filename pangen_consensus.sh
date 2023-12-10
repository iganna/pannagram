# This escript performs all stages of the alignment of genomes, 
# when the reference genome is already identified


# ./pipeline_consensus.sh -pref_global '../pan_test/anna_ly' -ref_set 'NT1_220222,TE11_final' -n_chr_ref 8 -path_in '../lyrata_updated/' -n_chr_query 8 -cores 20  -all_cmp F


# ./pipeline_consensus.sh -pref_global '../pan_test/anna_aln' -ref_set 'NT1_220222,TE11_final' -n_chr_ref 8 -path_in '../lyrata_updated/' -n_chr_query 8 -cores 20  -all_cmp F


# ./pipeline_consensus.sh -pref_global '../pan_test/tom' -ref_set '0,6046-v1.1,6191-v1.1' -n_chr_ref 5 -path_in '../pb_updated/' -n_chr_query 5 -cores 30 -acc_anal acc_tom.txt -all_cmp F

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
Usage: ${0##*/} [OPTIONS]

This script performs specific tasks based on various input parameters. 

Options:
    -ref_set REF_SET        Specify the reference set to be used.
    -pref_global PREFIX     Global prefix for the analysis.
    -cores CORES            Number of cores for processing.
    -path_consensus PATH    Path to the consensus directory.
    -path_chr_acc PATH      Path to chromosome accession files.

Additional parameters can be provided as needed. These will be passed on
to subsequent processes or commands used within the script.

Examples:
    ${0##*/} -pref_global 'output_folder' -ref_set '0,10024' -path_in 'thaliana_genomes' -n_chr_query 5 -n_chr_ref 5

EOF
}

# ----------------------------------------------------------------------------
#             MAIN
# ----------------------------------------------------------------------------


# ref_set=""
additional_params=""

# Iterate through all arguments
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -h | --help )  print_usage
                     exit
                     ;;
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
    -cores)
        # Get the value of the pref_global parameter
        cores="$2"
        additional_params+=" $1 $2"
        shift # Skip the value of the parameter
        shift # Skip the parameter itself
        ;;
    -path_consensus)
        # Get the value of the pref_global parameter
        path_consensus="$2"
        additional_params+=" $1 $2"
        shift # Skip the value of the parameter
        shift # Skip the parameter itself
        ;;
    -path_chr_acc)
        path_chr_acc="$2"
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


if [ -z "$cores" ]; then
    cores=1
fi
pokaz_message "Number of cores ${cores}"


# Split the value of ref_set into separate words
IFS=',' read -ra refs_all <<< "$ref_set"

# Iterate over each word in ref_set
for ref0 in "${refs_all[@]}"; do
    
    command="./pangen_ref.sh -ref_pref ${ref0} ${additional_params}"
    # echo "Executing command: ${command}"
    eval "${command}"

done

# ----------------------------------------------------------------------------
#             COMMON CONSENSUS
# ----------------------------------------------------------------------------


check_missing_variable "pref_global"
pref_global=$(add_symbol_if_missing "$pref_global" "/")


path_consensus="${path_consensus:-${pref_global}consensus/}"
path_consensus=$(add_symbol_if_missing "$path_consensus" "/")

pokaz_message "Path with consensus ${path_consensus}"



ref0=${refs_all[0]}
ref0=${ref0//_/$'-'}
for ((i = 1; i < ${#refs_all[@]}; i++)); do
    ref1=${refs_all[i]}

    ref1=${ref1//_/$'-'}
    
    Rscript pangen/comb_02_two_refs.R --path.cons ${path_consensus} --ref0 ${ref0} --ref1 ${ref1} --cores ${cores}

done


Rscript pangen/comb_03_find_gaps.R --path.cons ${path_consensus} --ref.pref ${ref0} --cores ${cores}



# exit 1
pref_mafftin="${pref_global}mafft_in/"
if [ ! -d "$pref_mafftin" ]; then
    mkdir -p "$pref_mafftin"
fi


path_chr_acc="${path_chr_acc:-${pref_global}chromosomes/}"
path_chr_acc=$(add_symbol_if_missing "$path_chr_acc" "/")


Rscript pangen/comb_04_prepare_aln.R --path.cons ${path_consensus} --ref.pref ${ref0} --cores ${cores} \
                  --path.chromosomes ${path_chr_acc} --path.mafft.in ${pref_mafftin}

pref_mafft_out="${pref_global}mafft_out/"
if [ ! -d "$pref_mafft_out" ]; then
    mkdir -p "$pref_mafft_out"
fi


./pangen/comb_05_run_mafft.sh  --cores ${cores} \
                  --path.mafft.in ${pref_mafftin} \
                  --path.mafft.out ${pref_mafft_out} \


Rscript pangen/comb_06_final_aln.R  --cores ${cores}  --ref.pref ${ref0} \
                  --path.mafft.in ${pref_mafftin} \
                  --path.mafft.out ${pref_mafft_out} \
                  --path.cons ${path_consensus} 


