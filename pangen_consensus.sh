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
#             PARAMETERS
# ----------------------------------------------------------------------------

# ref_set=""
additional_params=""
start_step=1

# Iterate through all arguments
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -h | -help | --help )  print_usage
                     exit
                     ;;
    -s | -stage )
        # Step, from which to start calculation
        start_step="$2"
        additional_params+=" $1 $2"
        shift 2
        ;;                     
    -ref_set)
        # Get the value of the ref_set parameter
        ref_set="$2"
        shift 2
        ;;
    -pref_global)
        # Get the value of the pref_global parameter
        pref_global="$2"
        additional_params+=" $1 $2"
        shift 2
        ;;
    -cores)
        # Get the value of the pref_global parameter
        cores="$2"
        additional_params+=" $1 $2"
        shift 2
        ;;
    -path_consensus)
        # Get the value of the pref_global parameter
        path_consensus="$2"
        additional_params+=" $1 $2"
        shift 2
        ;;
    -path_chr_acc)
        path_chr_acc="$2"
        additional_params+=" $1 $2"
        shift 2
        ;;
    *)
        # Add the remaining arguments as additional parameters
        additional_params+=" $1"
        shift
        ;;
  esac
done

pokaz_message "Sterting ste: ${start_step}"

# ----------------------------------------------------------------------------
#             CHECK PARAMETERS
# ----------------------------------------------------------------------------


# Check if ref_set is set
if [[ -z "$ref_set" ]]; then
  echo "Error: -ref_set parameter is not set" >&2
  exit 1
fi


if [ -z "$cores" ]; then
    cores=1
fi
pokaz_message "Number of cores ${cores}"


# Consensus paths

check_missing_variable "pref_global"
pref_global=$(add_symbol_if_missing "$pref_global" "/")


path_consensus="${path_consensus:-${pref_global}consensus/}"
path_consensus=$(add_symbol_if_missing "$path_consensus" "/")

pokaz_message "Path with consensus ${path_consensus}"

# Path with chromosomes
path_chr_acc="${path_chr_acc:-${pref_global}chromosomes/}"
path_chr_acc=$(add_symbol_if_missing "$path_chr_acc" "/")

# ----------------------------------------------------------------------------
#             CHECK STEPS
# ----------------------------------------------------------------------------

# Path with steps
path_flags="${pref_global}flags/"
path_flags=$(add_symbol_if_missing "$path_flags" "/")
if [ ! -d "$path_flags" ]; then
    mkdir -p "$path_flags"
fi

# Looping through and deleting files
for file_step in "$path_flags"step*_done; do
    if [ -f "$file_step" ]; then
        # Extracting step number from the file name
        step_tmp=$(echo "$file_step" | sed -e 's/.*step\([0-9]*\)_done/\1/')

        # Check if step number is greater or equal to start_step
        if [ "$step_tmp" -ge "$start_step" ]; then
            echo "Deleting file: $file_step"
            rm -f "$file_step"
        fi
    fi
done

# ----------------------------------------------------------------------------
#             ONE REFERENCE
# ----------------------------------------------------------------------------

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


# Run consensus for a pair of files
if [ $start_step -le 8 ] && [ ! -f "$path_flags/step8_done" ]; then

    ref0=${refs_all[0]}
    ref0=${ref0//_/$'-'}
    for ((i = 1; i < ${#refs_all[@]}; i++)); do
        ref1=${refs_all[i]}

        ref1=${ref1//_/$'-'}
        
        Rscript pangen/comb_02_two_refs.R --path.cons ${path_consensus} --ref0 ${ref0} --ref1 ${ref1} --cores ${cores}

    done

    touch "$path_flags/step8_done"
fi



if [ $start_step -le 9 ] && [ ! -f "$path_flags/step9_done" ]; then

    Rscript pangen/comb_03_find_gaps.R --path.cons ${path_consensus} --ref.pref ${ref0} --cores ${cores}

    touch "$path_flags/step9_done"
fi


# Create sequences to run MAFFT and perform some small alignments
if [ $start_step -le 10 ] && [ ! -f "$path_flags/step10_done" ]; then

    pref_mafftin="${pref_global}mafft_in/"
    if [ ! -d "$pref_mafftin" ]; then
        mkdir -p "$pref_mafftin"
    fi

    Rscript pangen/comb_04_prepare_aln.R --path.cons ${path_consensus} --ref.pref ${ref0} --cores ${cores} \
                      --path.chromosomes ${path_chr_acc} --path.mafft.in ${pref_mafftin}

    touch "$path_flags/step10_done"
fi


# Run MAFFT
if [ $start_step -le 11 ] && [ ! -f "$path_flags/step11_done" ]; then

    pref_mafft_out="${pref_global}mafft_out/"
    if [ ! -d "$pref_mafft_out" ]; then
        mkdir -p "$pref_mafft_out"
    fi


    ./pangen/comb_05_run_mafft.sh  --cores ${cores} \
                      --path.mafft.in ${pref_mafftin} \
                      --path.mafft.out ${pref_mafft_out}

    touch "$path_flags/step11_done"
fi


# Combine all together
if [ $start_step -le 12 ] && [ ! -f "$path_flags/step12_done" ]; then

    Rscript pangen/comb_06_final_aln.R  --cores ${cores}  --ref.pref ${ref0} \
                      --path.mafft.in ${pref_mafftin} \
                      --path.mafft.out ${pref_mafft_out} \
                      --path.cons ${path_consensus} 

    touch "$path_flags/step12_done"
fi


if [ $start_step -eq 0 ]; then
    rm -f "$FLAG_DIR"/.*
    echo "Скрипт успешно завершен"
fi