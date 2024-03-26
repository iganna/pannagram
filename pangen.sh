# This escript performs all stages of the alignment of genomes, 
# when the reference genome is already identified


# ./pipeline_consensus.sh -path_out '../pan_test/anna_ly' -ref_set 'NT1_220222,TE11_final' -n_chr_ref 8 -path_in '../lyrata_updated/' -n_chr_query 8 -cores 20  -one2one


# ./pipeline_consensus.sh -path_out '../pan_test/anna_aln' -ref_set 'NT1_220222,TE11_final' -n_chr_ref 8 -path_in '../lyrata_updated/' -n_chr_query 8 -cores 20  -one2one


# ./pipeline_consensus.sh -path_out '../pan_test/tom' -ref_set '0,6046-v1.1,6191-v1.1' -n_chr_ref 5 -path_in '../pb_updated/' -n_chr_query 5 -cores 30 -acc_anal acc_tom.txt -one2one

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
Usage: ${0##*/} [-h] [-s STAGE] [-cores CORES] 
                [-path_out OUTPUT_FOLDER] [-path_in INPUT_FOLDER] [-ref_set REF_NAMES] 
                [-nchr_ref N_CHR_REF] [-nchr_query N_CHR_QUERY] 
                [-path_ref PATH_CHR_REF] [-path_chrom PATH_CHROM] [-path_parts PATH_PARTS] [-path_cons PATH_CONSENSUS] 
                [-sort_len] [-one2one] [-accessions ACC_ANAL] 
                [-part_len PART_LEN] [-p_ident P_IDENT] [-purge_repeats]
                
This script performs alignment of query genomes to the reference genome.

Options:
    -h, --help                  Display this help message and exit.
    -s, -stage STAGE            Specify the stage from which to run. If not provided, the last interrupted stage will be re-run.
    -cores CORES                Number of cores for parallel processing. Default is 1.

    # Required parameters:
    -path_out OUTPUT_FOLDER     Folder where all results and intermediate files will appear.
    -path_in INPUT_FOLDER       Folder to the directory with query genomes. Possible file types: fasta, fna, fa, fas.
    -ref_set REF_NAMES          Names of the reference genomes.

    # Numbers of chromosomes:
    -nchr_ref N_CHR_REF         Number of chromosomes in the reference genome.
    -nchr_query N_CHR_QUERY     Number of chromosomes in the query genome.

    # Optional paths:
    -path_ref PATH_CHR_REF      Path where the reference genome is stored. Do not provide if it's the same folder as the path with query genomes.
    -path_chrom PATH_CHROM      Path to the folder with individual chromosomes in separate files. 
    -path_parts PATH_PARTS      Path to the folder with files of chromosomal parts.
    -path_cons PATH_CONSENSUS   Path to the consensus folder.

    # Optional Input Design Handling:
    -sort_len                   Flag to sort chromosomes by length.
    -one2one                    Flag for straightforward pairwise alignment of chromosomes, matching each chromosome sequentially with its corresponding one (NOT all vs all, as by default).
    -accessions ACC_ANAL        File with accessions to analyze. Accessions should be in rows.

    # Optional Tuning parameters:
    -p_ident P_IDENT            Percentage identity threshold (default: 85).
    -part_len PART_LEN          Fragments to which each chromosome should be cut (default value: 5000).
    -purge_repeats              Enable filtration of repeats (default: disabled).

    
Examples:
    ${0##*/} -path_out 'output_folder' -path_in 'input_genomes' -ref_set '0,10024' -nchr_query 5 -nchr_ref 5

EOF
}



# ----------------------------------------------------------------------------
#             PARAMETERS
# ----------------------------------------------------------------------------

# ref_set=""
additional_params=""

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
    -refs)
        # Get the value of the ref_set parameter
        ref_set="$2"
        shift 2
        ;;
    -path_out)
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
    -path_cons)
        # Get the value of the pref_global parameter
        path_consensus="$2"
        additional_params+=" $1 $2"
        shift 2
        ;;
    -path_chrom)
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
pokaz_message "Number of cores: ${cores}"


# Consensus paths

check_missing_variable "pref_global"
pref_global=$(add_symbol_if_missing "$pref_global" "/")


path_consensus="${path_consensus:-${pref_global}consensus/}"
path_consensus=$(add_symbol_if_missing "$path_consensus" "/")

pokaz_message "Path with consensus: ${path_consensus}"

# Path with chromosomes
path_chr_acc="${path_chr_acc:-${pref_global}chromosomes/}"
path_chr_acc=$(add_symbol_if_missing "$path_chr_acc" "/")

# ----------------------------------------------------------------------------
#             CHECK STEPS
# ----------------------------------------------------------------------------

# Path with steps
path_flags="${pref_global}.flags/"
path_flags=$(add_symbol_if_missing "$path_flags" "/")
if [ ! -d "$path_flags" ]; then
    mkdir -p "$path_flags"
fi


if [ -z "${start_step}" ]; then
    max_step_file=$(ls ${path_flags}step*_done 2> /dev/null | sort -V | tail -n 1)

    if [ -z "$max_step_file" ]; then
        start_step=1
    else
        start_step=${max_step_file##*step}
        start_step=${start_step%_done}
        start_step=$((start_step + 1))
    fi

fi

# Looping through and deleting files of stages, which are less then the current one
for file_step in "${path_flags}"step*_done*; do
    # echo "Processing file: $file_step"
    if [ -f "$file_step" ]; then
        # Extracting step number from the file name
        step_tmp=$(echo "$file_step" | sed -e 's/.*step\([0-9]*\)_done.*/\1/')

        # Check if step number is greater or equal to start_step
        if [ "$step_tmp" -ge "$start_step" ]; then
            # echo "remove ${file_step}"
            rm -f "$file_step"
        fi
    fi
done


if [ "$start_step" -eq 3 ]; then
    pokaz_message "Starting step: reference cycles 3-7"
else
    pokaz_message "Starting step: ${start_step}"
fi


# ----------------------------------------------------------------------------
#             ONE REFERENCE
# ----------------------------------------------------------------------------

# Split the value of ref_set into separate words
IFS=',' read -ra refs_all <<< "$ref_set"


# Iterate over each word in ref_set
for ref0 in "${refs_all[@]}"; do
    
    command="./pangen_ref.sh -ref ${ref0} ${additional_params}"
    echo "Executing command: ${command}"
    eval "${command}"

done

for tmp in {1..7}; do
    touch "${path_flags}/step${tmp}_done"
done

# ----------------------------------------------------------------------------
#             COMMON CONSENSUS
# ----------------------------------------------------------------------------

step_num=9

# ----------------------------------------------
# Run consensus for a pair of files
ref0=${refs_all[0]}

if [ $start_step -le ${step_num} ] && [ ! -f "$path_flags/step${step_num}_done" ]; then

    for ((i = 1; i < ${#refs_all[@]}; i++)); do
        ref1=${refs_all[i]}

        pokaz_stage "Step ${step_num}. Randomisation of alignments. Combine two references: ${ref0} and ${ref1}."
        # ref1=${ref1//_/$'-'}
        
        Rscript pangen/comb_02_two_refs.R --path.cons ${path_consensus} --ref0 ${ref0} --ref1 ${ref1} --cores ${cores}

    done

    touch "$path_flags/step${step_num}_done"
    pokaz_message "Done!"
fi

((step_num = step_num + 1))

# ----------------------------------------------
if [ $start_step -le ${step_num} ] && [ ! -f "$path_flags/step${step_num}_done" ]; then

    pokaz_stage "Step ${step_num}. Find Positions of Common Gaps in the Reference-Free Multiple Genome Alignment."

    Rscript pangen/comb_03_find_gaps.R --path.cons ${path_consensus} --ref.pref ${ref0} --cores ${cores}

    touch "$path_flags/step${step_num}_done"
    pokaz_message "Done!"
fi

((step_num = step_num + 1))

# ----------------------------------------------
# Create sequences to run MAFFT and perform some small alignments
pref_mafftin="${pref_global}mafft_in/"
if [ ! -d "$pref_mafftin" ]; then
    mkdir -p "$pref_mafftin"
fi

if [ $start_step -le ${step_num} ] && [ ! -f "$path_flags/step${step_num}_done" ]; then

    pokaz_stage "Step ${step_num}. Prepare sequences for MAFFT."

    Rscript pangen/comb_04_prepare_aln.R --path.cons ${path_consensus} --ref.pref ${ref0} --cores ${cores} \
                      --path.chromosomes ${path_chr_acc} --path.mafft.in ${pref_mafftin}

    touch "$path_flags/step${step_num}_done"
    pokaz_message "Done!"
fi

((step_num = step_num + 1))

# ----------------------------------------------
# Run MAFFT
pref_mafft_out="${pref_global}mafft_out/"
if [ ! -d "$pref_mafft_out" ]; then
    mkdir -p "$pref_mafft_out"
fi

if [ $start_step -le ${step_num} ] && [ ! -f "$path_flags/step${step_num}_done" ]; then

    pokaz_stage "Step ${step_num}. Run MAFFT."

    ./pangen/comb_05_run_mafft.sh  --cores ${cores} \
                      --path.mafft.in ${pref_mafftin} \
                      --path.mafft.out ${pref_mafft_out}

    touch "$path_flags/step${step_num}_done"
    pokaz_message "Done!"
fi

((step_num = step_num + 1))

# ----------------------------------------------
# Combine all together
if [ $start_step -le ${step_num} ] && [ ! -f "$path_flags/step${step_num}_done" ]; then

    pokaz_message "Step ${step_num}. Combine all alignments together into the final one."
    
    Rscript pangen/comb_06_final_aln.R  --cores ${cores}  --ref.pref ${ref0} \
                      --path.mafft.in ${pref_mafftin} \
                      --path.mafft.out ${pref_mafft_out} \
                      --path.cons ${path_consensus} 

    touch "$path_flags/step${step_num}_done"
    pokaz_message "Done!"
fi

((step_num = step_num + 1))

# if [ $start_step -eq 0 ]; then
#     rm -f "$FLAG_DIR"/.*
#     echo "Script completed successfully"
# fi