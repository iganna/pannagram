# This escript performs all stages of the alignment of genomes, 
# when the reference genome is already identified


# ./pipeline_consensus.sh -path_out '../pan_test/anna_ly' -ref_set 'NT1_220222,TE11_final' -n_chr_ref 8 -path_in '../lyrata_updated/' -n_chr_query 8 -cores 20  -one2one


# ./pipeline_consensus.sh -path_out '../pan_test/anna_aln' -ref_set 'NT1_220222,TE11_final' -n_chr_ref 8 -path_in '../lyrata_updated/' -n_chr_query 8 -cores 20  -one2one


# ./pipeline_consensus.sh -path_out '../pan_test/tom' -ref_set '0,6046-v1.1,6191-v1.1' -n_chr_ref 5 -path_in '../pb_updated/' -n_chr_query 5 -cores 30 -acc_anal acc_tom.txt -one2one


# ./pangen.sh  -path_in /Volumes/Samsung_T5/vienn/pannagram_test/symA  -path_out /Volumes/Samsung_T5/vienn/pannagram_test/symA_test -nchr_ref 1 -nchr_query 1 -nrefs 3

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


# Function to display help message
print_usage() {
    pokaz_help

    cat << EOF
Usage: ${0##*/}  -path_in INPUT_FOLDER  -path_out OUTPUT_FOLDER
                [-refs REF_NAME] [-nref NUM_OF_REFS] [-nchr NUM_OF_CHRS]
                [-path_ref PATH_CHR_REF] [-path_chrom PATH_CHROM] 
                [-path_parts PATH_PARTS] [-path_cons PATH_CONSENSUS] 
                [-sort_len] [-all2all] [-accessions ACC_FILE] 
                [-part_len PART_LEN] [-p_ident P_IDENT] [-purge_repeats]
                [-h] [-s STAGE] [-cores CORES] [-echo]
                
This script performs alignment of query genomes to the reference genome.

Options:
    -h, -help                  Display this help message and exit.
    -s, -stage STAGE            Specify the stage from which to run: a number from 1 to 12.
                                If not provided, the last interrupted stage will be re-run.
    -cores CORES                Number of cores for parallel processing. Default is 1.
    -echo                       Set up this flag to see the progress when number of cores is 1.

        
    * Required parameters:
        -path_in INPUT_FOLDER       Folder to the directory with query genomes. 
                                    Possible file types: fasta, fna, fa, fas.
        -path_out OUTPUT_FOLDER     Folder where all results and intermediate files will appear.


    * Optional parameters:
        -refs REF_NAMES             Names of reference genomes for randomisation.
                                    Names should NOT contain "_chr" as a substring.
        -nref NUM_OF_REFS           Number of reference genomes for randomising,
                                    when -refs is not specified.
        -nchr NUM_OF_CHRS           Number of chromosomes in both reference and query genomes.

    * Optional paths:
        -path_ref PATH_CHR_REF      Path where the reference genome is stored. 
                                    Do not provide if it's the same folder as the path with query genomes.
        -path_chrom PATH_CHROM      Path to the folder with individual chromosomes in separate files. 
        -path_parts PATH_PARTS      Path to the folder with files of chromosomal parts.
        -path_cons PATH_CONSENSUS   Path to the consensus folder.

    * Input Design Handling:
        -sort_len                   Flag to sort chromosomes by length.
        -one2one                    Flag for straightforward pairwise alignment of chromosomes,
                                    matching each chromosome sequentially with its corresponding one 
                                    (the default).
        -all2all                    Flag for aligning all chromosomes to all.
        -accessions ACC_FILE        File with accessions to analyze. Accessions should be in rows.
        -combinations COMB_FILE     File with combinations to analyze.

    * Tuning parameters: 
        -p_ident P_IDENT            Percentage identity threshold (default: 85).
        -part_len PART_LEN          Fragments to which each chromosome should be cut (default value: 5000).
        -purge_repeats              Enable filtration of repeats (default: disabled).

    
Examples:
    ${0##*/} -path_in 'input_genomes' -path_out 'output_folder'

EOF
}



# ----------------------------------------------------------------------------
#             PARAMETERS
# ----------------------------------------------------------------------------

# ref_set=""
ref_num=2
additional_params=""
# all2all=""
# one2one=""

# Iterate through all arguments
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -h | -help  )  print_usage
                     exit
                     ;;
    -s | -stage )
        # Step, from which to start calculation
        start_step="$2"
        additional_params+=" $1 $2"
        shift 2
        ;;
    -cores)
        cores="$2"
        additional_params+=" $1 $2"
        shift 2
        ;;  
    -log) 
        log_level=$2
        additional_params+=" $1 $2"
        shift 2   
        ;;      
    -refs)
        ref_set="$2"
        shift 2
        ;;
    -nref)
        ref_num="$2"
        shift 2
        ;;
    -nchr)
        chr_num="$2"
        additional_params+=" -nchr_ref ${chr_num} -nchr_query ${chr_num}"
        shift 2
        ;;
    -nchr_ref) # ignore this parameter
        pokaz_attention "Parameter -nchr_ref is ignored."
        pokaz_attention "Please use -nchr or omit the parameter for automatic chromosome count detection."
        shift 2
        ;;
    -nchr_query)
        pokaz_attention "Parameter -nchr_query is ignored."
        pokaz_attention "Please use -nchr or omit the parameter for automatic chromosome count detection."
        shift 2
        ;;
    -path_in)
        path_in="$2"
        additional_params+=" $1 $2"
        shift 2
        ;;
    -path_out)
        path_out="$2"
        path_out=$(add_symbol_if_missing "$path_out" "/")
        additional_params+=" $1 $2"
        shift 2
        ;;
    -path_cons)
        path_consensus="$2"
        additional_params+=" $1 $2"
        shift 2
        ;;
    -path_chrom)
        path_chr_acc="$2"
        additional_params+=" $1 $2"
        shift 2
        ;;
    -all2all)
        all2all="T"  # if this exists then pangen_ref should be run with the default mode not -one2one
        shift 1
        ;;
    -one2one)
        one2one="T"  # if this exists then pangen_ref should be run with -one2one mode
        shift 1
        ;;
    *)
        # Add the remaining arguments as additional parameters
        additional_params+=" $1"
        shift
        ;;
  esac
done


# Decide one2one or all2all mode

if [[ -n "$all2all" && -n "$one2one" ]]; then
    pokaz_error "Error: both parameters -all2all and -one2one cannot be set simultaneously."
    exit 1
elif [[ -z "$all2all" ]]; then
    additional_params+=" -one2one"  # Default, when no parameters are provided
fi


# ----------------------------------------------------------------------------
#           LOGS
# ----------------------------------------------------------------------------

log_level=${log_level:-1}  # Set the default value to 'pipeline'

if ! [[ "$log_level" =~ ^[0-3]$ ]]; then
    pokaz_error "Error: log_level must be a number between 0 and 3."
    exit 1
fi

# Hidden path with logs
path_logs="${path_out}logs/"
make_dir ${path_logs}

# Path for logs of this script
path_log="${path_logs}pangen/"
make_dir ${path_log}

# File with pipeline logs
file_log="${path_log}pipeline.log"
> "${file_log}"

# ----------------------------------------------------------------------------
#             CHECK PARAMETERS
# ----------------------------------------------------------------------------

# ----------------------------------------------
# Number of cores
if [ -z "$cores" ]; then
    cores=1
fi
log_message 1 "$log_level" "$file_log" \
    pokaz_message "Number of cores: ${cores}"

# ----------------------------------------------
# Check the number of reference genomes for randomisation
if ! [[ "$ref_num" =~ ^[0-9]+$ ]]; then
    log_message 0 "$log_level" "$file_log" \
        pokaz_error "Error: Number of reference genomes is not a number."
    exit 1
fi
log_message 1 "$log_level" "$file_log" \
    pokaz_message "Number of genomes for randomisation: ${ref_num}"


# ----------------------------------------------
# Get names of all genomes form the input folder
acc_set=($(find "$path_in" -type f -exec bash -c 'basename "$1" .${1##*.}' _ {} \;))

# # Print the genomes for verification
# for ref in "${acc_set[@]}"; do
#     echo "$ref"
# done

# Check if ref_num is greater than the number of genomes in acc_set
if (( ref_num > ${#acc_set[@]} )); then
    log_message 0 "$log_level" "$file_log" \
        pokaz_error "Error: ref_num ($ref_num) is greater than the number of available genomes (${#acc_set[@]}) in acc_set."
    exit 1
fi


# ----------------------------------------------
# Check if reference genomes are set up
if [[ -z "$ref_set" ]]; then

    # Take the first ref_num genomes
    refs_all=("${acc_set[@]:0:$ref_num}")

else 
    # Split the value of ref_set into separate words
    IFS=',' read -ra refs_all <<< "$ref_set"

    # Check if the number of reference genomes is sufficient
    if (( ${#refs_all[@]} < ref_num )); then
        genomes_needed=$((ref_num - ${#refs_all[@]}))
        log_message 0 "$log_level" "$file_log" \
            pokaz_attention "Not enough reference genomes. Adding $genomes_needed genome(s)."

        # Add genomes from acc_set to refs_all to satisfy the number of ref_num, ensuring no repeats
        for genome in "${acc_set[@]}"; do
            if (( ${#refs_all[@]} >= ref_num )); then
                break
            fi
            if [[ ! " ${refs_all[@]} " =~ " ${genome} " ]]; then
                refs_all+=("$genome")
            fi
        done
    fi

    if (( ${#refs_all[@]} > ref_num )); then
        log_message 0 "$log_level" "$file_log" \
            pokaz_attention "Parameter -nref is ignored; it doesn't match number of genome in -refs"
    fi
fi



# Print the selected reference genomes for verification
log_message 1 "$log_level" "$file_log" \
    pokaz_message "Names of genomes for randomisation: $(IFS=,; echo "${refs_all[*]}")"

# Check that all of the genomes are in the folder
for ref in "${refs_all[@]}"; do
    if [[ ! " ${acc_set[@]} " =~ " ${ref} " ]]; then
        log_message 0 "$log_level" "$file_log" \
            pokaz_attention "Error: Genome ${ref} is not in ${path_in}"
        exit 1
    fi
done

# ----------------------------------------------
# Check number of chromosomes
if [ -z "$chr_num" ]; then
    # Automatic chromosomal number detection
    file_ref=$(find "$path_in" -type f -name "${refs_all[0]}.*")
    chr_num=$(grep -c '^>' "$file_ref")

    log_message 0 "$log_level" "$file_log" \
        pokaz_attention "Number of chromosomes identified: ${chr_num}"
    additional_params+=" -nchr_ref ${chr_num} -nchr_query ${chr_num}"
else
    log_message 0 "$log_level" "$file_log" \
        pokaz_message "Number of chromosomes to analyse: ${chr_num}"    
fi


# ----------------------------------------------
# Consensus paths

check_missing_variable "path_out"
path_out=$(add_symbol_if_missing "$path_out" "/")

path_consensus="${path_consensus:-${path_out}consensus/}"
path_consensus=$(add_symbol_if_missing "$path_consensus" "/")

log_message 0 "$log_level" "$file_log" \
    pokaz_message "Path with consensus MSA: ${path_consensus}"

# ----------------------------------------------
# Path with chromosomes
path_chr_acc="${path_chr_acc:-${path_out}chromosomes/}"
path_chr_acc=$(add_symbol_if_missing "$path_chr_acc" "/")

# ----------------------------------------------------------------------------
#             CHECK STEPS
# ----------------------------------------------------------------------------

# Path with steps
path_flags="${path_out}.flags/"
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


# if [ "$start_step" -eq 3 ]; then
#     log_message 1 "$log_level" "$file_log" \
#         pokaz_message "Starting step: reference cycles 3-7"
# else
#     log_message 1 "$log_level" "$file_log" \
#         pokaz_message "Starting step: ${start_step}"
# fi


# ----------------------------------------------------------------------------
#             ONE REFERENCE
# ----------------------------------------------------------------------------


# Iterate over each word in ref_set
for ref0 in "${refs_all[@]}"; do
    
    command="./pangen_ref.sh -ref ${ref0} ${additional_params}"
    # echo "Executing command: ${command}"
    eval "${command}"


    # File with pipeline logs
    file_log_ref="${path_logs}pangen_ref_${ref0}/pipeline.log"

    cat ${file_log_ref} >> ${file_log}

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

        log_message 1 "$log_level" "$file_log" \
            pokaz_stage "Step ${step_num}. Randomisation of alignments. Combine two references: ${ref0} and ${ref1}."
        # ref1=${ref1//_/$'-'}
        
        Rscript pangen/comb_02_two_refs.R --path.cons ${path_consensus} --ref0 ${ref0} --ref1 ${ref1} --cores ${cores}

    done

    touch "$path_flags/step${step_num}_done"
    log_message 1 "$log_level" "$file_log" \
        pokaz_message "Done!"
fi

((step_num = step_num + 1))

# ----------------------------------------------
if [ $start_step -le ${step_num} ] && [ ! -f "$path_flags/step${step_num}_done" ]; then

    log_message 1 "$log_level" "$file_log" \
        pokaz_stage "Step ${step_num}. Find Positions of Common Gaps in the Reference-Free Multiple Genome Alignment."

    Rscript pangen/comb_03_find_gaps.R --path.cons ${path_consensus} --ref.pref ${ref0} --cores ${cores}

    touch "$path_flags/step${step_num}_done"
    log_message 1 "$log_level" "$file_log" \
        pokaz_message "Done!"
fi

((step_num = step_num + 1))

# ----------------------------------------------
# Create sequences to run MAFFT and perform some small alignments
pref_mafftin="${path_out}mafft_in/"
if [ ! -d "$pref_mafftin" ]; then
    mkdir -p "$pref_mafftin"
fi

if [ $start_step -le ${step_num} ] && [ ! -f "$path_flags/step${step_num}_done" ]; then

    log_message 1 "$log_level" "$file_log" \
        pokaz_stage "Step ${step_num}. Prepare sequences for MAFFT."

    Rscript pangen/comb_04_prepare_aln.R --path.cons ${path_consensus} --ref.pref ${ref0} --cores ${cores} \
                      --path.chromosomes ${path_chr_acc} --path.mafft.in ${pref_mafftin}

    touch "$path_flags/step${step_num}_done"
    log_message 1 "$log_level" "$file_log" \
        pokaz_message "Done!"
fi

((step_num = step_num + 1))

# ----------------------------------------------
# Run MAFFT
pref_mafft_out="${path_out}mafft_out/"
if [ ! -d "$pref_mafft_out" ]; then
    mkdir -p "$pref_mafft_out"
fi

if [ $start_step -le ${step_num} ] && [ ! -f "$path_flags/step${step_num}_done" ]; then

    log_message 1 "$log_level" "$file_log" \
        pokaz_stage "Step ${step_num}. Run MAFFT."

    ./pangen/comb_05_run_mafft.sh  --cores ${cores} \
                      --path.mafft.in ${pref_mafftin} \
                      --path.mafft.out ${pref_mafft_out}

    touch "$path_flags/step${step_num}_done"
    log_message 1 "$log_level" "$file_log" \
        pokaz_message "Done!"
fi

((step_num = step_num + 1))

# ----------------------------------------------
# Combine all together
if [ $start_step -le ${step_num} ] && [ ! -f "$path_flags/step${step_num}_done" ]; then

    log_message 1 "$log_level" "$file_log" \
        pokaz_message "Step ${step_num}. Combine all alignments together into the final one."

    Rscript pangen/comb_06_final_aln.R  --cores ${cores}  --ref.pref ${ref0} \
                      --path.mafft.in ${pref_mafftin} \
                      --path.mafft.out ${pref_mafft_out} \
                      --path.cons ${path_consensus} 

    touch "$path_flags/step${step_num}_done"
    log_message 1 "$log_level" "$file_log" \
        pokaz_message "Done!"
fi

((step_num = step_num + 1))


# ----------------------------------------------
# Get synteny blocks

aln_type='msa_'
if [ $start_step -le ${step_num} ] && [ ! -f "$path_flags/step${step_num}_done" ]; then

    log_message 1 "$log_level" "$file_log" \
        pokaz_message "Step ${step_num}. Get synteny blocks."

    Rscript analys/analys_01_blocks.R --path.cons ${path_consensus} --ref.pref  ${ref0} --cores ${cores} --aln.type ${aln_type}

    touch "$path_flags/step${step_num}_done"
    log_message 1 "$log_level" "$file_log" \
        pokaz_message "Done!"
fi

((step_num = step_num + 1))

# if [ $start_step -eq 0 ]; then
#     rm -f "$FLAG_DIR"/.*
#     echo "Script completed successfully"
# fi