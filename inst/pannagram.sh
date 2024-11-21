#!/bin/bash
# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------
INSTALLED_PATH=$(Rscript -e "cat(system.file(package = 'pannagram'))")
source $INSTALLED_PATH/utils/chunk_error_control.sh

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source $INSTALLED_PATH/utils/utils_bash.sh
source $INSTALLED_PATH/utils/utils_help.sh

# ----------------------------------------------------------------------------
#            PARAMETERS: parsing
# ----------------------------------------------------------------------------

aln_type_msa='msa_'
aln_type_ref='ref_'
aln_type_comb='comb_'

mode_pre="F"
mode_ref="F"
mode_msa="F"
clean="F"
one_step="F"
purge_reps="F"
rm_inter="F"  # remove intermediate files

unrecognized_options=()

step_start=0
step_end=100
while [ $# -gt 0 ]
do
    # echo $1
    case $1 in
        -h) print_usage_short; print_examples; exit ;;
        -help ) print_usage_detailed; print_examples; exit ;;
        -s | -stage ) step_start="$2"; shift 2 ;;  # stage from which to run, when the stage is not provided - the last interrupted stage withh be re-run
        -e | -end )   step_end="$2"; shift 2 ;;  # stage from which to run, when the stage is not provided - the last interrupted stage withh be re-run
        -log)         log_level=$2;    shift 2 ;;  # path to the output
        -cores)       cores=$2;        shift 2 ;;
        -clean    | -cleanup) clean="T";       shift 1 ;;
        -one_step | -1s)      one_step="T";    shift 1 ;;
        -rm_inter)            rm_inter="T";    shift 1 ;;
        
        # Required
        -path_out) path_out=$2; shift 2 ;;  # path to the output
        -path_in)  path_in=$2;  shift 2 ;;  # path with all genomes in fasta format

        # REF-based
        -pre)                   mode_pre="T"; shift 1 ;;  # shitch to the preliminary mode
        -ref)      ref_name=$2; mode_ref="T"; shift 2 ;;  # name of the reference genome
        -path_ref) path_ref=$2;               shift 2 ;;  # dont provide if it's the same folder as the path with query genomes

        # MSA
        -refs) ref_set="$2"; mode_msa="T"; shift 2 ;;
        -nref) ref_num="$2"; mode_msa="T"; shift 2 ;;

        # # Optional paths
        # -path_chrom) path_chrom=$2;  shift 2 ;;  # path to the folder with individual chromosomes in separate files
        # -path_parts) path_parts=$2;  shift 2 ;;  # path to the folder with chromosomal parts
        # -path_cons)  path_cons=$2;   shift 2 ;;  # path to the consensus folder

        # Number of chromosomes
        -nchr)      nchr=$2;     shift 2 ;;  # number of chromosomes
        -nchr_ref)  nchr_ref=$2; shift 2 ;;  # number of chromosome in the query genome

        -part_len)  part_len=$2; shift 2 ;;  # fragments to which each chromosome should be cut, has a default value 5000
        -p_ident)   p_ident=$2;  shift 2 ;;  # percent of identity
        -p_ident_gap)   p_ident=$2;  shift 2 ;;  # percent of identity
        -max_len_gap)  max_len_gap=$2; shift 2 ;;  # Max length that can be aligned with MAFFT

        -one2one)     one2one="T";      shift 1 ;;  # compare chromosomes one-to-one or not (default in REF and MSA modes)
        -all2all)     one2one="F";      shift 1 ;;  # compare chromosomes all-to-all or not (default in PRE mode)
        -purge_reps | -purge_repeats ) purge_reps="T";    shift 1 ;;  # filtration of repeats, default - not
        -rev )        flag_rev="T";      shift 1 ;;  # reverce parts
        -orf )        flag_orf="T";      shift 1 ;;  # reverce parts

        -accessions)   acc_file=$2;  shift 2 ;;  # file with accessions to analyse
        -combinations) comb_file=$2; shift 2 ;;  # file with chromosomal combinations to analyse: first column - query, second column - reference(base)
    
        *) unrecognized_options+=("$1"); shift 1 ;;
    esac
done

# Output of Unrecognized Parameters
if [[ ${#unrecognized_options[@]} -gt 0 ]]; then
    print_usage_short
    pokaz_message "Unrecognized options:"
    for option in "${unrecognized_options[@]}"; do
        echo "$option"
    done
    exit 1
fi


# ----------------------------------------------------------------------------
#            MODE
# ----------------------------------------------------------------------------

name_mode_pre='PRE'
name_mode_ref='REF'
name_mode_msa='MSA'

if   [ "$mode_pre" = "F" ] && [ "$mode_ref" = "F" ] && [ "$mode_msa" = "F" ] && [ -z "${path_ref}" ]; then  # Default MSA
    mode_pangen=${name_mode_msa}
elif [ "$mode_pre" = "F" ] && [ "$mode_ref" = "F" ] && [ "$mode_msa" = "T" ] && [ -z "${path_ref}" ]; then  # Only MSA parameters are setup
    mode_pangen=${name_mode_msa}
elif [ "$mode_pre" = "T" ] && [ "$mode_ref" = "T" ] && [ "$mode_msa" = "F" ]; then  # When -pre -ref are defined
    mode_pangen=${name_mode_pre}
elif [ "$mode_pre" = "F" ] && [ "$mode_ref" = "T" ] && [ "$mode_msa" = "F" ]; then  # When only -ref is defined
    mode_pangen=${name_mode_ref}
else
    pokaz_error "Error: Invalid combination of parameters to determine the pangen launch mode."
    print_usage_short
    exit 1
fi


if [[ ${#one2one} -gt 1 ]]; then
    pokaz_error "Error: -all2all and -one2one should not be set up together"
    exit 1
elif [[ ${#one2one} -eq 0 ]]; then
    # Parameters -all2all and -one2one should be taken with default values
    if [[ "$mode_pangen" == "${name_mode_pre}" ]]; then
        one2one="F"
    else
        one2one="T"
    fi
fi

# ----------------------------------------------------------------------------
#            PARAMETERS: checking
# ----------------------------------------------------------------------------

# ----------------------------------------------
# Required PATHS
if [ -z "${path_in}" ]; then
    pokaz_error "Error: The path to the genomes folder (-path_in) is not specified"
    exit 1
fi

if [ -z "${path_out}" ]; then
    pokaz_error "Error: The path to the output folder (-path_out) is not specified"
    exit 1
fi

path_in=$(add_symbol_if_missing "$path_in" "/")
path_out=$(add_symbol_if_missing "$path_out" "/")

# Optional paths paths
path_ref="${path_ref:-${path_in}}"
path_ref=$(add_symbol_if_missing "$path_ref" "/")

path_plots="${path_out}plots/"
path_inter="${path_out}intermediate/"
path_cons="${path_inter}consensus/"
path_chrom="${path_inter}chromosomes/"
path_parts="${path_inter}parts/"

# Make folders
mkdir -p "${path_out}"
mkdir -p "${path_plots}"
mkdir -p "${path_inter}"
mkdir -p "${path_cons}"
mkdir -p "${path_chrom}"
mkdir -p "${path_parts}"

# ----------------------------------------------
# Handling accessions
genome_extensions=('fasta' 'fna' 'fa' 'fas')

if [ ! -d "${path_in}" ]; then
    pokaz_error "Error: Directory ${path_in} does not exist."
    exit 1
fi

# All accessions from the folder
acc_set=()
for ext in "${genome_extensions[@]}"; do
    for file in "$path_in"/*."$ext"; do
        if [[ -f "$file" ]]; then
            basename=$(basename "$file" ".$ext")
            acc_set+=("$basename")
        fi
    done
done

# Check the input folder
if [ ${#acc_set[@]} -eq 0 ]; then
    pokaz_error "Error: No files with the specified extensions found in ${path_in}. Checked extensions: ${genome_extensions[*]}"
    exit 1
fi

# Target accessions from the file
acc_target=()
if [ -n "${acc_file}" ] && [ -f "${acc_file}" ]; then
    # Read names from the file into an array
    while IFS= read -r name; do
        acc_target+=("$name")
    done < "$acc_file"
fi

# Check if acc_target is not empty
if [ ${#acc_target[@]} -ne 0 ]; then
    # Find common elements between acc_set and acc_target
    common_acc=()
    for acc in "${acc_target[@]}"; do
        if [[ " ${acc_set[@]} " =~ " ${acc} " ]]; then
            common_acc+=("$acc")
        else
            pokaz_attention "Warning: ${acc} from acc_target is not found in acc_set"
        fi
    done

    # Update acc_set with common elements
    acc_set=("${common_acc[@]}")
fi

# Check Intersecton
if [ ${#acc_set[@]} -eq 0 ]; then
    pokaz_error "Error: No fcommon accessions in ${acc_file} file and ${path_in} folder"
    exit 1
fi

# ----------------------------------------------
# Handling reference genomes in PRE and REF modes

if [ "${mode_pangen}" != "${name_mode_msa}" ]; then

    # Chech that the reference-genome folder exists
    if [ ! -d "${path_ref}" ]; then
        pokaz_error "Error: Directory ${path_ref} does not exist."
        exit 1
    fi

    # Check that the reference file exists
    file_found=false
    for ext in "${genome_extensions[@]}"; do
        if [ -f "${path_ref}${ref_name}.${ext}" ]; then
            file_found=true
            break
        fi
    done

    if [ "$file_found" = false ]; then
        pokaz_error "Error: No reference genome ${ref_name} was found in ${path_ref}."
        exit 1  
    fi

    refs_all=("${ref_name}")
    
fi

# ----------------------------------------------
# Handling reference genomes in MSA modes
if [ "${mode_pangen}" = "${name_mode_msa}" ]; then

    # ----------------------------------------------
    # Check the number of reference genomes for randomisation

    if [ -z "${ref_num}" ]; then
        ref_num=2
    fi

    if ! [[ "$ref_num" =~ ^[0-9]+$ ]]; then
        with_level 0 pokaz_error "Error: Number of reference genomes is not a number."
        exit 1
    fi

    pokaz_message "Number of genomes for randomisation: ${ref_num}"

    # ----------------------------------------------
    # Check if ref_num is greater than the number of genomes in acc_set
    if (( ref_num > ${#acc_set[@]} )); then
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
            with_level 1 pokaz_attention "Not enough reference genomes. Adding $genomes_needed genome(s)."

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

    fi

    # Print the selected reference genomes for verification
    pokaz_message "Names of genomes for randomisation: $(IFS=,; echo "${refs_all[*]}")"

fi

# ----------------------------------------------
# Number of chromosomes

if [ -z "${nchr}" ] && [ -z "${nchr_ref}" ]; then  # Both nchr and nchr_ref are not defined.
    # Try to define options for number of chromosomes
    pokaz_stage "Define the number of chromosomes..."

    if [ "${mode_pangen}" == "${name_mode_pre}" ]; then  # PRE mode
        nchr=0
        nchr_ref=0
    else   # REF and MSA mode
        # Define from files
        # Count the number of chromosomes in files
        counts=()
        pokaz_message "Acc ChrNum"
        for ext in "${genome_extensions[@]}"; do
            # Loop over each element in acc_set
            for acc in "${acc_set[@]}"; do
                # Find the file with the current acc and extension
                file="$path_in/${acc}.${ext}"
                
                # Check if the file exists
                if [[ -f "$file" ]]; then
                    count=$(grep -c '^>' "$file")
                    pokaz_message "${acc}  ${count}"
                    counts+=("$count")
                fi
            done
        done

        # If the number of chromosomes in genomes vary, raise an error
        if [[ $(echo "${counts[@]}" | tr ' ' '\n' | uniq | wc -l) -eq 1 ]]; then
            nchr=${counts[0]}
            nchr_ref=${nchr}
        else
            pokaz_error "Genomes have different number of chromosomes. Please change files or specify the number -nchr."
            exit 1
        fi
    fi

elif [ -z "${nchr}" ] && [ ! -z "${nchr_ref}" ]; then  # nchr_ref is defined.
    pokaz_error "Error: -nchr not defined, when -nchr_ref is defined ."
    exit 1
elif [ ! -z "${nchr}" ] && [ -z "${nchr_ref}" ]; then  # nchr is defined.
    nchr_ref=${nchr}
fi

# # Check the number of chromosomes in the reference genome
# if [[ ${nchr_ref} -ne 0 ]]; then
#     for ref_name in "${refs_all[@]}"; do
#         for ext in "${genome_extensions[@]}"; do
#             # Find the file with the current acc and extension
#             file="${path_ref}/${ref_name}.${ext}"
#             # Check if the file exists
#             if [[ -f "$file" ]]; then
#                 count=$(grep -c '^>' "$file")
#                 if [[ ${count} -lt ${nchr_ref} ]]; then
#                     pokaz_error "Error: Number of chromosomes in the reference genome ${ref_name} is ${count}. Must be at least ${nchr_ref}."
#                     exit 1
#                 fi
#             fi
#         done
#     done
# fi

# # Check the number of chromosomes in the accessions genome
# if [[ ${nchr} -ne 0 ]]; then
#     for acc_name in "${acc_set[@]}"; do
#         for ext in "${genome_extensions[@]}"; do
#             # Find the file with the current acc and extension
#             file="${path_in}/${acc_name}.${ext}"
            
#             # Check if the file exists
#             if [[ -f "$file" ]]; then
#                 count=$(grep -c '^>' "$file")
#                 if [[ ${count} -lt ${nchr} ]]; then
#                     pokaz_error "Error: Number of chromosomes in the reference genome ${acc_name} is ${count}. Must be at least ${nchr}."
#                     exit 1
#                 fi
#             fi
#         done
#     done
# fi



# ----------------------------------------------
# File with combinations
file_combinations="${path_inter}combinations.txt"

# Free up this file
> ${file_combinations}

if [ -z "${comb_file}" ]; then

    if [[ "${one2one}" == 'T' ]]; then
        if [ "$nchr" != "$nchr_ref" ]; then
            pokaz_error "Error: ${nchr} and ${nchr_ref} should be equal"
            exit 1
        fi

        # Write combinations
        if [ "$nchr" -ne 0 ]; then
            # Write combinations
            for N in $(seq 1 $nchr); do
                echo -e "${N}\t${N}" >> ${file_combinations}
            done
        fi

    else

        # Write combinations
        if [ "$nchr" -ne 0 ] && [ "$nchr_ref" -ne 0 ]; then
            for N in $(seq 1 $nchr); do
                for M in $(seq 1 $nchr_ref); do
                    echo -e "${N}\t${M}" >> ${file_combinations}
                done
            done
        fi
    fi

    option_combinations=" --combinations ${file_combinations}"
else
    cp ${comb_file}  ${file_combinations}
    option_combinations=" --combinations ${file_combinations}"
fi

# ----------------------------------------------
# File with accessions

file_accessions="${path_inter}accessions.txt"
# pokaz_message "Accession file: ${file_accessions}"
printf "%s\n" "${acc_set[@]}" > "${file_accessions}"

# Add reference is it's in path_in, but not in accessions file.
if [[ "$path_ref" == "$path_in" ]]; then
    for s in "${ref_names[@]}"; do
        if [[ ! " ${ass_set[@]} " =~ " $s " ]]; then
            printf "%s\n" "$s" >> "${file_accessions}"
        fi
    done
fi

# ----------------------------------------------

# ----------------------------------------------
# Handling parts

# Default parameters
p_ident="${p_ident:-85}"  
p_ident_gap="${p_ident_gap:-85}"  
part_len="${part_len:-5000}"  
max_len_gap="${max_len_gap:-25000}"  

# Filter repeats

if [ "${purge_reps}" == "F" ]; then
    option_purge_reps=" "
elif [ "${purge_reps}" == "T" ]; then
    option_purge_reps=" --purge.reps T"
else
    pokaz_error "Error: purge_reps must be either 'F' or 'T'."
    exit 1
fi

if [ -z "${flag_rev}" ]; then
    option_mirror=" "
else
    path_parts="${path_inter}parts_mirror/"
    option_mirror=" --purge.reps T"
fi

# Reverse parts
if [ -z "${rev}" ]; then
    option_rev=" "
else
    option_rev=" --rev T "
fi


# Number of Cores

cores="${cores:-1}"  # Number of cores


# ----------------------------------------------------------------------------
#           LOGS
# ----------------------------------------------------------------------------

log_level=${log_level:-1}  # Set the default value to 'steps'

if ! [[ "$log_level" =~ ^[0-3]$ ]]; then
    pokaz_error "Error: log_level must be a number between 0 and 3."
    exit 1
fi

# Hidden path with logs
path_log="${path_out}logs/"
mkdir -p ${path_log}

# File with steps logs
file_log="${path_log}steps.log"
> "${file_log}"


# ----------------------------------------------
# Function for levels

with_level() {
    local log_level_command=$1
    local pokaz_command=$2
    shift 2
    local message="$*"

    # Print to the console
    if [ "${log_level_command}" -le "${log_level}" ]; then
        $pokaz_command "$message"
    fi

    $pokaz_command "$message" | sed 's/\x1b\[[0-9;]*m//g' >> "$file_log"
}


# ----------------------------------------------------------------------------
#            MESSAGES
# ----------------------------------------------------------------------------

if [ "${mode_pangen}" == "${name_mode_pre}" ]; then  # PRE mode
    message_mode="Pannagram runs in the PRELIMINARY mode."
elif [ "${mode_pangen}" == "${name_mode_ref}" ]; then  # PRE mode
    message_mode="Pannagram runs in the REFENRECE-based mode."
elif [ "${mode_pangen}" == "${name_mode_msa}" ]; then  # PRE mode
    message_mode="Pannagram runs in the Multiple Genome Alignment mode."
else 
    with_level 1 pokaz_error "Error: Wrong running mode"
    exit 1
fi

with_level 1 pokaz_attention ${message_mode}

if [ "${mode_pangen}" == "${name_mode_msa}" ]; then  # PRE mode
    with_level 1 pokaz_attention "Path with consensus MSA: ${path_cons}"
fi

with_level 2 pokaz_message "Number of chromosomes ${nchr}"
with_level 2 pokaz_message "Number of cores ${cores}"


# ----------------------------------------------------------------------------
#            Check previous command
# ----------------------------------------------------------------------------

file_params="${path_log}command.log"

if [[ -f "$file_params" ]]; then

    # Loading previous parameters
    source "$file_params"

    if [[ "${prev_mode_pangen}" != "${mode_pangen}" ]]; then
        rm -f ${path_log}/*_done
    fi

    if [[ "$prev_path_in" != "$path_in" || \
          "$prev_part_len" != "$part_len" || \
          "$prev_p_ident_gap" != "$p_ident_gap" || \
          "$prev_purge_reps" != "$purge_reps" || \
          "$prev_p_ident" != "$p_ident"  ]]; then
        pokaz_error "Error: One or more parameters have been changed!"
        
        [[ "$prev_path_in" != "$path_in" ]] && pokaz_attention "path_in: $prev_path_in -> $path_in"
        [[ "$prev_part_len" != "$part_len" ]] && pokaz_attention "part_len: $prev_part_len -> $part_len"
        [[ "$prev_p_ident_gap" != "$p_ident_gap" ]] && pokaz_attention "p_ident_gap: $prev_p_ident_gap -> $p_ident_gap"
        [[ "$prev_purge_reps" != "$purge_reps" ]] && pokaz_attention "purge_reps: $prev_purge_reps -> $purge_reps"
        [[ "$prev_p_ident" != "$p_ident" ]] && pokaz_attention "p_ident: $prev_p_ident -> $p_ident"
        
        exit 1
    fi
fi

echo "prev_mode_pangen=${mode_pangen}" > "$file_params"
echo "prev_path_in=${path_in}" >> "$file_params"
echo "prev_part_len=${part_len}" >> "$file_params"
echo "prev_purge_reps=${purge_reps}" >> "$file_params"
echo "prev_p_ident=${p_ident}" >> "$file_params"
echo "prev_p_ident_gap=${p_ident_gap}" >> "$file_params"


# ----------------------------------------------------------------------------
#            Check Steps
# ----------------------------------------------------------------------------

# Define the log directory
step_files=$(find "${path_log}" -type f -name "step*_done")

if [[ $step_start -eq 0 ]]; then

    # If there is at least one matching file
    if [[ -n "$step_files" ]]; then
        
        max_step=0  # Initialize to 0 
        for file in $step_files; do
            # Extract the number immediately following "step"
            step_num=$(echo "$file" | sed -n 's/.*step\([0-9]\{1,\}\)_.*/\1/p')

            # Check if the extracted number is the highest found so far
            if [[ $step_num -gt $max_step ]]; then
                max_step=$step_num
            fi
        done

        # Assign the maximum found value to step_start
        step_start=$((max_step + 1))
    else
        step_start=1
    fi

else
    # Remove all done files between start and end
    if [[ -n "$step_files" ]]; then
        # If the previous tep does not exist - error
        step_prev=$((step_start - 1))
        if [[ $step_prev -ne 0 ]] && ! ls "${path_log}step${step_prev}"*done 1> /dev/null 2>&1; then
            pokaz_error "Error: Step ${step_prev} was not done."
            exit 1
        fi

        for file in $step_files; do
            
            step_num=$(echo "$file" | sed -n 's/.*step\([0-9]\{1,\}\)_.*/\1/p')

            if [[ $step_num -ge $step_start && $step_num -le $step_end ]]; then
                # echo ${file}
                rm ${file}
            fi
        done

    fi 
fi

if [ "$step_end" -ne 100 ]; then
    pokaz_attention "Attention: the pipeline runs until step ${step_end}";
fi

if [ "$one_step" == "T" ]; then
    pokaz_attention "Attention: running only one step!";
    step_end=${step_start}
fi

pokaz_message "Start End: ${step_start} ${step_end}"

# Check start and end
if [ "$step_start" -gt "$step_end" ]; then
    echo "Error: step_start ($step_start) is greater than step_end ($step_end)"
    exit 1  # Exit the script with error status 1
fi


# ╔═══════════════════════════════════════════════════════════════════════════════════╗  /\_/\  
# ║                                                                                   ║ ( o.o ) 
# ║                                  MAIN PIPELINE                                    ║  > ^ <
# ║                                                                                   ║ /     \
# ╚═══════════════════════════════════════════════════════════════════════════════════╝ \_____/

# ┌───────────────────────────────────────────────────────────────────────────┐
# │                         CHROMOSOMES & PARTS                               │
# └───────────────────────────────────────────────────────────────────────────┘
    
step_num=1

# ----------------------------------------------
# Split query fasta into chromosomes

with_level 1 pokaz_stage "Step ${step_num}. Genomes into chromosomes."

# Logs
step_name="step${step_num}_query_01"
step_file="${path_log}${step_name}_done"
path_log_step="${path_log}${step_name}/"
make_dir ${path_log_step}

# Start
if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

    # Clean up the output folders
    if   [ "$clean" == "T" ]; then 
        touch ${path_chrom}fake.fasta
        touch ${path_chrom}fake.txt
        touch ${path_log_step}fake.log

        rm -f ${path_chrom}*fasta
        rm -f ${path_chrom}*txt
        rm -f ${path_log_step}*
    fi

    # Run the step
    Rscript $INSTALLED_PATH/pangen/query_01_to_chr.R --path.in ${path_in} --path.out ${path_chrom} \
            --cores ${cores}  \
            --n.chr ${nchr}  \
            --accessions ${file_accessions} \
            --path.log ${path_log_step} \
            --log.level ${log_level}

    # Done
    touch "${step_file}"
    
fi

source $INSTALLED_PATH/utils/chunk_step_done.sh

# ----------------------------------------------
# Split reference fasta into chromosomes if additionally needed

if [[ "${path_in}" != "$path_ref" ]]; then

    ((step_num = step_num - 1))

    with_level 1 pokaz_stage "Additional step. Reference genome into chromosomes."  # IT SHOULD BE INSIDE THE "IF"
    for ref0 in "${refs_all[@]}"; do

        # Logs
        step_name="step${step_num}_query_01_refpart_${ref0}/"
        step_file="${path_log}${step_name}_done"
        path_log_step="${path_log}${step_name}/"
        make_dir ${path_log_step}

        # Start
        if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

            with_level 1  pokaz_attention "Reference ${ref0}"

            # Temporary file to analyse only the reference genome from the folder
            file_acc_ref=${path_cons}ref_acc.txt
            echo "${ref0}" > ${file_acc_ref}

            # Run the step
            Rscript $INSTALLED_PATH/pangen/query_01_to_chr.R \
                    --path.in ${path_ref} --path.out ${path_chrom}   \
                    --cores ${cores} \
                    --accessions ${file_acc_ref} \
                    --n.chr ${nchr_ref}  \
                    --path.log ${path_log_step} --log.level ${log_level}

            # Remove the temporary file
            rm ${file_acc_ref}

            # Done
            touch "${step_file}"
        fi
    done

    source $INSTALLED_PATH/utils/chunk_step_done.sh
fi

# ----------------------------------------------
# ORF-Finder

if [ ! -z "${flag_orf}" ]; then

    ((step_num = step_num - 1))

    with_level 1 pokaz_stage "Additional step. Get all ORFs."

    # Logs
    step_name="step${step_num}_query_01_orf"
    step_file="${path_log}${step_name}_done"
    path_log_step="${path_log}${step_name}/"
    make_dir ${path_log_step}

    # Start
    if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

        # Make ORF folder
        path_orf="${path_inter}orf/"
        mkdir -p "${path_orf}"

        # Clean up the output folders
        if   [ "$clean" == "T" ]; then 
            touch ${path_orf}fake.fasta
            touch ${path_log_step}fake.log

            rm -f ${path_orf}*fasta
            rm -f ${path_log_step}*
        fi

        # Run the step
        Rscript $INSTALLED_PATH/pangen/query_01_to_orf.R  \
                --path.chr ${path_chrom}  \
                --path.orf ${path_orf} \
                --cores ${cores}  \
                --n.chr ${nchr}  \
                --accessions ${file_accessions} \
                --path.log ${path_log_step} --log.level ${log_level}

        # Done
        touch "${step_file}"
    fi

    source $INSTALLED_PATH/utils/chunk_step_done.sh

fi

# ----------------------------------------------
# Split query chromosomes into parts

with_level 1 pokaz_stage "Step ${step_num}. Chromosomes into parts."

# Logs
step_name="step${step_num}_query_02"
step_file="${path_log}${step_name}_done"
path_log_step="${path_log}${step_name}/"
make_dir ${path_log_step}
make_dir ${path_parts}

# Start
if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

    # Clean up the output folders
    if   [ "$clean" == "T" ]; then 
        touch ${path_parts}fake.fasta
        touch ${path_log_step}fake.log

        rm -f ${path_parts}*fasta
        rm -f ${path_log_step}*
    fi

    # Run the step
    Rscript $INSTALLED_PATH/pangen/query_02_to_parts.R --path.chr ${path_chrom} \
            --path.parts ${path_parts} \
            --part.len $part_len \
            --cores ${cores} \
            --accessions ${file_accessions} \
            --n.chr ${nchr}  \
            ${option_purge_reps}  \
            ${option_mirror}  ${option_rev}  \
            --path.log ${path_log_step} --log.level ${log_level}

    # Done
    touch "${step_file}"

fi

source $INSTALLED_PATH/utils/chunk_step_done.sh

# ┌───────────────────────────────────────────────────────────────────────────┐
# │                         REFERENCE-BASED                                   │
# └───────────────────────────────────────────────────────────────────────────┘

# ----------------------------------------------
# Blast parts on the reference genome

with_level 1 pokaz_stage "Step ${step_num}. BLAST of parts against the reference genome."
for ref0 in "${refs_all[@]}"; do

    # Paths
    path_blast_parts=${path_inter}blast_parts_${ref0}/
    make_dir ${path_blast_parts}
    
    # Logs
    step_name="step${step_num}_query_03_blast_${ref0}"
    step_file="${path_log}${step_name}_done"
    path_log_step="${path_log}${step_name}/"
    make_dir ${path_log_step}

    # Start
    if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

        with_level 1  pokaz_attention "Reference ${ref0}"

        # with_level 1 pokaz_message "NOTE: if this stage takes relatively long, use -purge_repeats -s 2 to mask highly repetative regions"

        if [[ $option_purge_reps =~ ^[[:space:]]*$ ]]; then
            with_level 1 pokaz_message "NOTE: if this stage takes relatively long, use -purge_repeats -s $((${step_num}-1)) to mask highly repetative regions"
        fi
    
        # ---- Clean up the output folders ----
        if   [ "$clean" == "T" ]; then 
            touch ${path_blast_parts}fake.txt
            touch ${path_log_step}fake.log

            rm -f ${path_blast_parts}*txt
            rm -f ${path_log_step}*
        fi    

        # ---- Create a database on the reference genome ----
        with_level 2 pokaz_message "Create BLAST databases."

        # Check if there are any files matching the pattern
        shopt -s nullglob
        files=(${path_chrom}${ref0}_chr*fasta)
        shopt -u nullglob

        if [ ${#files[@]} -eq 0 ]; then
          with_level 0 pokaz_error "Error: No files matching the pattern ${path_chrom}${ref0}_chr*fasta"
          exit 1
        fi

        # Iterate over each file and create the database
        for file in "${files[@]}"; do
          # Check if the BLAST database files already exist
          if [ ! -f "${file}.nin" ]; then
              makeblastdb -in ${file} -dbtype nucl > /dev/null
          fi
        done

        # ---- Run BLAST ----

        with_level 2 pokaz_message "Run BLAST."    

        # Blast parts on the reference genome
        $INSTALLED_PATH/pangen/query_03_blast_parts.sh \
                -path_chrom ${path_chrom} \
                -path_parts ${path_parts} \
                -path_result ${path_blast_parts} \
                -combinations ${file_combinations} \
                -accessions ${file_accessions} \
                -ref ${ref0} \
                -p_ident ${p_ident} \
                -cores ${cores} \
                -path_log ${path_log_step}

        # Done
        touch "${step_file}"
        
    fi

    # Remove reference-dependent variables
    unset path_blast_parts

done

source $INSTALLED_PATH/utils/chunk_step_done.sh

# ----------------------------------------------
# Search for the mirror universe

if [ ! -z "${flag_rev}" ]; then
    echo "Pannagram's Mirror mode is done."
    exit
fi

# ----------------------------------------------
# First round of alignments

with_level 1 pokaz_stage "Step ${step_num}. Alignment-1: Remaining syntenic (major) matches."
for ref0 in "${refs_all[@]}"; do
    
    # Paths
    path_blast_parts=${path_inter}blast_parts_${ref0}/
    path_alignment=${path_inter}alignments_${ref0}/
    make_dir ${path_alignment}

    # Logs
    step_name="step${step_num}_synteny_01_maj_${ref0}"
    step_file="${path_log}${step_name}_done"
    path_log_step="${path_log}${step_name}/"
    make_dir ${path_log_step}

    # Start
    if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

        with_level 1  pokaz_attention "Reference ${ref0}"

        # Clean up the output folders
        if   [ "$clean" == "T" ]; then 
            touch ${path_alignment}fake.rds
            touch ${path_log_step}fake.log

            rm -f ${path_alignment}*rds
            rm -f ${path_log_step}*
        fi    
        
        # Run the step
        Rscript $INSTALLED_PATH/pangen/synteny_01_majoir.R  \
                --path.blast ${path_blast_parts} \
                --path.aln ${path_alignment} \
                --ref ${ref0}   \
                --combinations ${file_combinations} \
                --accessions ${file_accessions} \
                --cores ${cores} \
                --path.log ${path_log_step} \
                --log.level ${log_level}

        # Done
        touch "${step_file}"

        # Clean up the output folders of previous stages
        # If the first round of alignment didn't have any errors - remove the blast which was needed for it
        # rm ${path_blast_parts}
    fi

    # Remove reference-dependent variables
    unset path_blast_parts
    unset path_alignment

done

source $INSTALLED_PATH/utils/chunk_step_done.sh

# ----------------------------------------------
# Plotting

with_level 1 pokaz_stage "Step ${step_num}. Plotting the results."
for ref0 in "${refs_all[@]}"; do

    # Paths
    path_alignment=${path_inter}alignments_${ref0}/
    path_plots_ref=${path_plots}plots_${ref0}/
    make_dir ${path_plots_ref}
    
    # Logs
    step_name="step${step_num}_synteny_02_plot_${ref0}"
    step_file="${path_log}${step_name}_done"
    path_log_step="${path_log}${step_name}/"
    make_dir ${path_log_step}

    # Step start
    if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

        with_level 1  pokaz_attention "Reference ${ref0}"

        # ---- Clean up the output folders ----
        if [ "$clean" == "T" ]; then 
            touch ${path_plots_ref}fake.pdf
            touch ${path_log_step}fake.log

            rm -f ${path_plots_ref}*pdf
            rm -f ${path_log_step}*
        fi    

        # Run the step
        Rscript $INSTALLED_PATH/pangen/synteny_02_plot.R \
                --ref ${ref0} \
                --path.chr ${path_chrom} \
                --path.plot ${path_plots_ref} \
                --path.aln ${path_alignment} \
                --accessions ${file_accessions} \
                --path.log ${path_log_step} \
                --log.level ${log_level} \
                --cores ${cores}

        # Done
        touch "${step_file}"
        
    fi

    # Remove reference-dependent variables
    unset path_alignment
    unset path_plots_ref

done

source $INSTALLED_PATH/utils/chunk_step_done.sh  

# ========== CHECK MODE ==========

if [ "${mode_pangen}" == "${name_mode_pre}" ]; then 
    with_level 0 pokaz_message "Pannagram's PREliminary mode is done."
    exit
fi 

# ┌───────────────────────────────────────────────────────────────────────────┐
# │                         FILL THE GAPS                                     │
# └───────────────────────────────────────────────────────────────────────────┘

# ----------------------------------------------
# Get sequences between the synteny blocks

with_level 1 pokaz_stage "Step ${step_num}. Get gaps between syntenic matches."
for ref0 in "${refs_all[@]}"; do        

    # Paths
    path_alignment=${path_inter}alignments_${ref0}/
    path_gaps=${path_inter}blast_gaps_${ref0}/
    make_dir ${path_gaps}

    # Logs
    step_name="step${step_num}_synteny_03_get_gaps_${ref0}"
    step_file="${path_log}${step_name}_done"
    path_log_step="${path_log}${step_name}/"
    make_dir ${path_log_step}

    # Step start
    if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

        with_level 1  pokaz_attention "Reference ${ref0}"

        # ---- Clean up the output folders ----
        if  [ "$clean" == "T" ]; then 
            touch ${path_gaps}fake.fasta
            touch ${path_log_step}fake.log

            rm -f ${path_gaps}*fasta
            rm -f ${path_log_step}*
        fi  

        # Run
        Rscript $INSTALLED_PATH/pangen/synteny_03_get_gaps.R \
                --path.chr ${path_chrom} \
                --path.aln ${path_alignment} \
                --path.gaps  ${path_gaps} \
                --ref ${ref0}   \
                --accessions ${file_accessions} \
                --combinations ${file_combinations} \
                --cores ${cores} \
                --path.log ${path_log_step} --log.level ${log_level}

        # Done
        touch "${step_file}"
    fi
    
    unset path_alignment
    unset path_gaps
done

source $INSTALLED_PATH/utils/chunk_step_done.sh

# ----------------------------------------------
# Blast regions between synteny blocks

with_level 1 pokaz_stage "Step ${step_num}. BLAST of gaps between syntenic matches."
for ref0 in "${refs_all[@]}"; do  

    # Paths
    path_alignment=${path_inter}alignments_${ref0}/
    path_gaps=${path_inter}blast_gaps_${ref0}/
    make_dir "${path_gaps}db/"

    # Logs
    step_name="step${step_num}_synteny_04_blast_gaps_${ref0}"
    step_file="${path_log}${step_name}_done"
    path_log_step="${path_log}${step_name}/"
    make_dir ${path_log_step}

    # Step start
    if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

        with_level 1  pokaz_attention "Reference ${ref0}"

        # ---- Clean up the output folders ----
        if 
            [ "$clean" == "T" ]; then 
            touch ${path_gaps}db/fake.db
            touch ${path_gaps}fake_out.txt
            touch ${path_log_step}fake.log

            rm -f ${path_gaps}db/*
            rm -f ${path_gaps}*out.txt
            rm -f ${path_log_step}*
        fi  

        # Run BLAST for gaps
        $INSTALLED_PATH/pangen/synteny_04_blast_gaps.sh \
                -path_gaps ${path_gaps} \
                -cores ${cores} \
                -log_path ${path_log_step} \
                -p_ident ${p_ident_gap}

        # Done
        touch "${step_file}"
    fi

    unset path_alignment
    unset path_gaps
done

source $INSTALLED_PATH/utils/chunk_step_done.sh

# ----------------------------------------------
# Second round of alignments

with_level 1 pokaz_stage "Step ${step_num}. Alignment-2: Fill the gaps between synteny blocks."
for ref0 in "${refs_all[@]}"; do  
    
    # Paths
    path_alignment=${path_inter}alignments_${ref0}/
    path_gaps=${path_inter}blast_gaps_${ref0}/

    # Logs
    step_name="step${step_num}_synteny_05_full_${ref0}"
    step_file="${path_log}${step_name}_done"
    path_log_step="${path_log}${step_name}/"
    make_dir ${path_log_step}

    # Step start
    if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

        with_level 1  pokaz_attention "Reference ${ref0}"

        # ---- Clean up the output folders ----
        if [ "$clean" == "T" ]; then 
            touch ${path_alignment}fake_full.rds
            touch ${path_log_step}fake.log

            rm -f ${path_alignment}*full.rds
            rm -f ${path_log_step}*
        fi  

        # Run the step
        Rscript $INSTALLED_PATH/pangen/synteny_05_merge_gaps.R  \
                --path.chr ${path_chrom}\
                --path.aln ${path_alignment} \
                --path.gaps ${path_gaps}   \
                --ref ${ref0} \
                --accessions ${file_accessions} \
                --combinations ${file_combinations} \
                --cores ${cores} \
                --path.log ${path_log_step} \
                --log.level ${log_level}

        # # If the second round of alignment didn't have any errors - remove the blast which was needed for it
        # rm -rf ${path_gaps}
        # ls ${path_alignment}*maj*
        # rm -rf ${path_alignment}*maj*

        # Done
        touch "${step_file}"
    fi

    unset path_alignment
    unset path_gaps
done

source $INSTALLED_PATH/utils/chunk_step_done.sh

# ----------------------------------------------
# Create a consensus

with_level 1 pokaz_stage "Step ${step_num}. Combine reference-based alignments by chromosomes."
for ref0 in "${refs_all[@]}"; do  

    # Paths
    path_alignment=${path_inter}alignments_${ref0}/

    # Logs
    step_name="step${step_num}_comb_01_${ref0}"
    step_file="${path_log}${step_name}_done"
    path_log_step="${path_log}${step_name}/"
    make_dir ${path_log_step}

    # Step start
    if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

        with_level 1  pokaz_attention "Reference ${ref0}"

        # ---- Clean up the output folders ----
        if   [ "$clean" == "T" ]; then 
            touch ${path_cons}ref_fake_${ref0}_.h5
            touch ${path_log_step}fake.log

            rm -f ${path_cons}ref*${ref0}_.h5
            rm -f ${path_log_step}*
        fi  

        # Run the step
        Rscript $INSTALLED_PATH/pangen/comb_01_one_ref.R \
                --path.cons ${path_cons} \
                --path.aln ${path_alignment} \
                --path.chr ${path_chrom} \
                --ref ${ref0}  \
                --accessions ${file_accessions} \
                --combinations ${file_combinations} \
                --cores ${cores} \
                --path.log ${path_log_step} \
                --log.level ${log_level}

        # Done
        touch "${step_file}"

    fi
    unset path_alignment

done

source $INSTALLED_PATH/utils/chunk_step_done.sh

# ----------------------------------------------------------------------------
#            CHECK MODE
# ----------------------------------------------------------------------------
if [ "${mode_pangen}" != "${name_mode_msa}" ]; then 
    echo "Pannagram's REFRENCE-based mode is done."
    exit
fi

# ./pannagram.sh -path_in /Volumes/Samsung_T5/vienn/genomes/symA -path_out /Volumes/Samsung_T5/vienn/test/symA_test_0 -log 2 -cores 2

# ╔═══════════════════════════════════════════════════════════════════════════════════╗       ______
# ║                                                                                   ║      /|_||_\`.__
# ║                                  MSA                                              ║     (   _    _ _\
# ║                                                                                   ║     =`-(_)--(_)-'
# ╚═══════════════════════════════════════════════════════════════════════════════════╝ 

with_level 1 pokaz_stage "MGA is strting..  (^>^)"

# ----------------------------------------------
# Run consensus for a pair of files

with_level 1 pokaz_stage "Step ${step_num}. Randomisation of references.." 

ref0=${refs_all[0]}
for ((i = 1; i < ${#refs_all[@]}; i++)); do
    ref1=${refs_all[i]}

    # Logs
    step_name="step${step_num}_comb_02_${ref0}_${ref1}"
    step_file="${path_log}${step_name}_done"
    path_log_step="${path_log}${step_name}/"
    make_dir ${path_log_step}
        
    # Start
    if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

        with_level 1 pokaz_attention "Combine two references: ${ref0} and ${ref1}.."  # SHOULD BE INSIDE THE LOOP
        
        # ---- Clean up the output folders ----
        if   [ "$clean" == "T" ]; then 
            touch ${path_cons}comb_fake.h5
            touch ${path_log_step}fake.log

            rm -f ${path_cons}comb*h5
            rm -f ${path_log_step}*
        fi  

        Rscript $INSTALLED_PATH/pangen/comb_02_two_refs.R \
                --path.cons ${path_cons} \
                --ref0 ${ref0} \
                --ref1 ${ref1} \
                --cores ${cores} \
                --path.log ${path_log_step} \
                --log.level ${log_level}

        touch "${step_file}"
    fi

done

source $INSTALLED_PATH/utils/chunk_step_done.sh

# # ----------------------------------------------
# Remain only the trustable positions
with_level 1 pokaz_stage "Step ${step_num}. Remain only the trustable syntenic positions.."

# Logs
step_name="step${step_num}_comb_03_cleanup"
step_file="${path_log}${step_name}_done"
path_log_step="${path_log}${step_name}/"
make_dir ${path_log_step}

# Start
if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

    # ---- Clean up the output folders ----
    if   [ "$clean" == "T" ]; then 
        touch ${path_cons}clean_fake.h5
        touch ${path_log_step}fake.log

        rm -f ${path_cons}clean*h5
        rm -f ${path_log_step}*
    fi  

    Rscript $INSTALLED_PATH/pangen/comb_03_cleanup.R \
            --path.cons ${path_cons} \
            --cores ${cores} \
            --path.log ${path_log_step} \
            --log.level ${log_level} 

    # Done
    touch "${step_file}"
fi

source $INSTALLED_PATH/utils/chunk_step_done.sh

# ----------------------------------------------
# Prepare breaks

with_level 1 pokaz_stage "Step ${step_num}. Prepare breakes for an additional alignment"

# Logs
step_name="step${step_num}_comb_04_prepare_breaks"
step_file="${path_log}${step_name}_done"
path_log_step="${path_log}${step_name}/"
make_dir ${path_log_step}

# Start
if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

    # ---- Clean up the output folders ----
    if [ "$clean" == "T" ]; then 
        touch ${path_cons}breaks_ws_fake.RData
        touch ${path_log_step}fake.log

        rm -f  ${path_cons}breaks_ws_*.RData
        rm -f ${path_log_step}*
    fi  

    Rscript $INSTALLED_PATH/pangen/comb_04_prepare_breaks.R \
            --path.cons "${path_cons}" \
            --cores "${cores}" \
            --path.chromosomes "${path_chrom}" \
            --path.log "${path_log_step}" \
            --log.level "${log_level}" \
            --max.len.gap "${max_len_gap}"

    # Done
    touch "${step_file}"
fi

source $INSTALLED_PATH/utils/chunk_step_done.sh
# ----------------------------------------------
# Create sequences

with_level 1 pokaz_stage "Step ${step_num}. Prepare sequences for alignments."

# Logs
step_name="step${step_num}_comb_04_prepare_seqs"
step_file="${path_log}${step_name}_done"
path_log_step="${path_log}${step_name}/"
make_dir ${path_log_step}

# Paths for MAFFT, common for the next code too
path_mafft_in="${path_inter}mafft_in/"
path_mafft_out="${path_inter}mafft_out/"
if [ ! -d "$path_mafft_in" ]; then
    mkdir -p "$path_mafft_in"
fi
if [ ! -d "$path_mafft_out" ]; then
    mkdir -p "$path_mafft_out"
fi

# Start
if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

    # ---- Clean up the output folders ----
    if [ "$clean" == "T" ]; then 
        touch ${path_mafft_in}fake.fasta
        touch ${path_mafft_in}fake.tree
        touch ${path_log_step}fake.log
        touch ${path_cons}small_ws_fake.RData

        find ${path_mafft_in} -name "*.fasta" -type f -exec rm -f {} +
        find ${path_mafft_in} -name "*.tree" -type f -exec rm -f {} +
        # rm -f ${path_mafft_in}*fasta
        rm -f ${path_log_step}*
        rm -f ${path_cons}small_ws_*.RData
    fi  

    Rscript $INSTALLED_PATH/pangen/comb_04_prepare_seqs.R \
            --path.cons "${path_cons}" \
            --cores "${cores}" \
            --path.chromosomes "${path_chrom}" \
            --path.mafft.in "${path_mafft_in}" \
            --path.log "${path_log_step}" \
            --log.level "${log_level}" \
            --max.len.gap "${max_len_gap}"

    # Done
    touch "${step_file}"
fi

source $INSTALLED_PATH/utils/chunk_step_done.sh

# ----------------------------------------------
# Perform some small alignments

with_level 1 pokaz_stage "Step ${step_num}. Align short sequences."

# Logs
step_name="step${step_num}_comb_05_small"
step_file="${path_log}${step_name}_done"
path_log_step="${path_log}${step_name}/"
make_dir ${path_log_step}

# Start
if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

    # ---- Clean up the output folders ----
    if [ "$clean" == "T" ]; then 
        touch ${path_log_step}fake.log
        rm -f ${path_log_step}*
    fi  

    Rscript $INSTALLED_PATH/pangen/comb_05_small.R \
            --path.cons "${path_cons}" \
            --cores "${cores}" \
            --path.log "${path_log_step}" \
            --log.level "${log_level}" \
            --max.len.gap "${max_len_gap}"

    # Done
    touch "${step_file}"
fi

source $INSTALLED_PATH/utils/chunk_step_done.sh

# ----------------------------------------------
# Run MAFFT

with_level 1 pokaz_stage "Step ${step_num}. Run MAFFT."

# Logs
step_name="step${step_num}_comb_05_mafft"
step_file="${path_log}${step_name}_done"
path_log_step="${path_log}${step_name}/"
make_dir "${path_log_step}"

# Start
if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

    # ---- Clean up the output folders ----
    if   [ "$clean" == "T" ]; then 
        touch ${path_mafft_out}fake_aligned.fasta
        touch ${path_log_step}fake.log

        find ${path_mafft_out} -name "*aligned*.fasta" -type f -exec rm -f {} +
        # rm -f ${path_mafft_out}*aligned.fasta
        rm -f ${path_log_step}*
    fi

    "$INSTALLED_PATH/pangen/comb_05_mafft.sh" \
            -cores "${cores}" \
            -path_mafft_in "${path_mafft_in}" \
            -path_mafft_out "${path_mafft_out}" \
            -log_path "${path_log_step}"

    # Done
    touch "${step_file}"
fi


source $INSTALLED_PATH/utils/chunk_step_done.sh

# ----------------------------------------------
# Remove bad Mafft alignments

with_level 1 pokaz_stage "Step ${step_num}. Remove bad mafft."

# Logs
step_name="step${step_num}_comb_06_bad_mafft"
step_file="${path_log}${step_name}_done"
path_log_step="${path_log}${step_name}/"
make_dir ${path_log_step}

# Start
if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

    # ---- Clean up the output folders ----
    if   [ "$clean" == "T" ]; then 
        touch ${path_log_step}fake.log
        rm -f ${path_log_step}*
    fi

    Rscript $INSTALLED_PATH/pangen/comb_06_bad_mafft.R \
            --cores ${cores} \
            --path.mafft.out ${path_mafft_out} \
            --path.log ${path_log_step} \
            --log.level ${log_level}

    # Done
    touch "${step_file}"
fi

source $INSTALLED_PATH/utils/chunk_step_done.sh


# ----------------------------------------------
# Additional MAFFT

with_level 1 pokaz_stage "Step ${step_num}. Run ADDITIONAL MAFFT."

# Logs
step_name="step${step_num}_comb_06_mafft2"
step_file="${path_log}${step_name}_done"
path_log_step="${path_log}${step_name}/"
make_dir ${path_log_step}

# Start
if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

    # ---- Clean up the output folders ----
    if   [ "$clean" == "T" ]; then 
        touch ${path_mafft_out}fake_aligned2.fasta
        touch ${path_log_step}fake.log

        rm -f ${path_mafft_out}*aligned2.fasta
        rm -f ${path_log_step}*
    fi

    Rscript $INSTALLED_PATH/pangen/comb_06_mafft2.R \
            --cores ${cores} \
            --path.mafft.in ${path_mafft_in} \
            --path.mafft.out ${path_mafft_out} \
            --path.log ${path_log_step} \
            --log.level ${log_level}

    # Done
    touch "${step_file}"
fi

source $INSTALLED_PATH/utils/chunk_step_done.sh

# ----------------------------------------------
# Combine all together

with_level 1 pokaz_stage "Step ${step_num}. Combine all alignments together into the final one."

# Logs
step_name="step${step_num}_comb_07"
step_file="${path_log}${step_name}_done"
path_log_step="${path_log}${step_name}/"
make_dir ${path_log_step}

# Start
if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

    # ---- Clean up the output folders ----
    if   [ "$clean" == "T" ]; then 
        touch ${path_cons}msa_fake_h5
        touch ${path_log_step}fake.log

        rm -f ${path_cons}msa*h5
        rm -f ${path_log_step}*
    fi  

    Rscript $INSTALLED_PATH/pangen/comb_07_final_aln.R  \
            --cores ${cores} \
            --path.mafft.out ${path_mafft_out} \
            --path.cons ${path_cons} \
            --path.log ${path_log_step} \
            --log.level ${log_level}

    # Done
    touch "${step_file}"
fi

source $INSTALLED_PATH/utils/chunk_step_done.sh

# ----------------------------------------------
# Add synteny positions which were lost

with_level 1 pokaz_stage "Step ${step_num}. Add synteny positions which were lost."

# Paths
path_extra="${path_inter}extra_regions/"
if [ ! -d "$path_extra" ]; then
    mkdir -p "$path_extra"
fi

# Logs
step_name="step${step_num}_comb_08"
step_file="${path_log}${step_name}_done"
path_log_step="${path_log}${step_name}/"
make_dir ${path_log_step}

# Start
if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

    Rscript $INSTALLED_PATH/pangen/comb_08_insert_back.R  \
            --cores ${cores} \
            --path.cons ${path_cons} \
            --path.log ${path_log_step} \
            --log.level ${log_level}

    # Done
    touch "${step_file}"
fi

source $INSTALLED_PATH/utils/chunk_step_done.sh

# ----------------------------------------------
# Alignment of extra long fragments

# with_level 1 pokaz_stage "Step ${step_num}. Alignment of extra long fragments."

# # Paths
# path_extra="${path_inter}extra_regions/"
# if [ ! -d "$path_extra" ]; then
#     mkdir -p "$path_extra"
# fi

# # Logs
# step_name="step${step_num}_comb_09"
# step_file="${path_log}${step_name}_done"
# path_log_step="${path_log}${step_name}/"
# make_dir ${path_log_step}

# # Start
# if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

#     Rscript $INSTALLED_PATH/pangen/comb_09_extra_seqs.R  \
#             --cores ${cores} \
#             --path.chromosomes "${path_chrom}" \
#             --path.extra ${path_extra} \
#             --path.cons ${path_cons} \
#             --path.log ${path_log_step} \
#             --log.level ${log_level}

#     # Done
#     touch "${step_file}"
# fi

# source $INSTALLED_PATH/utils/chunk_step_done.sh



# # ----------------------------------------------
# # Get synteny blocks

# with_level 1 pokaz_stage "Step ${step_num}. Get synteny blocks."

# # Logs
# step_name="step${step_num}_analys_01_blocks"
# step_file="${path_log}${step_name}_done"
# path_log_step="${path_log}${step_name}/"
# make_dir ${path_log_step}

# # Start
# if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

#     Rscript $INSTALLED_PATH/analys/analys_01_blocks.R \
#             --path.cons ${path_cons} \
#             --cores ${cores}
#     # Done
#     touch "${step_file}"
# fi

# source $INSTALLED_PATH/utils/chunk_step_done.sh

# with_level 1 pokaz_message "* The pipeline is done."

# # if [ $step_start -eq 0 ]; then
# #     rm -f "$FLAG_DIR"/.*
# #     echo "Script completed successfully"
# # fi









