
# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

source utils/chunk_error_control.sh

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source utils/utils_bash.sh
source utils/utils_help.sh

# ----------------------------------------------------------------------------
#            PARAMETERS: parsing
# ----------------------------------------------------------------------------

mode_pre="F"
mode_ref="F"
mode_msa="F"
clean_dir="F"
one_step="F"

unrecognized_options=()

start_step=100
while [ $# -gt 0 ]
do
    # echo $1
    case $1 in
        -h) print_usage_short; print_examples; exit ;;
        -help ) print_usage_detailed; print_examples; exit ;;
        -s | -stage ) start_step="$2"; shift 2 ;;  # stage from which to run, when the stage is not provided - the last interrupted stage withh be re-run
        -log)         log_level=$2;    shift 2 ;;  # path to the output
        -cores)       cores=$2;        shift 2 ;;
        -clean_dir)   clean_dir="T"    shift 1 ;;
        -one_step)    one_step="T"    shift 1 ;;
        
        # Required
        -path_out) path_out=$2; shift 2 ;;  # path to the output
        -path_in)  path_in=$2;  shift 2 ;;  # path with all genomes in fasta format

        # REF-based
        -pre)      mode_pre="T"; shift 1 ;;  # shitch to the preliminary mode
        -ref)      ref_name=$2; mode_ref="T"; shift 2 ;;  # name of the reference genome
        -path_ref) path_ref=$2; shift 2 ;;  # dont provide if it's the same folder as the path with query genomes

        # MSA
        -refs) ref_set="$2"; mode_msa="T"; shift 2 ;;
        -nref) ref_num="$2"; mode_msa="T"; shift 2 ;;

        # # Optional paths
        # -path_chrom) path_chrom=$2;  shift 2 ;;  # path to the folder with individual chromosomes in separate files
        # -path_parts) path_parts=$2;  shift 2 ;;  # path to the folder with chromosomal parts
        # -path_cons)  path_cons=$2;   shift 2 ;;  # path to the consensus folder

        # Number of chromosomes
        -nchr)      nchr=$2;     shift 2 ;;  # number of chromosomes
        -nchr_ref)  nchr_acc=$2; shift 2 ;;  # number of chromosome in the query genome

        -part_len)  part_len=$2; shift 2 ;;  # fragments to which each chromosome should be cut, has a default value 5000
        -p_ident)   p_ident=$2;  shift 2 ;;  # percent of identity

        -sort_by_len) sort_chr_len="T"; shift 1 ;;  # flag whether to sort chromosomes by length or not
        -one2one)     one2one="${one2one}T";      shift 1 ;;  # compare chromosomes one-to-one or not (default in REF and MSA modes)
        -all2all)     one2one="${one2one}F";      shift 1 ;;  # compare chromosomes all-to-all or not (default in PRE mode)
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

if [[ "${one2one}" == 'T' ]]; then
    option_one2one=" -one2one "
    option_one2one_r=" --one2one T "
else
    option_one2one=" "
fi

# ----------------------------------------------------------------------------
#            PARAMETERS: checking
# ----------------------------------------------------------------------------

# ----------------------------------------------
# Required parameters
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

# ----------------------------------------------
# Check the input folder

genome_extensions=('fasta' 'fna' 'fa' 'fas')

# If the input folder exists
if [ -d "${path_in}" ]; then
    found_files=false
    for ext in "${genome_extensions[@]}"; do
        if ls "${path_in}"/*.${ext} > /dev/null 2>&1; then
            found_files=true
            break
        fi
    done

    if [ "$found_files" = false ]; then
        pokaz_error "Error: No files with the specified extensions found in ${path_in}. Checked extensions: ${genome_extensions[*]}"
        exit 1
    fi

else
    pokaz_error "Error: Directory ${path_in} does not exist."
    exit 1
fi


# ----------------------------------------------
# Optional paths paths

path_ref="${path_ref:-${path_in}}"
path_ref=$(add_symbol_if_missing "$path_ref" "/")

# path_cons="${path_cons:-${path_out}consensus/}"
# path_cons=$(add_symbol_if_missing "$path_cons" "/")

# path_chrom="${path_chrom:-${path_out}chromosomes/}"
# path_chrom=$(add_symbol_if_missing "$path_chrom" "/")

# path_parts="${path_parts:-${path_out}parts/}"
# path_parts=$(add_symbol_if_missing "$path_parts" "/")

path_cons="${path_out}consensus/"
path_chrom="${path_out}chromosomes/"
path_parts="${path_out}parts/"

# ----------------------------------------------
# Make folders

mkdir -p "${path_out}"
mkdir -p "${path_cons}"
mkdir -p "${path_chrom}"
mkdir -p "${path_parts}"

# ----------------------------------------------
# Flags
path_flags="${path_out}flags/"
mkdir -p "${path_flags}"

# ----------------------------------------------
# Number pf chromosomes

if [ -z "${nchr}" ] && [ -z "${nchr_ref}" ]; then  # Both nchr and nchr_ref are not defined.
    # Try to define options for number of chromosomes

    if [ "${mode_pangen}" == "${name_mode_pre}" ]; then  # PRE mode
        option_nchr=" --all.chr T"
        option_nchr_ref=" --all.chr T "
    else   # REF and MSA mode
        # Define from files

        # Initialize an array to hold the counts
        counts=()

        # Loop through all extensions
        for ext in "${genome_extensions[@]}"; do
            # Find all files with the current extension
            for file in "$path_in"/*."$ext"; do
                # Check if the file exists
                if [[ -f "$file" ]]; then
                    # Count the number of lines starting with >
                    count=$(grep -c '^>' "$file")
                    # Append the count to the array
                    counts+=("$count")
                fi
            done
        done

        # If the number of chromosomes in genomes vary, raise an error
        if [[ $(echo "${counts[@]}" | tr ' ' '\n' | uniq | wc -l) -eq 1 ]]; then
            nchr=${counts[0]}
            option_nchr=" --n.chr ${nchr} "
            option_nchr_ref=" --n.chr ${nchr} "
        else
            pokaz_error "Genomes have different number of chromosomes. Please change files or specify the number -nchr."
            exit 1
        fi
    fi

elif [ -z "${nchr}" ] && [ ! -z "${nchr_ref}" ]; then  # nchr_ref is defined.
    pokaz_error "Error: -nchr should be defined."
    exit 1

elif [ ! -z "${nchr}" ] && [ -z "${nchr_ref}" ]; then  # nchr is defined.
    option_nchr=" --n.chr ${nchr} "
    option_nchr_ref=" --n.chr ${nchr} "

else  # Both are defined
    option_nchr=" --n.chr ${nchr} "
    option_nchr_ref=" --n.chr ${nchr_ref} "
fi

# ----------------------------------------------
# Handling parts

# Default parameters
p_ident="${p_ident:-85}"  
part_len="${part_len:-5000}"  

# Filter repeats

if [ -z "${purge_reps}" ]; then
    option_purge_reps=" "
else
    option_purge_reps=" --purge.reps T"
fi


if [ -z "${flag_rev}" ]; then
    option_mirror=" "
else
    path_parts="${path_out}parts_mirror/"
    option_mirror=" --purge.reps T"
fi


# Reverse parts
if [ -z "${rev}" ]; then
    option_rev=" "
else
    option_rev=" --rev "
fi

# ----------------------------------------------
# Files for Filtration

if [ -z "${comb_file}" ]; then
    option_combinations=" "
else
    option_combinations=" --combinations ${comb_file}"
fi

if [ -z "${acc_file}" ]; then
    option_accessions=" "
else
    option_accessions=" --accessions ${acc_file} "
fi

# ----------------------------------------------
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
path_logs="${path_out}logs/"
mkdir -p ${path_logs}

# ----------------------------------------------
# Define the uniquie path for the current run

# Search for directories named "pangenX" and extract the numbers X
max_number=-1
for dir in ${path_logs}pangen*; do
    if [[ -d $dir && $dir =~ pangen([0-9]+) ]]; then
        number=${BASH_REMATCH[1]}
        if (( number > max_number )); then
            max_number=$number
        fi
    fi
done

# Output the largest number
if [[ $max_number -eq -1 ]]; then
    pangen_run_num=0
else
    pangen_run_num=$((max_number + 1))
fi

# if [[ $start_step -eq 0 ]]; then
#     pangen_run_num=$((max_number + 1))
# else
#     pangen_run_num=${max_number}
# fi

# Path for logs of this script

pangen_run_num=0

path_log="${path_logs}pangen${pangen_run_num}/"
echo ${path_log}
mkdir -p ${path_log}


# File with steps logs
file_log="${path_log}steps.log"
> "${file_log}"

#TODO: which command was run



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
#            Check reference genomes for MSA
# ----------------------------------------------------------------------------


if [ "${mode_pangen}" == "${name_mode_msa}" ]; then


    # ----------------------------------------------
    # Check the number of reference genomes for randomisation

    if [ -z "${ref_num}" ]; then
        ref_num=2
    fi

    if ! [[ "$ref_num" =~ ^[0-9]+$ ]]; then
        with_level 0 pokaz_error "Error: Number of reference genomes is not a number."
        exit 1
    fi
    with_level 2 pokaz_message "Number of genomes for randomisation: ${ref_num}"


    # ----------------------------------------------
    # Get names of all genomes form the input folder
    acc_set=($(find "$path_in" -type f -exec bash -c 'basename "$1" .${1##*.}' _ {} \;))

    # # Print the genomes for verification
    # for ref in "${acc_set[@]}"; do
    #     echo "$ref"
    # done

    # Check if ref_num is greater than the number of genomes in acc_set
    if (( ref_num > ${#acc_set[@]} )); then
        with_level 0 pokaz_error "Error: ref_num ($ref_num) is greater than the number of available genomes (${#acc_set[@]}) in acc_set."
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

        # if (( ${#refs_all[@]} > ref_num )); then
        #     with_level 1 pokaz_attention "Parameter -nref is ignored; it doesn't match number of genome in -refs"
        # fi
    fi

    # Print the selected reference genomes for verification
    with_level 2 pokaz_message "Names of genomes for randomisation: $(IFS=,; echo "${refs_all[*]}")"

    # Check that all of the genomes are in the folder
    for ref in "${refs_all[@]}"; do
        if [[ ! " ${acc_set[@]} " =~ " ${ref} " ]]; then
            with_level 0 pokaz_attention "Error: Genome ${ref} is not in ${path_in}"
            exit 1
        fi
    done
else
    refs_all=("${ref_name}")
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

step_file="$path_flags/step${step_num}_done"
if [ $start_step -le ${step_num} ] || [ ! -f ${step_file} ]; then

    with_level 1 pokaz_stage "Step ${step_num}. Genomes into chromosomes."

    # Clean up the output folders
    if   [ "$clean_dir" = "T" ]; then 
        rm -f ${path_chrom}*fasta
    fi

    # Path for logging
    path_log_step="${path_log}step${step_num}_query_01/"

    # Run the step
    Rscript pangen/query_01_to_chr.R --path.in ${path_in} --path.out ${path_chrom} \
            --cores ${cores}  \
            ${option_nchr} \
            ${option_accessions} \
            --path.log ${path_log_step} --log.level ${log_level}

    # Done
    touch "${step_file}"
    with_level 1 pokaz_message "Step is done."

    # If only one step
    if   [ "$one_step" = "T" ]; then 
        exit 0
    fi
    
fi

((step_num = step_num + 1))


# ----------------------------------------------
# If to find ORFs
if [ ! -z "${flag_orf}" ]; then

    with_level 1 pokaz_stage "Additional step. Get all ORFs."

    # Make ORF folder
    path_orf="${path_out}orf/"
    mkdir -p "${path_orf}"

    # Clean up the output folders
    if   [ "$clean_dir" = "T" ]; then 
        rm -f ${path_orf}*fasta
    fi

    # Path for logging
    path_log_step="${path_log}step${step_num}_query_01_orf/"

    # Run the step
    Rscript pangen/query_01_to_orf.R --path.in ${path_in} --path.orf ${path_orf} \
            --cores ${cores}  \
            ${option_nchr} \
            ${option_accessions} \
            --path.log ${path_log_step} --log.level ${log_level}

    # Done
    touch "${step_file}"
    with_level 1 pokaz_message "Step is done."

    # If only one step
    if   [ "$one_step" = "T" ]; then 
        exit 0
    fi
        

fi



# ----------------------------------------------
# Split query chromosomes into parts
step_file="$path_flags/step${step_num}_done"
if [ $start_step -le ${step_num} ] || [ ! -f ${step_file} ]; then

    with_level 1 pokaz_stage "Step ${step_num}. Chromosomes into parts."

    # Clean up the output folders
    if   [ "$clean_dir" = "T" ]; then 
        rm -f ${path_parts}*fasta
    fi    

    # Path for logging
    path_log_step="${path_log}step${step_num}_query_02/"

    # Run the step
    Rscript pangen/query_02_to_parts.R --path.chr ${path_chrom} \
            --path.parts ${path_parts} --part.len $part_len --cores ${cores} \
            ${option_purge_reps} ${option_rev} \
            ${option_nchr} \
            ${option_mirror} \
            --path.log ${path_log_step} --log.level ${log_level}

    # Done
    touch "${step_file}"
    with_level 1 pokaz_message "Step is done."

    # If only one step
    if   [ "$one_step" = "T" ]; then 
        exit 0
    fi
    
fi

((step_num = step_num + 1))


# ┌───────────────────────────────────────────────────────────────────────────┐
# │                         REFERENCE-BASED                                   │
# └───────────────────────────────────────────────────────────────────────────┘

for ref0 in "${refs_all[@]}"; do


    with_level 1  pokaz_attention "Reference ${ref0}"

    # New paths
    path_blast_parts=${path_out}blast_parts_${ref0}/
    path_alignment=${path_out}alignments_${ref0}/
    path_gaps=${path_out}blast_gaps_${ref0}/

    # ----------------------------------------------
    # Split reference fasta into chromosomes if additionally needed
    if [[ "${path_in}" != "$path_ref" ]]; then

        step_file="$path_flags/step${step_num}_${ref0}_done"
        if [ $start_step -le ${step_num} ] || [ ! -f ${step_file} ]; then
            with_level 1 pokaz_stage "Step ${step_num}. Reference genome into chromosomes."

            # Temporary file to analyse only the reference genome from the folder
            file_acc_ref=${path_cons}ref_acc.txt
            echo "${ref0}" > ${file_acc_ref}

            # Path for logging
            path_log_step="${path_log}step${step_num}_query_01_ref/"

            echo ${nchr_ref}
            
            Rscript pangen/query_01_to_chr.R --all.chr T \
                    --path.in ${path_ref} --path.out ${path_chrom}   \
                    --cores ${cores} \
                    --accessions ${file_acc_ref} \
                    ${nchr_ref} \
                    --path.log ${path_log_step} --log.level ${log_level}

            # Remove the temporary file
            rm ${file_acc_ref}

            # Done
            touch "${step_file}"
            with_level 1 pokaz_message "Step is done."

            # If only one step
            if   [ "$one_step" = "T" ]; then 
                exit 0
            fi
        fi

        ((step_num = step_num + 1))
    fi

    # ----------------------------------------------
    # Blast parts on the reference genome

    step_file="$path_flags/step${step_num}_${ref0}_done"
    if [ $start_step -le ${step_num} ] || [ ! -f "$step_file" ]; then

        with_level 1 pokaz_stage "Step ${step_num}. BLAST of parts against the reference genome."
        # with_level 1 pokaz_message "NOTE: if this stage takes relatively long, use -purge_repeats -s 2 to mask highly repetative regions"

        if [[ $option_purge_reps =~ ^[[:space:]]*$ ]]; then
            with_level 1 pokaz_message "NOTE: if this stage takes relatively long, use -purge_repeats -s $((${step_num}-1)) to mask highly repetative regions"
        fi
    
        # ---- Clean up the output folders ----
        if   [ "$clean_dir" = "T" ]; then 
            rm -f ${path_blast_parts}*txt
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

        # Iterate over each file
        for file in "${files[@]}"; do
          # Check if the BLAST database files already exist
          if [ ! -f "${file}.nin" ]; then
              makeblastdb -in ${file} -dbtype nucl > /dev/null
          fi
        done

        # ---- Run BLAST ----

        with_level 2 pokaz_message "Run BLAST."    

        # Logs for the BLAST
        path_log_step="${path_log}step${step_num}_query_03/"
        make_dir ${path_log_step}

        # Blast parts on the reference genome
        ./pangen/query_03_blast_parts.sh \
                -path_chrom ${path_chrom} \
                -path_parts ${path_parts} \
                -path_result ${path_blast_parts} \
                -ref ${ref0} \
                ${option_one2one} \
                -p_ident ${p_ident} \
                -cores ${cores} \
                -log_path ${path_log_step}

        # Done
        touch "${step_file}"
        with_level 1 pokaz_message "Step is done."

        # If only one step
        if   [ "$one_step" = "T" ]; then 
            exit 0
        fi
    fi

    ((step_num = step_num + 1))

    # ----------------------------------------------
    # Search for the mirror universe

    if [ ! -z "${flag_rev}" ]; then
        echo "Pannagram's Mirror mode is done."
        exit
    fi


    # ----------------------------------------------
    # First round of alignments

    step_file="$path_flags/step${step_num}_${ref0}_done"
    if [ $start_step -le ${step_num} ] || [ ! -f "$step_file" ]; then

        with_level 1 pokaz_stage "Step ${step_num}. Alignment-1: Remaining syntenic (major) matches."

        # Clean up the output folders
        if   [ "$clean_dir" = "T" ]; then 
            rm -f "${path_alignment}*rds"
            rm -f "${path_gaps}*fasta"
        fi    
        
        # Logs for the BLAST
        path_log_step="${path_log}step${step_num}_synteny_01/"
        make_dir ${path_log_step}
        
        # Run the step
        Rscript pangen/synteny_01_majoir.R --path.blast ${path_blast_parts} --path.aln ${path_alignment} \
                --ref ${ref0}   \
                --path.gaps  ${path_gaps} --path.chr ${path_chrom} \
                --cores ${cores} \
                ${option_one2one_r} \
                --path.log ${path_log_step} --log.level ${log_level}

        # Done
        touch "${step_file}"
        with_level 1 pokaz_message "Step is done."

        # If only one step
        if   [ "$one_step" = "T" ]; then 
            exit 0
        fi

        # Clean up the output folders of previous stages
        # If the first round of alignment didn't have any errors - remove the blast which was needed for it
        # rm ${path_blast_parts}
    fi

    ((step_num = step_num + 1))


    # ----------------------------------------------
    # Step 6: Plotting

    step_file="$path_flags/step${step_num}_${ref0}_done"
    if [ $start_step -le ${step_num} ] || [ ! -f "$step_file" ]; then

        with_level 1 pokaz_stage "Step ${step_num}. Plotting the results."

        # Logs for the Plot
        path_log_step="${path_log}step${step_num}_synteny_01_plot_${ref0}/"
        make_dir ${path_log_step}

        Rscript pangen/synteny_01_plot.R \
                --ref ${ref0} \
                --path_ref ${path_ref} \
                --path_chr ${path_chrom} \
                --path_out ${path_out} \
                --algn_path ${path_alignment} \
                --path.log ${path_log_step} \
                --log.level ${log_level} \
                --cores ${cores}

        # Done
        touch "${step_file}"
        with_level 1 pokaz_message "Step is done."

        # If only one step
        if   [ "$one_step" = "T" ]; then 
            exit 0
        fi
    fi

    ((step_num = step_num + 1))


    # ========== CHECK MODE ==========
    if [ "${mode_pangen}" == "${name_mode_pre}" ]; then 
        with_level 0 pokaz_message "Pannagram's PREliminary mode is done."
        exit
    fi 
        
    # ----------------------------------------------
    # Blast regions between synteny blocks

    step_file="$path_flags/step${step_num}_${ref0}_done"
    if [ $start_step -le ${step_num} ] || [ ! -f "$step_file" ]; then

        with_level 1 pokaz_stage "Step ${step_num}. BLAST of gaps between syntenic matches."

        # Logs for the BLAST
        path_log_step="${path_log}step${step_num}_synteny_02_${ref0}/"
        make_dir ${path_log_step}

        # Run BLAST for gaps
        ./pangen/synteny_02_blast_gaps.sh \
                -path_gaps ${path_gaps} \
                -cores ${cores} \
                -log_path ${path_log_step}

        # Done
        touch "${step_file}"
        with_level 1 pokaz_message "Step is done."

        # If only one step
        if   [ "$one_step" = "T" ]; then 
            exit 0
        fi
    fi

    ((step_num = step_num + 1))

    # ----------------------------------------------
    # Second round of alignments

    step_file="$path_flags/step${step_num}_${ref0}_done"
    if [ $start_step -le ${step_num} ] || [ ! -f "$step_file" ]; then

        with_level 1 \
            pokaz_stage "Step ${step_num}. Alignment-2: Fill the gaps between synteny blocks."

        # Logs for the merging gaps
        path_log_step="${path_log}step${step_num}_synteny_03_${ref0}/"
        make_dir ${path_log_step}

        Rscript pangen/synteny_03_merge_gaps.R --ref ${ref0} \
                --path.aln ${path_alignment}  --path.chr ${path_chrom}\
                --path.gaps ${path_gaps}   \
                --cores ${cores} \
                --path.log ${path_log_step} \
                --log.level ${log_level}

        # # If the second round of alignment didn't have any errors - remove the blast which was needed for it
        # rm -rf ${path_gaps}
        # ls ${path_alignment}*maj*
        # rm -rf ${path_alignment}*maj*

        # Done
        touch "${step_file}"
        with_level 1 pokaz_message "Step is done."

        # If only one step
        if   [ "$one_step" = "T" ]; then 
            exit 0
        fi
    fi

    ((step_num = step_num + 1))

    # ----------------------------------------------
    # Create a consensus

    step_file="$path_flags/step${step_num}_${ref0}_done"
    if [ $start_step -le ${step_num} ] || [ ! -f "$step_file" ]; then

        with_level 1 \
            pokaz_stage "Step ${step_num}. Combine reference-based alignments by chromosomes."

        # Logs for creating the consensus
        path_log_step="${path_log}step${step_num}_comb_01_${ref0}/"
        make_dir ${path_log_step}

        Rscript pangen/comb_01_one_ref.R \
                --path.cons ${path_cons} \
                --path.aln ${path_alignment} \
                --pref ${ref0}  \
                --cores ${cores} \
                --path.chr ${path_chrom} \
                --path.log ${path_log_step} \
                --log.level ${log_level}

        # Done
        touch "${step_file}"
        with_level 1 pokaz_message "Step is done."

        # If only one step
        if   [ "$one_step" = "T" ]; then 
            exit 0
        fi
    fi

    ((step_num = step_num + 1))


done
# ----------------------------------------------------------------------------
#            CHECK MODE
# ----------------------------------------------------------------------------
if [ "${mode_pangen}" != "${name_mode_msa}" ]; then 
    echo "Pannagram's REFRENCE-based mode is done."
    exit
fi



#./pangen_new.sh -path_in /Volumes/Samsung_T5/vienn/pannagram_test/symA -path_out /Volumes/Samsung_T5/vienn/pannagram_test/symA_test2 -log 2


# ╔═══════════════════════════════════════════════════════════════════════════════════╗       ______
# ║                                                                                   ║      /|_||_\`.__
# ║                                  MSA                                              ║     (   _    _ _\
# ║                                                                                   ║     =`-(_)--(_)-'
# ╚═══════════════════════════════════════════════════════════════════════════════════╝ 

with_level 1 pokaz_stage "* MGA is strting * *"

# ----------------------------------------------
# Run consensus for a pair of files
ref0=${refs_all[0]}

step_file="$path_flags/step${step_num}_done"
if [ $start_step -le ${step_num} ] || [ ! -f "$step_file" ]; then

    for ((i = 1; i < ${#refs_all[@]}; i++)); do
        ref1=${refs_all[i]}

        with_level 1 pokaz_stage "Step ${step_num}. Randomisation of alignments. Combine two references: ${ref0} and ${ref1}."
        # ref1=${ref1//_/$'-'}

        # Logging
        path_log_step="${path_log}step${step_num}_comb_02_${ref0}_${ref1}/"
        make_dir ${path_log_step}
        
        Rscript pangen/comb_02_two_refs.R \
                --path.cons ${path_cons} \
                --ref0 ${ref0} \
                --ref1 ${ref1} \
                --cores ${cores} \
                --path.log ${path_log_step} \
                --log.level ${log_level}

    done

    touch "${step_file}"
    with_level 1 pokaz_message "Step is done."
fi

((step_num = step_num + 1))


# ----------------------------------------------
# Common gaps
step_file="$path_flags/step${step_num}_done"
if [ $start_step -le ${step_num} ] || [ ! -f "$step_file" ]; then

    with_level 1 pokaz_stage "Step ${step_num}. Find Positions of Common Gaps in the Reference-Free Multiple Genome Alignment."

    # Logging
    path_log_step="${path_log}step${step_num}_comb_03/"
    make_dir ${path_log_step}

    Rscript pangen/comb_03_find_gaps.R \
            --path.cons ${path_cons} \
            --ref.pref ${ref0} \
            --cores ${cores} \
            --path.log ${path_log_step} \
            --log.level ${log_level}

    # Done
    touch "${step_file}"
    with_level 1 pokaz_message "Step is done."

    # If only one step
    if   [ "$one_step" = "T" ]; then 
        exit 0
    fi
fi

((step_num = step_num + 1))


# ----------------------------------------------
# Create sequences to run MAFFT and perform some small alignments
path_mafft_in="${path_out}mafft_in/"
if [ ! -d "$path_mafft_in" ]; then
    mkdir -p "$path_mafft_in"
fi

step_file="$path_flags/step${step_num}_done"
if [ $start_step -le ${step_num} ] || [ ! -f "$step_file" ]; then

    with_level 1 pokaz_stage "Step ${step_num}. Prepare sequences for MAFFT."

    # Logging
    path_log_step="${path_log}step${step_num}_comb_04/"
    make_dir ${path_log_step}

    Rscript pangen/comb_04_prepare_aln.R \
            --path.cons ${path_cons} \
            --ref.pref ${ref0} \
            --cores ${cores} \
            --path.chromosomes ${path_chrom} \
            --path.mafft.in ${path_mafft_in} \
            --path.log ${path_log_step} \
            --log.level ${log_level}

    # Done
    touch "${step_file}"
    with_level 1 pokaz_message "Step is done."

    # If only one step
    if   [ "$one_step" = "T" ]; then 
        exit 0
    fi
fi

((step_num = step_num + 1))

# ----------------------------------------------
# Run MAFFT
path_mafft_out="${path_out}mafft_out/"
if [ ! -d "$path_mafft_out" ]; then
    mkdir -p "$path_mafft_out"
fi

step_file="$path_flags/step${step_num}_done"
if [ $start_step -le ${step_num} ] || [ ! -f "$step_file" ]; then

    with_level 1 pokaz_stage "Step ${step_num}. Run MAFFT."

    # Logging
    path_log_step="${path_log}step${step_num}_comb_05/"
    make_dir ${path_log_step}

    ./pangen/comb_05_run_mafft.sh \
            -cores ${cores} \
            -path_mafft_in ${path_mafft_in} \
            -path_mafft_out ${path_mafft_out} \
            -log_path ${path_log_step}

    # Done
    touch "${step_file}"
    with_level 1 pokaz_message "Step is done."

    # If only one step
    if   [ "$one_step" = "T" ]; then 
        exit 0
    fi
fi

((step_num = step_num + 1))
echo ${step_num}
# ----------------------------------------------
# Combine all together

step_file="$path_flags/step${step_num}_done"
if [ $start_step -le ${step_num} ] || [ ! -f "$step_file" ]; then

    with_level 1 pokaz_stage "Step ${step_num}. Combine all alignments together into the final one."

    # Logging
    path_log_step="${path_log}step${step_num}_comb_06/"
    make_dir ${path_log_step}

    Rscript pangen/comb_06_final_aln.R  \
            --cores ${cores} \
            --ref.pref ${ref0} \
            --path.mafft.in ${path_mafft_in} \
            --path.mafft.out ${path_mafft_out} \
            --path.cons ${path_cons} \
            --path.log ${path_log_step} \
            --log.level ${log_level}

    # Done
    touch "${step_file}"
    with_level 1 pokaz_message "Step is done."

    # If only one step
    if   [ "$one_step" = "T" ]; then 
        exit 0
    fi
fi

((step_num = step_num + 1))


# ----------------------------------------------
# Get synteny blocks

aln_type='msa_'
step_file="$path_flags/step${step_num}_done"
if [ $start_step -le ${step_num} ] || [ ! -f "$step_file" ]; then

    with_level 1 pokaz_stage "Step ${step_num}. Get synteny blocks."

    Rscript analys/analys_01_blocks.R \
            --path.cons ${path_cons} \
            --ref.pref  ${ref0} \
            --cores ${cores} \
            --aln.type ${aln_type}

    # Done
    touch "${step_file}"
    with_level 1 pokaz_message "Step is done."

    # If only one step
    if   [ "$one_step" = "T" ]; then 
        exit 0
    fi
fi

((step_num = step_num + 1))

# if [ $start_step -eq 0 ]; then
#     rm -f "$FLAG_DIR"/.*
#     echo "Script completed successfully"
# fi









