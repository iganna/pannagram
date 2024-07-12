# This escript performs all stages of the alignment of genomes, 
# when the reference genome is already identified

# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

source utils/chunk_error_control.sh

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source utils/utils_bash.sh

# Function to display help message
print_usage() {
    pokaz_help

    cat << EOF
Usage: ${0##*/} [-h] [-s STAGE] [-cores CORES] 
                [-ref REF_NAME] [-path_in INPUT_FOLDER] [-path_out OUTPUT_FOLDER]
                [-path_ref PATH_REF] [-path_chrom PATH_CHROM] [-path_parts PATH_PARTS] [-path_cons PATH_CONSENSUS] 
                [-sort_len] [-one2one] [-accessions ACC_ANAL] 
                [-part_len PART_LEN] [-p_ident P_IDENT] [-purge_repeats]
                
This script performs alignment of query genomes to the reference genome.

Options:
    -h, -help                  Display this help message and exit.
    -s, -stage STAGE            Specify the stage from which to run. If not provided, the last interrupted stage will be re-run.
    -cores CORES                Number of cores for parallel processing. Default is 1.

    # Required parameters
    -ref REF_NAME               Name of the reference genome: REF_NAME.fasta.
    -path_in INPUT_FOLDER       Folder to the directory with query genomes. Possible file types: fasta, fna, fa, fas
    -path_out OUTPUT_FOLDER     Folder where all results and intermediate files will appear.
    
    # Optional paths
    -path_ref PATH_REF          Path where the reference genome is stored. Do not provide if it's the same folder as the path with query genomes.
    -path_chrom PATH_CHROM      Path to the folder with individual chromosomes in separate files. 
    -path_parts PATH_PARTS      Path to the folder with files of chromosomal parts.
    -path_cons PATH_CONSENSUS   Path to the consensus folder.
    -nchr N_CHR                 Number of chromosomes in every.
    -nchr_ref N_CHR_REF         Number of chromosomes in the reference genome.
    -nchr_acc N_CHR_QUERY     Number of chromosomes in the query genome.

    # Input Design Handling
    -sort_len                   Flag to sort chromosomes by length.
    -one2one                    Flag for straightforward pairwise alignment of chromosomes, matching each chromosome sequentially with its corresponding one (NOT all vs all, as by default).
    -accessions ACC_ANAL        File with accessions to analyze. Accessions should be in rows.

    # Tuning parameters    
    -p_ident P_IDENT            Percentage identity threshold (default: 85).
    -part_len PART_LEN          Fragments to which each chromosome should be cut (default value: 5000).
    -purge_repeats              Enable filtration of repeats (default: disabled).

    
Examples:
    ${0##*/} -ref '0' -path_in 'input_genomes' -path_out 'output_folder'

EOF
}


# ----------------------------------------------------------------------------
#            PARAMETERS
# ----------------------------------------------------------------------------

# Default values of parameters
start_step=0

unrecognized_options=()

while [ $# -gt 0 ]
do
    case $1 in
        -h | -help ) print_usage; exit ;;
        -s | -stage ) start_step="$2"; shift 2 ;;  # stage from which to run, when the stage is not provided - the last interrupted stage withh be re-run
        -log)         log_level=$2;    shift 2 ;;  # path to the output
        -cores)       cores=$2;        shift 2 ;;

        -path_out) pref_global=$2; shift 2 ;;  # path to the output
        -path_in)  path_in=$2;     shift 2 ;;  # path with all genomes in fasta format

        -ref)      ref_name=$2; shift 2 ;;  # name of the reference genome
	    -path_ref) path_ref=$2; shift 2 ;;  # dont provide if it's the same folder as the path with query genomes

        -path_chrom) path_chrom=$2;     shift 2 ;;  # path to the folder with individual chromosomes in separate files
        -path_parts) path_parts=$2;     shift 2 ;;  # path to the folder with chromosomal parts
        -path_cons)  path_consensus=$2; shift 2 ;;  # path to the consensus folder

        -nchr)      nchr=$2;     shift 2 ;;  # number of chromosomes
        -nchr_ref)  nchr_ref=$2; shift 2 ;;  # number of chromosomes in the reference genome
        -nchr_acc)  nchr_acc=$2; shift 2 ;;  # number of chromosome in the query genome

        -part_len)  part_len=$2; shift 2 ;;  # fragments to which each chromosome should be cut, has a default value 5000
        -p_ident)   p_ident=$2;  shift 2 ;;  # percent of identity

        -sort_by_len)   sort_chr_len="T"; shift 1 ;;  # flag whether to sort chromosomes by length or not
        -one2one)       one2one="T";      shift 1 ;;  # conpare all to all or not
        -purge_reps )   purge_reps="T"    shift 1 ;;  # filtration of repeats, default - not
        -rev )          flag_rev="T"      shift 1 ;;  # reverce parts

        -accessions) acc_anal=$2;    shift 2 ;;  # file with accessions to analyse
        -combinations) comb_anal=$2; shift 2 ;;  # file with chromosomal combinations to analyse: first column - query, second column - reference(base)
    
        *) 
            unrecognized_options+=("$1"); shift 1 ;;
    esac
done


# Output of Unrecognized Parameters
if [[ ${#unrecognized_options[@]} -gt 0 ]]; then
    print_usage
    echo "Unrecognized options:"
    for option in "${unrecognized_options[@]}"; do
        echo "$option"
    done
fi

# ---- Check of missimg parameters

check_missing_variable "pref_global"
check_missing_variable "path_in"
check_missing_variable "ref_name"

# ---- if some parameters zre not given - set the default value

# Basic parameters

cores="${cores:-1}"  # Number of cores
p_ident="${p_ident:-85}"  
part_len="${part_len:-5000}"  

sort_chr_len="${sort_chr_len:-F}"
purge_reps="${purge_reps:-F}"
flag_rev="${flag_rev:-0}"

one2one="${one2one:-F}"
if [ "$one2one" = "T" ]; then
    all2all="F"
else
    all2all="T"
fi

acc_anal="${acc_anal:-NULL}"   # Set of accessions to analyse

# ---- Paths ----
# Required

path_in=$(add_symbol_if_missing "$path_in" "/")
pref_global=$(add_symbol_if_missing "$pref_global" "/")

# Could be defined
path_chrom="${path_chrom:-${pref_global}chromosomes/}"
path_chrom=$(add_symbol_if_missing "$path_chrom" "/")

path_parts="${path_parts:-${pref_global}parts/}"
path_parts=$(add_symbol_if_missing "$path_parts" "/")

path_ref="${path_ref:-${path_in}}"
path_ref=$(add_symbol_if_missing "$path_ref" "/")

path_consensus="${path_consensus:-${pref_global}consensus/}"
path_consensus=$(add_symbol_if_missing "$path_consensus" "/")
# echo "consensus${path_consensus}"
if [ ! -d "$path_consensus" ]; then
    mkdir -p "$path_consensus"
fi

# New paths
path_blast_parts=${pref_global}blast_parts_${ref_name}/
path_alignment=${pref_global}alignments_${ref_name}/
path_gaps=${pref_global}blast_gaps_${ref_name}/


# ----------------------------------------------------------------------------
#           STEPS
# ----------------------------------------------------------------------------

source utils/chunk_steps.sh

echo "Pangen-pre starts from step ${start_step}"

# ----------------------------------------------------------------------------
#           LOGS
# ----------------------------------------------------------------------------

source utils/chunk_logging_ref.sh

> "${file_log_ref}"  # Clean up the logging file

# ----------------------------------------------------------------------------
#           NUMBER OF CHROMOSOMES
# ----------------------------------------------------------------------------

# Check conditions and process values
if [[ -n $nchr ]]; then
    if [[ -n $nchr_ref || -n $nchr_acc ]]; then
        pokaz_attention "WARNING: If -nchr is set, -nchr_ref and -nchr_acc should not be set. -nchr has the priority."
    fi
    nchr_ref=$nchr
    nchr_acc=$nchr
fi


# Number of chromosomes handling
if [ -n "${nchr_acc}" ]; then
    nchr_acc_option=" --n.chr ${nchr_acc} "
else
    nchr_acc_option=" --all.chr T "
fi


if [ -n "${nchr_ref}" ]; then
    nchr_ref_option=" --n.chr ${nchr_ref} "
else
    nchr_ref_option=" --all.chr T "
fi


# ----------------------------------------------------------------------------
#           MAIN PIPELINE
# ----------------------------------------------------------------------------

step_num=1
# ----------------------------------------------
# Split query fasta into chromosomes

if [ $start_step -le ${step_num} ] || [ ! -f "$path_flags/step${step_num}_done" ]; then

    log_message 1 "$log_level" "$file_log_ref" pokaz_stage "Step ${step_num}. Genomes into chromosomes."

    # Clean up the output folders
    rm -rf ${path_chrom}

    # Path for logging
    path_log_step="${path_log_ref}step${step_num}_query_01/"

    # Run the step
    Rscript pangen/query_01_to_chr.R --path.in ${path_in} --path.out ${path_chrom} \
            --sort ${sort_chr_len} --cores ${cores} --acc.anal ${acc_anal} \
            ${nchr_acc_option} \
            --path.log ${path_log_step} --log.level ${log_level}

    # Done
    touch "$path_flags/step${step_num}_done"
    log_message 1 "$log_level" "$file_log_ref" pokaz_message "Step is done."
    
fi

((step_num = step_num + 1))

# ----------------------------------------------
# Split query chromosomes into parts
if [ $start_step -le ${step_num} ] || [ ! -f "$path_flags/step${step_num}_done" ]; then

    log_message 1 "$log_level" "$file_log_ref" pokaz_stage "Step ${step_num}. Chromosomes into parts."

    # Clean up the output folders
    rm -rf ${path_parts}

    # Path for logging
    path_log_step="${path_log_ref}step${step_num}_query_02/"

    # Run the step
    Rscript pangen/query_02_to_parts.R --path.chr ${path_chrom} \
            --path.parts ${path_parts} --part.len $part_len --cores ${cores} \
            --purge_reps ${purge_reps} --rev ${flag_rev} \
            ${nchr_acc_option} \
            --path.log ${path_log_step} --log.level ${log_level}

    # Done
    touch "$path_flags/step${step_num}_done"
    log_message 1 "$log_level" "$file_log_ref" pokaz_message "Step is done."
    
fi

((step_num = step_num + 1))

# ----------------------------------------------
# Split reference fasta into chromosomes if additionally needed
if [[ "${path_in}" != "$path_ref" ]]; then

    if [ $start_step -le ${step_num} ] || [ ! -f "$path_flags/step${step_num}_done_${ref_name}" ]; then
        log_message 1 "$log_level" "$file_log_ref" pokaz_stage "Step ${step_num}. Reference genome into chromosomes."

        # Temporary file to analyse only the reference genome from the folder
        file_acc_ref=${path_consensus}ref_acc.txt
        echo "${ref_name}" > ${file_acc_ref}

        # Path for logging
        path_log_step="${path_log_ref}step${step_num}_query_01_ref/"
        
        Rscript pangen/query_01_to_chr.R --all.chr T \
                --sort ${sort_chr_len}
                --path.in ${path_ref} --path.out ${path_chrom}   \
                --cores ${cores} --acc.anal ${file_acc_ref} \
                ${nchr_ref} \
                --path.log ${path_log_step} --log.level ${log_level}

        # Remove the temporary file
        rm ${file_acc_ref}

        # Done
        touch "$path_flags/step${step_num}_done_${ref_name}"
        log_message 1 "$log_level" "$file_log_ref" pokaz_message "Step is done."

    fi

    ((step_num = step_num + 1))
fi

# ----------------------------------------------
# Blast parts on the reference genome

if [ $start_step -le ${step_num} ] || [ ! -f "$path_flags/step${step_num}_done_${ref_name}" ]; then

    log_message 1 "$log_level" "$file_log_ref" pokaz_stage "Step ${step_num}. BLAST of parts against the reference genome."
    log_message 1 "$log_level" "$file_log_ref" pokaz_message "NOTE: if this stage takes relatively long, use -purge_repeats -s 2 to mask highly repetative regions"

    # ---- Clean up the output folders ----
    rm -rf ${path_blast_parts}

    # ---- Create a database on the reference genome ----
    log_message 2 "$log_level" "$file_log_ref" pokaz_message "Create BLAST databases."

    for file in ${path_chrom}${ref_name}_chr*fasta ; do
      # Check if the BLAST database files already exist
      if [ ! -f "${file}.nin" ]; then
          makeblastdb -in ${file} -dbtype nucl > /dev/null
      fi
    done

    # ---- Run BLAST ----

    log_message 2 "$log_level" "$file_log_ref" pokaz_message "Run BLAST."    

    # Logs for the BLAST
    path_log_step="${path_log_ref}step${step_num}_query_03/"
    make_dir ${path_log_step}

    # Blast parts on the reference genome
    ./pangen/query_03_blast_parts.sh -path_ref ${path_chrom} -path_parts ${path_parts} \
            -path_result ${path_blast_parts} -ref_name ${ref_name} \
            -all_vs_all ${all_cmp} -p_ident ${p_ident} -cores ${cores} -log_path ${path_log_step}

    # Done
    touch "$path_flags/step${step_num}_done_${ref_name}"
    log_message 1 "$log_level" "$file_log_ref" pokaz_message "Step is done."
fi

((step_num = step_num + 1))

# ----------------------------------------------
# First round of alignments

if [ $start_step -le ${step_num} ] || [ ! -f "$path_flags/step${step_num}_done_${ref_name}" ]; then

    log_message 1 "$log_level" "$file_log_ref" pokaz_stage "Step ${step_num}. Alignment-1: Remaining syntenic (major) matches."

    # Clean up the output folders
    rm -rf ${path_alignment}
    rm -rf ${path_gaps}

    # Logs for the BLAST
    path_log_step="${path_log_ref}step${step_num}_synteny_01/"
    make_dir ${path_log_step}
    
    # Run the step
    Rscript pangen/synteny_01_majoir.R --path.blast ${path_blast_parts} --path.aln ${path_alignment} \
            --ref ${ref_name}   \
            --path.gaps  ${path_gaps} --path.chr ${path_chrom} \
            --cores ${cores} \
            --path.log ${path_log_step} --log.level ${log_level}

    # Done
    touch "$path_flags/step${step_num}_done_${ref_name}"
    log_message 1 "$log_level" "$file_log_ref" pokaz_message "Step is done."

    # Clean up the output folders of previous stages
    # If the first round of alignment didn't have any errors - remove the blast which was needed for it
    # rm -rf ${path_blast_parts}
fi

((step_num = step_num + 1))


# ----------------------------------------------
# Step 6: Plotting
if [ $start_step -le ${step_num} ] || [ ! -f "$path_flags/step${step_num}_done_${ref_name}" ]; then

    pokaz_stage "Step ${step_num}. Plotting the results."

    Rscript pangen/plot_genome_synteny.R \
            --ref ${ref_name} \
            --path_ref ${path_ref} \
            --path_in ${path_in} \
            --path_out ${pref_global} \
            --algn_path ${path_alignment}

    touch "$path_flags/step${step_num}_done_${ref_name}"
fi

((step_num = step_num + 1))


