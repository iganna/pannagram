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


#./pangen_ref.sh -path_out '../pan_test/tom/' -ref_name '0'  -n_chr_ref 5 -path_in '../pb_updated/' -n_chr_query 5 -all_cmp F -acc_anal 'acc_tom.txt'
#./pangen_ref.sh -path_out '../pan_test/tom/' -ref_name '6046-v1.1'  -n_chr_ref 5 -path_in '../pb_updated/' -n_chr_query 5 -all_cmp F -acc_anal 'acc_tom.txt'

#./pangen_ref.sh -path_out ../pan_test/ly_th/ -ref_name 0 -path_chr_ref ../pan_test/p27/chromosomes/ -n_chr_ref 5 -path_in ../lyrata/ -n_chr_query 8 -all_cmp T -cores 30


# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source utils/utils_bash.sh

# Function to display help message
print_usage() {
    pokaz_help

    cat << EOF
Usage: ${0##*/} -ref REF_NAME -path_in INPUT_FOLDER -path_out OUTPUT_FOLDER
                -nchr_ref N_CHR_REF -nchr_query N_CHR_QUERY
                [-path_ref PATH_CHR_REF] [-path_chrom PATH_CHROM] 
                [-path_parts PATH_PARTS] [-path_cons PATH_CONSENSUS] 
                [-sort_len] [-one2one] [-accessions ACC_FILE] 
                [-part_len PART_LEN] [-p_ident P_IDENT] [-purge_repeats]
                [-h] [-s STAGE] [-cores CORES] [-echo] [-log LOG_LEVEL]
                
This script performs alignment of query genomes to the reference genome.

Options:
    -h, -help                   Display this help message and exit.
    -s, -stage STAGE            Specify the stage from which to run: a number from 1 to 13.
                                If not provided, the last interrupted stage will be re-run.
    -cores CORES                Number of cores for parallel processing. Default is 1.
    -log LOG_LEVEL              The level of logging to be shown on the screen. 
                                Options: 
                                    0 - off
                                    1 - steps (default),
                                    2 - progress
                                    3 - all

        
    * Paths and names (required):
        -path_in INPUT_FOLDER       Folder to the directory with query genomes. 
                                    Possible file types: fasta, fna, fa, fas.
        -path_out OUTPUT_FOLDER     Folder where all results and intermediate files will appear.
        -ref REF_NAME               Name of the reference genome: REF_NAME.fasta. 
                                    Names should NOT contain "_chr" as a substring.
        
    * Numbers of chromosomes (required):
        -nchr_ref N_CHR_REF         Number of chromosomes in the reference genome.
        -nchr_query N_CHR_QUERY     Number of chromosomes in the query genome.

    * Paths (optional):
        -path_ref PATH_CHR_REF      Path where the reference genome is stored. 
                                    Do not provide if it's the same folder as the path with query genomes.
        -path_chrom PATH_CHROM      Path to the folder with individual chromosomes in separate files. 
        -path_parts PATH_PARTS      Path to the folder with files of chromosomal parts.
        -path_cons PATH_CONSENSUS   Path to the consensus folder.

    * Input Design Handling (optional):
        -sort_len                   Flag to sort chromosomes by length.
        -one2one                    Flag for straightforward pairwise alignment of chromosomes,
                                    matching each chromosome sequentially with its corresponding one 
                                    (NOT all vs all, as by default).
        -accessions ACC_FILE        File with accessions to analyze. Accessions should be in rows.
        -combinations COMB_FILE     File with combinations to analyze.

    * Tuning parameters (optional): 
        -p_ident P_IDENT            Percentage identity threshold (default: 85).
        -part_len PART_LEN          Fragments to which each chromosome should be cut (default value: 5000).
        -purge_repeats              Enable filtration of repeats (default: disabled).

    
Examples:
    ${0##*/} -path_out 'output_folder' -path_in 'input_genomes' -ref '0' -nchr_query 5 -nchr_ref 5

EOF
}


# ----------------------------------------------------------------------------
#            PARAMETERS
# ----------------------------------------------------------------------------

unrecognized_options=()

while [ $# -gt 0 ]
do
    case $1 in
        -h | -help ) print_usage; exit ;;
        -s | -stage ) start_step="$2"; shift ;;  # stage from which to run, when the stage is not provided - the last interrupted stage withh be re-run
        -log) log_level=$2; shift ;;  # path to the output
        -cores) cores=$2; shift ;;

        -path_out) path_out=$2; shift ;;  # path to the output
        -path_in) path_in=$2; shift ;;  # path with all genomes in fasta format

        -ref) ref_name=$2; shift ;;  # name of the reference genome
	    -path_ref) path_chr_ref=$2; shift ;;  # dont provide if it's the same folder as the path with query genomes

        -nchr_ref) n_chr_ref=$2; shift ;;  # number of chromosomes in the reference genome
        -nchr_query) n_chr_query=$2; shift ;;  # number of chromosome in the query genome

        -path_chrom) path_chr_acc=$2; shift ;;  # path to the folder with individual chromosomes in separate files
        -path_parts) path_parts=$2; shift ;;  # path to the folder with chromosomal parts
        -path_cons) path_consensus=$2; shift ;;  # path to the consensus folder

        -part_len) part_len=$2; shift ;;  # fragments to which each chromosome should be cut, has a default value 5000

        -sort_by_len) sort_chr_len="T" ;;  # flag whether to sort chromosomes by length or not
        -one2one) all_cmp="F" ;;  # compare all to all or not

        -p_ident) p_ident=$2; shift ;;  # percent of identity
        -purge_repeats ) filter_rep=1 ;;  # filtration of repeats, default - not

        -rev ) flag_rev=1 ;;  # filtration of repeats, default - not

        -accessions) acc_anal=$2; shift ;;  # file with accessions to analyse
        -combinations) comb_anal=$2; shift ;;  # file with chromosomal combinations to analyse: first column - query, second column - reference(base)

        
    
        *) print_usage
            unrecognized_options+=("$1"); shift ;;
    esac
    shift
done

# Output of Unrecognized Parameters
if [[ ${#unrecognized_options[@]} -gt 0 ]]
then
  pokaz_error "Unrecognized options:"
  for option in "${unrecognized_options[@]}"
  do
    pokaz_error "$option"
  done
  exit 1
fi


# ---- Check of missimg parameters


check_missing_variable "path_out"
check_missing_variable "path_in"
check_missing_variable "ref_name"
check_missing_variable "n_chr_ref"
check_missing_variable "n_chr_query"


# ---- if some parameters zre not given - set the default value

# Basic parameters

start_step="${start_step:-100}"  # Starting step
cores="${cores:-1}"  # Number of cores
p_ident="${p_ident:-85}"  
part_len="${part_len:-5000}"  
all_cmp="${all_cmp:-T}"
sort_chr_len="${sort_chr_len:-F}"
filter_rep="${filter_rep:-0}"
flag_rev="${flag_rev:-0}"

acc_anal="${acc_anal:-NULL}"   # Set of accessions to analyse
comb_anal="${comb_anal:-NULL}"   # Set of combinations to analyse

# Rename the reference genome prefix
# Rename the reference, it sould not contain any '_' symbol, because it is used later for splitting
# ref_name_true=${ref_name}
# ref_name=${ref_name//_/$'-'}


#---- Paths
# Required

path_in=$(add_symbol_if_missing "$path_in" "/")
path_out=$(add_symbol_if_missing "$path_out" "/")

# Could be defined
path_chr_acc="${path_chr_acc:-${path_out}chromosomes/}"
path_chr_acc=$(add_symbol_if_missing "$path_chr_acc" "/")

path_parts="${path_parts:-${path_out}parts/}"
path_parts=$(add_symbol_if_missing "$path_parts" "/")

path_chr_ref="${path_chr_ref:-${path_chr_acc}}"
path_chr_ref=$(add_symbol_if_missing "$path_chr_ref" "/")

path_consensus="${path_consensus:-${path_out}consensus/}"
path_consensus=$(add_symbol_if_missing "$path_consensus" "/")
# echo "consensus${path_consensus}"
if [ ! -d "$path_consensus" ]; then
    mkdir -p "$path_consensus"
fi

# New paths
path_blast_parts=${path_out}blast_parts_${ref_name}/
path_alignment=${path_out}alignments_${ref_name}/
path_gaps=${path_out}blast_gaps_${ref_name}/


# Path with stages
path_flags="${path_out}.flags/"
if [ ! -d "$path_flags" ]; then
    mkdir -p "$path_flags"
fi

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
make_dir ${path_logs}

# Path for logs of this script
path_log_ref="${path_logs}pangen_ref_${ref_name}/"
make_dir ${path_log_ref}

# File with steps logs
file_log_ref="${path_log_ref}steps.log"
> "${file_log_ref}"

# ----------------------------------------------------------------------------
#           ATTENTION MESSAGES
# ----------------------------------------------------------------------------

# Attention: If you use -one2one option: please be sure that all chromosomes in files are sorted in the same order

# if(!sort.by.lengths){
#   msg = 'If you use -one2one option: please be sure that all chromosomes in files are sorted in the same order' # or use \"-s T\" flag'
#   pokazAttention(msg)
# } else {
#   msg = 'Chromosomes will be sorted by their length'
#   pokazAttention(msg)
# }


# ----------------------------------------------------------------------------
#           MAIN PIPELINE
# ----------------------------------------------------------------------------

step_num=1
# ----------------------------------------------
# Split quiery fasta into chromosomes

if [ $start_step -le ${step_num} ] || [ ! -f "$path_flags/step${step_num}_done" ]; then

    log_message 1 "$log_level" "$file_log_ref" pokaz_stage "Step ${step_num}. Genomes into chromosomes."

    rm -rf ${path_chr_acc}

    path_log_step="${path_log_ref}step${step_num}_query_01/"

    Rscript pangen/query_01_to_chr.R --n.chr ${n_chr_query}  \
            --path.in ${path_in} --path.out ${path_chr_acc} --sort ${sort_chr_len} --cores ${cores} --acc.anal ${acc_anal} \
            --path.log ${path_log_step} --log.level ${log_level}
    
    touch "$path_flags/step${step_num}_done"
    log_message 1 "$log_level" "$file_log_ref" pokaz_message "Step is done."
fi

((step_num = step_num + 1))

# ----------------------------------------------
# Split quiery chromosomes into parts
if [ $start_step -le ${step_num} ] || [ ! -f "$path_flags/step${step_num}_done" ]; then

    log_message 1 "$log_level" "$file_log_ref" pokaz_stage "Step ${step_num}. Chromosomes into parts."

    rm -rf ${path_parts}

    path_log_step="${path_log_ref}step${step_num}_query_02/"

    Rscript pangen/query_02_to_parts.R --n.chr ${n_chr_query}  --path.chr  ${path_chr_acc}  \
            --path.parts ${path_parts} --part.len $part_len --cores ${cores} \
            --filter_rep ${filter_rep} --rev ${flag_rev} \
            --path.log ${path_log_step} --log.level ${log_level}

    touch "$path_flags/step${step_num}_done"
    log_message 1 "$log_level" "$file_log_ref" pokaz_message "Step is done."
fi

((step_num = step_num + 1))


# ----------------------------------------------
# Split reference fasta into chromosomes if additionally needed

if [[ "${path_chr_acc}" != "$path_chr_ref" ]]; then

    if [ $start_step -le ${step_num} ] || [ ! -f "$path_flags/step${step_num}_done_${ref_name}" ]; then
        log_message 1 "$log_level" "$file_log_ref" pokaz_stage "Step ${step_num}. Reference genome into chromosomes."

        file_acc_ref=${path_chr_acc}ref_acc.txt
        echo "${ref_name}" > ${file_acc_ref}
        Rscript pangen/query_01_to_chr.R --n.chr ${n_chr_ref}  \
                --path.in ${path_chr_ref} --path.out ${path_chr_acc}   \
                --cores ${cores} --acc.anal ${file_acc_ref}

        rm ${file_acc_ref}

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

    # ---- Cleanup the blast folder ----

    rm -rf ${path_blast_parts}

    # ---- Create a database on the reference genome ----

    log_message 2 "$log_level" "$file_log_ref" pokaz_message "Create BLAST databases."

    for file in ${path_chr_acc}${ref_name}_chr*fasta ; do
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
    ./pangen/query_03_blast_parts.sh -path_ref ${path_chr_acc} -path_parts ${path_parts} \
            -path_result ${path_blast_parts} -ref_name ${ref_name} \
            -all_vs_all ${all_cmp} -p_ident ${p_ident} -cores ${cores} -log_path ${path_log_step}

    touch "$path_flags/step${step_num}_done_${ref_name}"
    log_message 1 "$log_level" "$file_log_ref" pokaz_message "Step is done."
fi

((step_num = step_num + 1))

# ----------------------------------------------
# First round of alignments
if [ $start_step -le ${step_num} ] || [ ! -f "$path_flags/step${step_num}_done_${ref_name}" ]; then

    log_message 1 "$log_level" "$file_log_ref" pokaz_stage "Step ${step_num}. Alignment-1: Remaining syntenic (major) matches."

    rm -rf ${path_alignment}
    rm -rf ${path_gaps}

    path_log_step="${path_log_ref}step${step_num}_synteny_01/"
    
    Rscript pangen/synteny_01_majoir.R --path.blast ${path_blast_parts} --path.aln ${path_alignment} \
            --ref ${ref_name}   \
            --path.gaps  ${path_gaps} --path.chr ${path_chr_acc} \
            --cores ${cores} \
            --path.log ${path_log_step} --log.level ${log_level}

    touch "$path_flags/step${step_num}_done_${ref_name}"
    log_message 1 "$log_level" "$file_log_ref" pokaz_message "Step is done."

    # # If the first round of alignment didn't have any errors - remove the blast which was needed for it
    # rm -rf ${path_blast_parts}
fi

((step_num = step_num + 1))

# ----------------------------------------------
# Blast regions between synteny blocks
if [ $start_step -le ${step_num} ] || [ ! -f "$path_flags/step${step_num}_done_${ref_name}" ]; then

    log_message 1 "$log_level" "$file_log_ref" pokaz_stage "Step ${step_num}. BLAST of gaps between syntenic matches."

    # Logs for the BLAST
    path_log_step="${path_log_ref}step${step_num}_synteny_02/"
    make_dir ${path_log_step}

    # Run BLAST for gaps
    ./pangen/synteny_02_blast_gaps.sh -path_gaps ${path_gaps} -cores ${cores} #-log_path ${path_log_step}

    touch "$path_flags/step${step_num}_done_${ref_name}"
    log_message 1 "$log_level" "$file_log_ref" pokaz_message "Step is done."
fi

((step_num = step_num + 1))

# ----------------------------------------------
# Second round of alignments
if [ $start_step -le ${step_num} ] || [ ! -f "$path_flags/step${step_num}_done_${ref_name}" ]; then

    log_message 1 "$log_level" "$file_log_ref" \
        pokaz_stage "Step ${step_num}. Alignment-2: Fill the gaps between synteny blocks."

    Rscript pangen/synteny_03_merge_gaps.R --ref ${ref_name} \
            --path.aln ${path_alignment}  --path.chr ${path_chr_acc}\
            --path.gaps ${path_gaps}   \
            --cores ${cores}

    # # If the second round of alignment didn't have any errors - remove the blast which was needed for it
    # rm -rf ${path_gaps}
    # ls ${path_alignment}*maj*
    # rm -rf ${path_alignment}*maj*

    touch "$path_flags/step${step_num}_done_${ref_name}"
    log_message 1 "$log_level" "$file_log_ref" \
        pokaz_message "Done."
fi

((step_num = step_num + 1))

# ----------------------------------------------
# Creaete a consensus
if [ $start_step -le ${step_num} ] || [ ! -f "$path_flags/step${step_num}_done_${ref_name}" ]; then

    log_message 1 "$log_level" "$file_log_ref" \
        pokaz_stage "Step ${step_num}. Combine reference-based alignments by chromosomes."

    Rscript pangen/comb_01_one_ref.R --path.cons ${path_consensus} --path.aln ${path_alignment} \
    --pref ${ref_name}  --cores ${cores} --path.chr ${path_chr_acc}

    touch "$path_flags/step${step_num}_done_${ref_name}"
    log_message 1 "$log_level" "$file_log_ref" \
        pokaz_message "Done."
fi

((step_num = step_num + 1))


