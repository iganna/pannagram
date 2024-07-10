# This escript performs all stages of the alignment of genomes, 
# when the reference genome is already identified

# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

source utils/error_block.sh

# ----------------------------------------------------------------------------
#             USAGE
# ----------------------------------------------------------------------------


#./pangen_ref.sh -path_out '../pan_test/tom/' -ref_name '0'  -n_chr_ref 5 -path_in '../pb_updated/' -n_chr_query 5 -all2all F -acc_anal 'acc_tom.txt'
#./pangen_ref.sh -path_out '../pan_test/tom/' -ref_name '6046-v1.1'  -n_chr_ref 5 -path_in '../pb_updated/' -n_chr_query 5 -all2all F -acc_anal 'acc_tom.txt'

#./pangen_ref.sh -path_out ../pan_test/ly_th/ -ref_name 0 -path_ref ../pan_test/p27/chromosomes/ -n_chr_ref 5 -path_in ../lyrata/ -n_chr_query 8 -all2all T -cores 30


# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source utils/utils_bash.sh

# Function to display help message
print_usage() {
    pokaz_help

    cat << EOF
Usage: ${0##*/} -ref REF_NAME -path_in INPUT_FOLDER -path_out OUTPUT_FOLDER
                -nchr_ref N_CHR_REF -nchr_acc N_CHR_QUERY
                [-path_ref PATH_REF] [-path_chrom PATH_CHROM] 
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
        -nchr_acc N_CHR_QUERY     Number of chromosomes in the query genome.

    * Paths (optional):
        -path_ref PATH_REF          Path where the reference genome is stored. 
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
    ${0##*/} -path_out 'output_folder' -path_in 'input_genomes' -ref '0' -nchr_acc 5 -nchr_ref 5

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

        -accessions) acc_anal=$2;    shift 1 ;;  # file with accessions to analyse
        -combinations) comb_anal=$2; shift 1 ;;  # file with chromosomal combinations to analyse: first column - query, second column - reference(base)
    
        *) print_usage
            unrecognized_options+=("$1"); shift 1 ;;
    esac
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
one2one="${one2one:-F}"
sort_chr_len="${sort_chr_len:-F}"
purge_reps="${purge_reps:-F}"
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
path_chrom="${path_chrom:-${path_out}chromosomes/}"
path_chrom=$(add_symbol_if_missing "$path_chrom" "/")

path_parts="${path_parts:-${path_out}parts/}"
path_parts=$(add_symbol_if_missing "$path_parts" "/")

path_ref="${path_ref:-${path_in}}"
path_ref=$(add_symbol_if_missing "$path_ref" "/")


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
#           STAGES
# ----------------------------------------------------------------------------



# ----------------------------------------------------------------------------
#           LOGS
# ----------------------------------------------------------------------------

source utils/logging_ref.sh

# ----------------------------------------------------------------------------
#           MAIN PIPELINE
# ----------------------------------------------------------------------------

# ----------------------------------------------
# Run the pangen_pre


if [ "$one2one" = "F" ]; then
    one2one_option=" "
else
    one2one_option=" -one2one "
fi

if [ "$purge_reps" = "F" ]; then
    purge_reps_option=" "
else
    purge_reps_option=" -purge_reps "
fi

if [ "$purge_reps" = "F" ]; then
    purge_reps_option=" "
else
    purge_reps_option=" -purge_reps "
fi

./pangen_pre.sh -path_in ${path_in} -path_out ${path_out} \
        -ref ${ref_name} -path_ref ${path_ref} \
        -path_chrom ${path_chrom} -path_parts ${path_parts} -path_cons ${path_consensus} \
        -nchr_ref ${nchr_ref} -nchr_acc ${nchr_acc} \
        -part_len ${part_len}  -p_ident ${-p_ident} \
        ${one2one_option} ${purge_reps_option} \
        -accessions ${acc_anal} -combinations) ${comb_anal}

# ----------------------------------------------
# Fing the number of last successful stage        


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
            --path.aln ${path_alignment}  --path.chr ${path_chrom}\
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
    --pref ${ref_name}  --cores ${cores} --path.chr ${path_chrom}

    touch "$path_flags/step${step_num}_done_${ref_name}"
    log_message 1 "$log_level" "$file_log_ref" \
        pokaz_message "Done."
fi

((step_num = step_num + 1))


