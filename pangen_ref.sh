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


#./work.sh -pref_global '../pan_test/tom/' -ref_pref '0'  -n_chr_ref 5 -path_in '../pb_updated/' -n_chr_query 5 -all_cmp F -acc_anal 'acc_tom.txt'
#./work.sh -pref_global '../pan_test/tom/' -ref_pref '6046-v1.1'  -n_chr_ref 5 -path_in '../pb_updated/' -n_chr_query 5 -all_cmp F -acc_anal 'acc_tom.txt'


# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source utils/utils_bash.sh

# Function to display help message
print_usage() {

pokaz_help

    cat << EOF
Usage: ${0##*/} [-pref_global PREFIX] [-ref_pref REF_PREFIX] [-path_chr_ref PATH_CHR_REF]
                [-n_chr_ref N_CHR_REF] [-path_in PATH_IN] [-n_chr_query N_CHR_QUERY]
                [-path_parts PATH_PARTS] [-path_chr_acc PATH_CHR_ACC] [-path_consensus PATH_CONSENSUS]
                [-sort_chr_len SORT_CHR_LEN] [-part_len PART_LEN] [-all_cmp ALL_CMP]
                [-p_ident P_IDENT] [-fasta_type FASTA_TYPE] [-acc_anal ACC_ANAL] [-cores CORES]

This script performs genomic analysis with several options to specify inputs, paths, and parameters.

Options:
    -pref_global PREFIX        Global prefix for the analysis.
    -ref_pref REF_PREFIX       Reference genome prefix: REF_PREFIX.fasta
    -path_in PATH_IN           Path to directory with genomes.
    -n_chr_query N_CHR_QUERY   Number of chromosomes in the query genome.
    -n_chr_ref N_CHR_REF       Number of chromosomes in the reference genome.
    -path_chr_ref PATH_CHR_REF Path to reference chromosomes.
    -path_parts PATH_PARTS     Path to parts directory.
    -path_chr_acc PATH_CHR_ACC Path to chromosome accession files.
    -path_consensus PATH_CONSENSUS
                                Path to consensus directory.
    -sort_chr_len SORT_CHR_LEN Flag to sort chromosomes by length.
    -part_len PART_LEN         Length of part, which should be searched on the first step (default: 5000).
    -all_cmp ALL_CMP           Compare all vs all (default: "T").
    -p_ident P_IDENT           Percentage identity threshold.
    -fasta_type FASTA_TYPE     Type of FASTA file used.
    -acc_anal ACC_ANAL         File with accessions to analyse. 
                               Accessions should be in rows.
    -cores CORES               Number of cores for parallel processing.

Examples:
    ${0##*/} -pref_global 'output_folder' -ref_pref '0' -path_in 'thaliana_genomes' -n_chr_query 5 -n_chr_ref 5

EOF
}


# ----------------------------------------------------------------------------
#            PARAMETERS
# ----------------------------------------------------------------------------


unrecognized_options=()

while [ $# -gt 0 ]
do
    case $1 in
  -h | --help )  print_usage
                       exit
                       ;;
    # for options with required arguments, an additional shift is required

  -s | -stage ) start_step="$2"; shift ;;                     
  -pref_global) pref_global=$2; shift ;;
	
	-ref_pref) ref_pref=$2; shift ;;
	-path_chr_ref) path_chr_ref=$2; shift ;;  # dont provide if it's the same as the path with query genomes
	-n_chr_ref) n_chr_ref=$2; shift ;;

	-path_in) path_in=$2; shift ;;
	-n_chr_query) n_chr_query=$2; shift ;;

  -path_parts) path_parts=$2; shift ;;
  -path_chr_acc) path_chr_acc=$2; shift ;;

	-path_consensus) path_consensus=$2; shift ;;

	-sort_chr_len) sort_chr_len=$2; shift ;;
	-part_len) part_len=$2; shift ;;  # has a default value 5000

	-all_cmp) all_cmp=$2; shift ;;   # has a defaul value "T": compare all vs all
	-p_ident) p_ident=$2; shift ;;

  -fasta_type) fasta_type=$2; shift ;;

  -acc_anal) acc_anal=$2; shift ;;  # accessions to analyse

	-cores) cores=$2; shift ;;
    
    *) print_usage
       # echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
      unrecognized_options+=("$1")
      shift
      ;;
    esac
    shift
done


# Вывод нераспознанных параметров
if [[ ${#unrecognized_options[@]} -gt 0 ]]
then
  echo "Unrecognized options:"
  for option in "${unrecognized_options[@]}"
  do
    echo "$option"
  done
fi


# ---- Check of missimg parameters


check_missing_variable "pref_global"

check_missing_variable "path_in"
check_missing_variable "ref_pref"
check_missing_variable "n_chr_ref"
check_missing_variable "n_chr_query"



# ---- if some parameters zre not given - set the default value

# Basic parameters

start_step="${start_step:-1}"  # Starting step
cores="${cores:-1}"  # Number of cores
p_ident="${p_ident:-85}"  
part_len="${part_len:-5000}"  
all_cmp="${all_cmp:-T}"
sort_chr_len="${sort_chr_len:-F}"
fasta_type="${fasta_type:-fasta}"


acc_anal="${acc_anal:-NULL}"   # Set of accessions to analyse

# Rename the reference genome prefix
# Rename the reference, it sould not contain any '_' symbol, because it is used later for splitting
ref_pref_true=${ref_pref}
ref_pref=${ref_pref//_/$'-'}


#---- Paths
# Required

path_in=$(add_symbol_if_missing "$path_in" "/")
pref_global=$(add_symbol_if_missing "$pref_global" "/")


# Could be defined
path_chr_acc="${path_chr_acc:-${pref_global}chromosomes/}"
path_chr_acc=$(add_symbol_if_missing "$path_chr_acc" "/")

path_parts="${path_parts:-${pref_global}parts/}"
path_parts=$(add_symbol_if_missing "$path_parts" "/")

path_chr_ref="${path_chr_ref:-${path_chr_acc}}"
path_chr_ref=$(add_symbol_if_missing "$path_chr_ref" "/")

path_consensus="${path_consensus:-${pref_global}consensus/}"
path_consensus=$(add_symbol_if_missing "$path_consensus" "/")
# echo "consensus${path_consensus}"
if [ ! -d "$path_consensus" ]; then
    mkdir -p "$path_consensus"
fi

# New paths
path_blast_parts=${pref_global}blast_parts_${ref_pref}/
path_alignment=${pref_global}alignments_${ref_pref}/
path_gaps=${pref_global}blast_gaps_${ref_pref}/


# Path with stages
path_flags="${pref_global}flags/"
if [ ! -d "$path_flags" ]; then
    mkdir -p "$path_flags"
fi


# ----------------------------------------------------------------------------
#           MAIN PIPELINE
# ----------------------------------------------------------------------------


# Split quiery fasta into chromosomes
if [ $start_step -le 1 ] && [ ! -f "$path_flags/step1_done" ]; then
    Rscript pangen/query_01_to_chr.R -n ${n_chr_query} -t ${fasta_type} --path.in ${path_in} --path.out ${path_chr_acc} -s ${sort_chr_len} -c ${cores} --acc.anal ${acc_anal}
    touch "$path_flags/step1_done"
fi


# Split quiery chromosomes into parts
if [ $start_step -le 2 ] && [ ! -f "$path_flags/step2_done" ]; then
    Rscript pangen/query_02_to_parts.R -n ${n_chr_query} -t ${fasta_type} --path.chr  ${path_chr_acc} --path.parts ${path_parts} --part.len $part_len -c ${cores}
    touch "$path_flags/step2_done"
fi


# Blast parts on the reference genome
if [ $start_step -le 3 ] && [ ! -f "$path_flags/step3_done" ]; then

    # Create a database on the reference genome
    for file in ${path_chr_ref}${ref_pref}_chr*${fasta_type} ; do
      # Check if the BLAST database files already exist
      if [ ! -f "${file}.nin" ]; then
          makeblastdb -in ${file} -dbtype nucl > /dev/null
      fi
    done

    # Blast parts on the reference genome
    ./pangen/query_03_blast_parts.sh -path_ref ${path_chr_ref} -path_parts ${path_parts} -path_result ${path_blast_parts} \
     -ref_pref ${ref_pref}_chr -ref_type ${fasta_type} -all_vs_all ${all_cmp} -p_ident ${p_ident} -cores ${cores}

    touch "$path_flags/step3_done"
fi


# First round of alignments
if [ $start_step -le 4 ] && [ ! -f "$path_flags/step4_done" ]; then
    Rscript pangen/synteny_01_majoir.R --path.blast ${path_blast_parts} --path.aln ${path_alignment} \
                        --type ${fasta_type} --pref ${ref_pref} --path.ref  ${path_chr_ref}  \
                        --path.gaps  ${path_gaps} --path.query ${path_chr_acc} \
                        --n.chr.ref ${n_chr_ref} --n.chr.acc ${n_chr_query}  --all.vs.all ${all_cmp} -c ${cores}

    touch "$path_flags/step4_done"

    # # If the first round of alignment didn't have any errors - remove the blast which was needed for it
    rm -rf ${path_blast_parts}
fi


# Blast regions between synteny blocks
if [ $start_step -le 5 ] && [ ! -f "$path_flags/step5_done" ]; then
    ./pangen/synteny_02_blast_gaps.sh -path_gaps ${path_gaps} -cores ${cores}
    touch "$path_flags/step5_done"
fi




# Second round of alignments
if [ $start_step -le 6 ] && [ ! -f "$path_flags/step6_done" ]; then

    Rscript pangen/synteny_03_merge_gaps.R --path.aln ${path_alignment} \
    --type ${fasta_type} --pref ${ref_pref} --path.ref  ${path_chr_ref}  \
    --path.gaps ${path_gaps}  --path.query ${path_chr_acc} \
    --n.chr.ref ${n_chr_ref} --n.chr.acc ${n_chr_query}  --all.vs.all ${all_cmp} -c ${cores}

    # # If the second round of alignment didn't have any errors - remove the blast which was needed for it
    rm -rf ${path_gaps}
    # ls ${path_alignment}*maj*
    rm -rf ${path_alignment}*maj*

    touch "$path_flags/step6_done"
fi



# -----------------------------------
# Creaete a consensus
if [ $start_step -le 7 ] && [ ! -f "$path_flags/step7_done" ]; then

    Rscript pangen/comb_01_one_ref.R --path.cons ${path_consensus} --path.aln ${path_alignment} \
    --type ${fasta_type} --pref ${ref_pref} --path.ref  ${path_chr_ref}  \
    --n.chr.ref ${n_chr_ref} --n.chr.acc ${n_chr_query}  --all.vs.all ${all_cmp} -c ${cores}

    touch "$path_flags/step7_done"
fi


