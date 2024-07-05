# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

source utils/error_block.sh

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source utils/utils_bash.sh

print_usage() {
  echo "-path_ref"
  echo "-path_parts"
  echo "-path_result"
  echo "-ref_name"
  echo "-all_vs_all"
  echo "-p_ident"
  echo "-cores"
  echo "-penalty"
  echo "-gapopen"
  echo "-gapextend"
  echo "-max_hsps"
  echo "-log_path"
}

# ----------------------------------------------------------------------------
#                 PARAMETERS
# ----------------------------------------------------------------------------

all_vs_all=""

while [ $# -gt 0 ]
do
    case $1 in
    # for options with required arguments, an additional shift is required
    -path_ref) path_ref=$2; shift ;;
    -path_parts) parts=$2; shift ;;
    -path_result) blastres=$2; shift ;;
    -ref_name) ref_name=$2; shift ;;
    -all_vs_all) all_vs_all=$2; shift ;;
    -p_ident) p_ident=$2; shift ;;
    -cores) cores=$2; shift ;;
    -penalty) penalty=$2; shift ;;
    -gapopen) gapopen=$2; shift ;;
    -gapextend) gapextend=$2; shift ;;
    -max_hsps) max_hsps=$2; shift ;;
    -log_path) log_path=$2; shift ;;
    *) print_usage
       echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    esac
    shift
done



# -penalty -2 -gapopen 10 -gapextend 2 -max_hsps 5

penalty="${penalty:--2}"
gapopen="${gapopen:-10}"
gapextend="${gapextend:-2}"
max_hsps="${max_hsps:-1}"
cores="${cores:-30}"

if [ -z ${all_vs_all} ]; then
    echo "Error: Variable 'all_vs_all' is not set."
    exit 1
fi


# ----------------------------------------------------------------------------
#                 MAIN
# ----------------------------------------------------------------------------

mkdir -p $blastres

# BLAST-search function
run_blast() {
    part_file=$1
    ref_file=$2
    blastres=$3
    p_ident=$4
    penalty=$5
    gapopen=$6
    gapextend=$7
    max_hsps=$8
    all_vs_all=$9
    log_path=${10}

    p_filename=$(basename "$part_file" .fasta)
    p_prefix=${p_filename%_*}
    part_chr=${p_filename##*chr}

    r_filename=$(basename "$ref_file" .fasta)
    r_prefix=${r_filename%_*}
    ref_chr=${r_filename##*chr}

    # Check -one2one
    if [[ "$p_prefix" == "$r_prefix" ]] || { [[ "$part_chr" != "$ref_chr" ]] && [[ ${all_vs_all} == "F" ]]; } || [[ -f "$outfile" ]]; then
        return
    fi

  	p_filename=$(echo "$p_filename" | sed 's/_chr\(.*\)$/_\1/')
    outfile=${blastres}${p_filename}_${ref_chr}.txt

    # Create a log file
    if [ -d "$log_path" ]; then
        file_log="${log_path}${p_filename}_${ref_chr}.log"
        > "$file_log"
    else
        file_log="/dev/null"
    fi

    # Run BLAST
    blastn -db "${ref_file}" -query "${part_file}" -out "${outfile}" \
           -outfmt "6 qseqid qstart qend sstart send pident length qseq sseq sseqid" \
           -perc_identity "${p_ident}" -penalty "$penalty" -gapopen "$gapopen" -gapextend "$gapextend" \
           -max_hsps "$max_hsps" >> "$file_log" 2>&1 # -word_size 50 

    if [ -d "$log_path" ]; then
        echo "Done." >> "$file_log"
    fi

}

export -f run_blast

# Run BLAST in parallel
parallel -j $cores run_blast ::: ${parts}*.fasta ::: $path_ref${ref_name}_chr*.fasta ::: $blastres ::: $p_ident ::: $penalty ::: $gapopen ::: $gapextend ::: $max_hsps ::: $all_vs_all ::: ${log_path}

