#tair="../ref/"
#parts="../parts/"
#blastres="../blast_res_tair/"
#ref_pref="TAIR10_chr"
#ref_type='fas'


echo -e "\e[38;2;52;252;252m* BLAST of parts on the reference genome\e[0m"


print_usage() {
  echo "-path_ref"
  echo "-path_parts"
  echo "-path_result"
  echo "-ref_pref"
  echo "-ref_type"
  echo "-all_vs_all"
  echo "-p_ident"
  echo "-cores"
  echo "-penalty"
  echo "-gapopen"
  echo "-gapextend"
  echo "-max_hsps"
}

while [ $# -gt 0 ]
do
    case $1 in
    # for options with required arguments, an additional shift is required
    -path_ref) tair=$2; shift ;;
    -path_parts) parts=$2; shift ;;
    -path_result) blastres=$2; shift ;;
    -ref_pref) ref_pref=$2; shift ;;
    -ref_type) ref_type=$2; shift ;;
    -all_vs_all) all_vs_all=$2; shift ;;
    -p_ident) p_ident=$2; shift ;;
    -cores) cores=$2; shift ;;
    -penalty) penalty=$2; shift ;;
    -gapopen) gapopen=$2; shift ;;
    -gapextend) gapextend=$2; shift ;;
    -max_hsps) max_hsps=$2; shift ;;
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


#echo $blastres
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

    p_filename=$(basename "$part_file" .fasta)
    p_prefix=${p_filename%_*}
    part_chr=${p_filename##*_}

    r_filename=$(basename "$ref_file" .fasta)
    r_prefix=${r_filename%_*}
    ref_chr=${r_filename##*chr}

    outfile=${blastres}${p_filename}_${ref_chr}.txt

    if [[ "$p_prefix" == "$r_prefix" ]] || { [[ "$part_chr" != "$ref_chr" ]] && [[ ${all_vs_all} == "F" ]]; } || [[ -f "$outfile" ]]; then
        return
    fi

    blastn -db ${ref_file} -query ${part_file} -out ${outfile} \
           -outfmt "7 qseqid qstart qend sstart send pident length qseq sseq sseqid" \
           -perc_identity ${p_ident} -penalty $penalty -gapopen $gapopen -gapextend $gapextend -max_hsps $max_hsps > /dev/null 2>> log_err.txt 
}

export -f run_blast

# Run the parallel
parallel run_blast ::: ${parts}*.fasta ::: $tair${ref_pref}*.$ref_type ::: $blastres ::: $p_ident ::: $penalty ::: $gapopen ::: $gapextend ::: $max_hsps ::: $all_vs_all

echo "  Done!"
