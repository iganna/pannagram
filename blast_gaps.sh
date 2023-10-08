echo -e "\e[38;2;52;252;252m* BLAST of gaps between syn blocks\e[0m"


print_usage() {
  echo "-path_gaps"
}

while [ $# -gt 0 ]
do
    case $1 in
    # for options with required arguments, an additional shift is required
    -path_gaps) path_gaps=$2; shift ;;
    *) print_usage
       echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    esac
    shift
done


for file in ${path_gaps}*_base*.txt; do

tmp="${file}.nhr"
if [ -f $tmp ]; then
    continue
fi
#echo "DTABASE ${file}"
    makeblastdb -in $file -dbtype nucl &> /dev/null
done



cores=30

for file_q in ${path_gaps}*_query*.txt; do
    
    file_b="${file_q/query/base}"
    file_out="${file_q/query/out}"

if [ ! -f ${path_gaps}${file_b} ]
then
#	echo "File not found"
continue
fi

#    echo $file_b
    blastn -db  ${path_gaps}${file_b} -query  ${path_gaps}${file_q} \
     -out  ${path_gaps}${file_out} \
     -outfmt "7 qseqid qstart qend sstart send pident length qseq sseq sseqid" &  # -perc_identity 85 -penalty -2 -gapopen 5 -gapextend 2
     #-max_hsps 5  \
    
   pids="$pids $!"
blast_number=$(pgrep -c blastn)
#echo $blast_number
#    background=( $(jobs -p) )
    if (( ${blast_number} > $cores )); then
        wait -n
    fi
done

wait $pids

echo "  Done!"
#Rscript ../pipeline/get_alignments_fixed.R --path.blast ../blast_res_ref/ --path.aln ../alignments_ref/ --type fasta --pref TAIR10 --path.ref  ../ref/ --path.gaps ../gaps_ref/
