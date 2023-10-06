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
max_hsps="${max_hsps:-30}"
cores="${cores:-30}"


#echo $blastres
mkdir -p $blastres

pids=""


for partfile in ${parts}*.fasta; do

  for tairfile in $tair${ref_pref}*.$ref_type; do
  #echo $partfile $tairfile 
  if [[ "$tairfile" == "$partfile" ]] 
  then
      continue
  fi

	b=$(basename $partfile)
    	name="${b%.fasta}"
	part_chr="${name: -1}"
	t=$(basename $tairfile)
	t="${t#${ref_pref}}"
	tair_chr="${t%.${ref_type}}"
	outfile=${blastres}${name}_${tair_chr}.txt
        
	if [[ "$part_chr" != "$tair_chr" ]] && [[ ${all_vs_all} == "F" ]]
  then
	    continue
	fi

	if [[ ! -f "$outfile" ]]
	then
    	
#    echo "BLAST is running with output ${outfile} with reference ${tairfile}"
    blastn -db ${tairfile} -query ${partfile} -out ${outfile} \
           -outfmt "7 qseqid qstart qend sstart send pident length qseq sseq sseqid" \
           -perc_identity ${p_ident} -penalty $penalty -gapopen $gapopen -gapextend $gapextend -max_hsps $max_hsps & 

 	  blast_number=$(pgrep -c blastn)
    pids="$pids $!"
#    echo $blast_number
 # else
#    echo "file ${outfile} exists, no blast running"
  fi


  blast_number=$(pgrep -c blastn)
  if (( ${blast_number} > $cores )); then
      wait -n
  fi
  done

done

echo "  Waiting till all BLAST of chromosomal parts is over.."
wait $pids
echo "  Done!"
