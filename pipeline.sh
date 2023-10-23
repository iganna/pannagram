# This escript performs all stages of the alignment of genomes, 
# when the reference genome is already identified

# ==============================================================================
#      PARAMETERS

# ./pipeline.sh -pref_global 'rhiz' -ref_pref 'ref_1021' -n_chr_ref 1 -path_in '../rhizobia/' -n_chr_query 1 -sort_chr_len T 

# ./pipeline.sh -pref_global 'ly' -ref_pref '0' -path_chr_ref "../pb_chromosomes/" -n_chr_ref 5 -path_in '../lyrata/' -n_chr_query 8


# ./pipeline.sh -pref_global 'ly2' -ref_pref '0' -path_chr_ref "../pb_chromosomes/" -n_chr_ref 5 -path_in '../lyrata/' -n_chr_query 8

# ./pipeline.sh -pref_global 'toy' -ref_pref '0'  -n_chr_ref 5 -path_in '../pb_genomes/' -n_chr_query 5 -all_cmp F -acc_anal 'acc_analysis.txt'


print_usage() {
  echo "-pref_global"
  
  echo "-ref_pref"
  echo "-path_chr_ref"
  echo "-n_chr_ref"

  echo "-path_in"
  echo "-n_chr_query"

  echo "-path_parts"
  echo "-path_chr_acc"

  echo "-path_consensus"

  echo "-sort_chr_len"
  echo "-part_len"  

  echo "-all_cmp"
  echo "-p_ident"

  echo "-cores"
}

unrecognized_options=()

while [ $# -gt 0 ]
do
    case $1 in
    # for options with required arguments, an additional shift is required
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


acc_anal="${acc_anal:-NULL}"

cores="${cores:-30}"
p_ident="${p_ident:-85}"
part_len="${part_len:-5000}"
all_cmp="${all_cmp:-T}"
sort_chr_len="${sort_chr_len:-F}"

path_chr_acc="${path_chr_acc:-../${pref_global}_chromosomes/}"
path_parts="${path_parts:-../${pref_global}_parts/}"

path_chr_ref="${path_chr_ref:-${path_chr_acc}}"

path_consensus="${path_consensus:-../${pref_global}_consensus/}"


# echo $pref_global
# echo $ref_pref
# echo $path_chr_ref
# echo $n_chr_ref
# echo $path_in
# echo $n_chr_query
# echo $part_len
# echo $all_cmp
# echo $p_ident
# echo $cores
# echo $sort_chr_len
# echo ${acc_anal}

# exit 1


# ==============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
#trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

trap 'catch $?' EXIT
catch() {
if [ $1 -ne 0 ]; then
        echo "\"${last_command}\" command filed with exit code $1."
fi
#  echo "\"${last_command}\" command filed with exit code $1."
}

# ==============================================================================
# Split quiery fasta into chromosomes
Rscript query_to_chr.R -n ${n_chr_query} -t fasta --path.in ${path_in} --path.out ${path_chr_acc} -s ${sort_chr_len} -c ${cores} --acc.anal ${acc_anal}
# Split quiery chromosomes into parts
Rscript query_to_parts.R -n ${n_chr_query} -t fasta --path.chr  ${path_chr_acc} --path.parts ${path_parts} --part.len $part_len -c ${cores}


# Rename the reference, it sould not contain any '_' symbol, because it is used later for splitting
ref_pref=${ref_pref//_/$'-'}

# Create a database on the reference genome
for file in ${path_chr_ref}${ref_pref}_chr*fasta ; do
  # Check if the BLAST database files already exist
  if [ ! -f "${file}.nin" ]; then
      makeblastdb -in ${file} -dbtype nucl > /dev/null
  fi
done

# Blast parts on the reference genome
./blast_parts.sh -path_ref ${path_chr_ref} -path_parts ${path_parts} -path_result ../${pref_global}_blast_res_${ref_pref}/ \
 -ref_pref ${ref_pref}_chr -ref_type fasta -all_vs_all ${all_cmp} -p_ident ${p_ident} -cores 30


# First round of alignments
Rscript synteny_majoir.R --path.blast ../${pref_global}_blast_res_${ref_pref}/ --path.aln ../${pref_global}_alignments_${ref_pref}/ \
--type fasta --pref ${ref_pref} --path.ref  ${path_chr_ref}  \
--path.gaps ../${pref_global}_gaps_${ref_pref}/  --path.query ${path_chr_acc} \
--n.chr.ref ${n_chr_ref} --n.chr.acc ${n_chr_query}  --all.vs.all ${all_cmp} -c ${cores}

# # If the first round of alignment didn't have any errors - remove the blast which was needed for it
rm -rf ../${pref_global}_blast_res_${ref_pref}/

# Blast regions between synteny blocks
./blast_gaps.sh -path_gaps ../${pref_global}_gaps_${ref_pref}/ 

# # Second round of alignments
Rscript synteny_gaps.R --path.blast ../${pref_global}_blast_res_${ref_pref}/ --path.aln ../${pref_global}_alignments_${ref_pref}/ \
--type fasta --pref ${ref_pref} --path.ref  ${path_chr_ref}  \
--path.gaps ../${pref_global}_gaps_${ref_pref}/  --path.query ${path_chr_acc} \
--n.chr.ref ${n_chr_ref} --n.chr.acc ${n_chr_query}  --all.vs.all ${all_cmp} -c ${cores}

# # If the second round of alignment didn't have any errors - remove the blast which was needed for it
rm -rf ../${pref_global}_gaps_${ref_pref}/

# # -----------------------------------
# # Creaete a consensus


# if [ ! -d "${path_consensus}" ]; then
#   mkdir ${path_consensus}
# fi

Rscript  comb_alignment.R --path.cons ${path_consensus} --path.aln ../${pref_global}_alignments_${ref_pref}/ \
--type fasta --pref ${ref_pref} --path.ref  ${path_chr_ref}  \
--n.chr.ref ${n_chr_ref} --n.chr.acc ${n_chr_query}  --all.vs.all ${all_cmp} -c ${cores}