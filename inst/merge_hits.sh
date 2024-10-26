#!/bin/bash
INSTALLED_PATH=$(Rscript -e "cat(system.file(package = 'pannagram'))")
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
#             FUNCTIONS
# ----------------------------------------------------------------------------

# Function to display a short help message
print_usage() {
    cat << EOF

Usage:
${0##*/} -i F_GFF -g F_GENOME -o PATH_OUT \
		[-s SIMILARITY] [-d DISTANCE] [-p PATTERNS] [-n COPY_NUM] [-k] [-h]

Options:
    -h, --help         Display this help message and exit.
    
    -i, --file_gff F_GFF         Path to the input gff file with TEs.
    -g, --file_genome F_GENOME   Path to the genome file that will be used for the analysis.
    -o, --path_out PATH_OUT      Path to the directory where all results and intermediate files will be saved.
    
    -s, --sim  SIMILARITY        Similarity cutoff value used in the analysis. Default: 90.
    -c, --covegare  COVERAGE     Similarity cutoff value used in the analysis. Default: equal to sim.
    -d, --distance  DISTANCE     Distance between two hits. Default: 1000.
    -p, --patterns PATTERNS      Patterns of repeats to analyse. Default: LTR.
    -n, --copy_number COPY_NUM   Minimum number of copies per genome. Default: 4.
    -r, --max_rounds MAX_ROUND   Maximum number of rounds of merging. Default: 20.
    -k, --keepblast              Keep the intermediate BLAST files

EOF
}

# ----------------------------------------------------------------------------
#            PARAMETERS: parsing
# ----------------------------------------------------------------------------

unrecognized_options=()

sim_sutoff=90
distance=1000
patterns="LTR"
copy_number=4
max_rounds=20
keepblast=""


while [ $# -gt 0 ]
do
    # echo $1
    case $1 in
        -h | --help) print_usage; exit ;;
        
        # Required
        -i | --file_gff)     file_gff=$2;      shift 2 ;;  
		-g | --file_genome)  file_genome=$2;  shift 2 ;;  
        -o | --path_out)     path_out=$2;     shift 2 ;;  
		
		# Optional
		-s | --sim) 		 sim_sutoff=$2;    shift 2 ;;
        -c | --covegare)     covegare=$2;    shift 2 ;;
		-d | --distance) 	 distance=$2;      shift 2 ;;  
        -p | --patterns)     patterns=$2;      shift 2 ;;
		-n | --copy_number)  copy_number=$2;   shift 2 ;;
        -r | --max_rounds)   max_rounds=$2;    shift 2 ;;

        -k | --keepblast)   keepblast=' -keepblast ';    shift 1 ;;

        *) unrecognized_options+=("$1"); shift 1 ;;
    esac
done


# Output of Unrecognized Parameters
if [[ ${#unrecognized_options[@]} -gt 0 ]]; then
    print_usage
    echo "Unrecognized options:"
    for option in "${unrecognized_options[@]}"; do
        echo "$option"
    done
    exit 1
fi


# Check if coverage parameter is provided. If not - set qeual to sim
if [ -z "$coverage" ]; then
    coverage=${sim_threshold}
fi

# ----------------------------------------------------------------------------
#            PARAMETERS: checking
# ----------------------------------------------------------------------------


# Check for required parameters
if [ -z "$file_gff" ] || [ -z "$file_genome" ] || [ -z "$path_out" ]; then
    echo "Error: Missing required parameters (-i, -g, -o)."
    print_usage
    exit 1
fi

# Check the output directory
if [ -d "$path_out" ]; then   # Exist ?
    if [ "$(ls -A "$path_out")" ]; then  # Empty ?
        echo "Warning: The directory '$path_out' exists and is not empty."
    fi
else
    mkdir -p "$path_out"  

    # if the directory was not created successfully
    if [ ! -d "$path_out" ]; then
        echo "Error: Failed to create the directory '$path_out'."
        exit 1
    fi
fi

# Ensure $path_out ends with /
if [ "${path_out: -1}" != "/" ]; then
    path_out="$path_out/"
fi

# ----------------------------------------------------------------------------
#            MAIN
# ----------------------------------------------------------------------------

# ----------------------------------------
# Read the gff file, and get sequences

file_merged_seqs="${path_out}merged_seqs_1.fasta"

Rscript $INSTALLED_PATH/merge/merge_01_extract_hits.R \
        --file.gff=${file_gff} \
        --file.genome=${file_genome} \
        --file.seqs=${file_merged_seqs} \
        --patterns=${patterns} \
        --len.gap=${distance}

# # ----------------------------------------
# # Simrearch and Merge

file_merged_seqs_fixed="${path_out}merged_seqs_fixed.txt"

for ((i=1; i<=max_rounds; i++))
do
	file_merged_seqs="${path_out}merged_seqs_${i}.fasta"

	if [ ! -f "${file_merged_seqs}" ]; then
        break
    fi
    
    path_simsearch="${path_out}simseqrch_seqs_${i}/"

	# Run simsearch

	$INSTALLED_PATH/simsearch.sh \
    -in_seq ${file_merged_seqs}    \
    -on_genome ${file_genome} \
    -out "${path_out}simseqrch_seqs_${i}/" \
    -sim ${sim_sutoff} \
    -cov ${covegare} \
    ${keepblast}
    

	# Get Collapsed sequences - neighbours only

    file_merged_seqs_next="${path_out}merged_seqs_$((i+1)).fasta"

    file_cnt=$(find ${path_simsearch} -type f -name "*${sim_sutoff}_${covegare}.cnt")

    if [ -z "$file_cnt" ]; then
        echo "Error: No file ending with ${sim_sutoff}_${covegare}.cnt found." >&2
        exit 1
    elif [ $(echo "$file_cnt" | wc -l) -ne 1 ]; then
        echo "Error: More than one file matching the pattern *${sim_sutoff}_${covegare}.cnt found." >&2
        exit 1
    fi

    Rscript $INSTALLED_PATH/merge/merge_02_new_hits.R \
        --file.cnt ${file_cnt} \
        --file.genome ${file_genome} \
        --file.seqs ${file_merged_seqs_next} \
        --file.fix ${file_merged_seqs_fixed} \
        --copy.number=${copy_number}

done

echo ${file_merged_seqs_fixed}
if [ ! -s ${file_merged_seqs_fixed} ];  then
    exit 0
fi


file_fix_seqs="${path_out}seqs_fix.fasta"

Rscript $INSTALLED_PATH/merge/merge_03_get_hits.R \
    --path.out ${path_out} \
    --file.fix ${file_merged_seqs_fixed} \
    --file.fix.seqs=${file_fix_seqs}


if grep -q "^>" ${file_fix_seqs}; then
    $INSTALLED_PATH/simsearch.sh \
        -in_seq ${file_fix_seqs}    \
        -on_genome ${file_genome} \
        -out "${path_out}simseqrch_seqs_fix/" \
        -sim ${sim_sutoff} \
        -cov ${covegare} \
        ${keepblast}
else
    echo "Nothing to merge"
fi


Rscript $INSTALLED_PATH/merge/merge_04_visualisation.R \
        --path.out  ${path_out} \
        --file.genome ${file_genome} \
        --file.gff=${file_gff} \
        --plot FALSE









