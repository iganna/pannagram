#!/bin/bash

# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

source utils/chunk_error_control.sh

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source utils/utils_bash.sh

show_help() {
    cat << EOF

╔════════════════════════════════════════╗
║   S e a r c h   f o r   s i m i l a r  ║
║            s e q u e n c e s           ║
╚════════════════════════════════════════╝

This script performs a BLAST search on a given FASTA file against a specified genome 
and processes the results based on similarity thresholds.

Usage: ${0##*/}  -in_seq FASTA_FILE -out OUTPUT_FILE
                [mode_option] [options]

Mode (at least one option is required:
    -on_seq SEQUENCE_FILE   Fasta-file containing sequences for comparison.
    -on_genome GENOME_FILE  Fasta-file containing genomes for comparison.
    -on_path SENOME_FOLDER  Path to the folder containing fasta-files with genomes.

Options:
    -in_seq FASTA_FILE     Path to the input FASTA file containing DNA sequences to be processed.
    -out OUTPUT_FILE       Path to the output file where the processed results will be saved.
    
    -sim SIMILARITY_CUTOFF (Optional) Similarity cutoff for sequence comparison.
    -afterblast            (Optional) Flag to process existing BLAST results.
    -keepblast             (Optional) Flag to keep intermediate BLAST results.
    -strandfree            (Optional) Use both strands for coverage. This option is used together with -on_seq.
    -h                     Show this help message and exit.

Examples:
    ${0##*/} -in_seq input.fasta -on_genome genome.fasta -out out/
    ${0##*/} -in_seq input.fasta -on_genome genome.fasta -out out_90.txt -sim 90 -keepblast
    ${0##*/} -in_seq input.fasta -on_genome genome.fasta -out out_95.txt -sim 95 -afterblast 

    ${0##*/} -in_seq input.fasta -on_seq sequences.fasta -out out/

    ${0##*/} -in_seq input.fasta -on_path folder_with_genomes -out out_90.txt -sim 90 -keepblast

EOF
}

# ----------------------------------------------------------------------------
#            PARAMETERS
# ----------------------------------------------------------------------------

#!/bin/bash

after_blast_flag=0
keep_blast_flag=0
use_strand=T

# Read arguments
# Read arguments
while [ "$1" != "" ]; do
    case $1 in
        -h | --help ) show_help; exit ;;
        -in_seq )    file_input=$2; shift 2 ;;
        -out )       output_pref=$2; shift 2 ;;
        -sim )       sim_threshold=$2; shift 2 ;;

        -on_seq )    file_seq=$2; shift 2 ;;
        -on_genome ) file_genome=$2; shift 2 ;;
        -on_path )   path_genome=$2; shift 2 ;;

        -afterblast ) after_blast_flag=1; shift ;;
        -keepblast )  keep_blast_flag=1; shift ;;

        -strandfree ) use_strand=F; shift ;;
        * ) echo "Invalid parameter: $1"; show_help; exit 1 ;;
    esac
done

# Ensure only one of -on_seq, -on_genome, -on_path is set
count=0
[ ! -z "$file_seq" ] && count=$((count + 1))
[ ! -z "$file_genome" ] && count=$((count + 1))
[ ! -z "$path_genome" ] && count=$((count + 1))

if [ $count -ne 1 ]; then
    echo "Error: You must specify exactly one of -on_seq, -on_genome, or -on_path."
    exit 1
fi

# Check if FASTA file parameter is provided
if [ -z "$file_input" ]; then
    echo "FASTA file not specified"
    exit 1
fi

# Check if the FASTA file exists
if [ ! -f "$file_input" ]; then
    echo "Input FASTA file not found: $file_input"
    exit 1
fi

# Check if output file parameter is provided
if [ -z "$output_pref" ]; then
    echo "Output file not specified"
    exit 1
fi

# Check if similarity threshold parameter is provided
if [ -z "$sim_threshold" ]; then
    sim_threshold=85
    pokaz_message "Similarity threshold not specified, default: ${sim_threshold}"
fi

# Your script code goes here

# ----------------------------------------------------------------------------
#            MAIN
# ----------------------------------------------------------------------------

# ---------------------------------------------
# Fix the output file result
if [[ "${output_pref}" == */ ]]; then
    if [ ! -d "${output_pref}" ]; then
      mkdir -p "${output_pref}"
    fi

    output_pref="${output_pref}simsearch"
    pokaz_message "Prefex for the output file was changed to ${output_pref}"
  
fi

# ---------------------------------------------
# Files for the blast

# Add all .fasta files from path_genome to db_files if path_genome is not empty
if [ ! -z "$path_genome" ]; then
    db_files=()
    for genome_file in "$path_genome"/*.fasta; do
        # Get the base name of the file without path and extension
        base_name=$(basename "$genome_file" .fasta)
        db_files+=("$base_name")
    done
fi

# Add file_seq to db_files if it's not empty
if [ ! -z "$file_seq" ]; then
    path_genome="$(dirname "$file_seq")/"
    db_files=($(basename "$file_seq" .fasta))
fi

# Add file_genome to db_files if it's not empty
if [ ! -z "$file_genome" ]; then
    path_genome="$(dirname "$file_genome")/"
    db_files=($(basename "$file_genome" .fasta))
fi

# ---------------------------------------------
# Run the pileline

for db_file in "${db_files[@]}"; do

    # ---------------------------------------------
    # Check if the BLAST database exists for the current file
    db_file_full="${path_genome}$db_file.fasta"
    if [ ! -f "${db_file_full}.nhr" ]; then
        pokaz_stage "Creating database for $db_file..."
        makeblastdb -in  ${db_file_full} -dbtype nucl > /dev/null
    fi

    # ---------------------------------------------
    # Run BLAST
    # Define the temporary file for storing BLAST results
    db_name=$(basename -- "$db_file")
    db_name="${db_name%.*}"
    blast_res="${output_pref}.${db_name}.blast.tmp"

    # Check if BLAST results should be used from an existing file
    if [ "$after_blast_flag" -eq 1 ]; then
        if [ ! -f "${blast_res}" ]; then
            pokaz_error "Blast results file not found: ${file_input}"
            exit 1
        fi
    else
        # Perform BLAST search
        pokaz_stage "BLAST search in $db_file..."
        blastn  -db ${db_file_full} \
                -query ${file_input} \
                -out ${blast_res} \
                -outfmt "6 qseqid qstart qend sstart send pident length sseqid qlen slen" 
                # -perc_identity ${sim_threshold}
    fi

    # Check if the BLAST results file is empty
    if [ ! -s "${blast_res}" ]; then
        pokaz_message "Blast result is empty for ${db_file}"
        continue
    fi
    
    # ---------------------------------------------
    # Proceed to similarity search
    pokaz_stage "Similarity search in ${db_file}, cutoff ${sim_threshold}..."

    # Determine if the search is on a set of sequences or a genome
    if [ -n "$file_seq" ]; then
        # On a set of sequences
        Rscript sim/sim_in_seqs.R \
                --in_file $file_input \
                --res $blast_res \
                --out ${output_pref}.${db_name}.rds \
                --sim $sim_threshold \
                --use_strand $use_strand \
                --db_file $db_file
    else
        # On a genome
        Rscript sim/sim_in_genome.R \
                --in_file $file_input \
                --res $blast_res \
                --out $output_pref.${db_name} \
                --sim $sim_threshold
    fi

    # Remove the BLAST temporary file if not needed
    if [ "$keep_blast_flag" -ne 1 ] && [ "$after_blast_flag" -ne 1 ]; then
        rm ${blast_res}
    fi

done

# Combine all files to the total count file
if [ ! -z "$path_genome" ]; then
    Rscript sim/sim_in_genome_combine.R  \
            --out $output_pref \
            --sim $sim_threshold
fi

pokaz_message "Done!"


