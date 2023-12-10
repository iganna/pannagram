#!/bin/bash



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

source utils_bash.sh

show_help() {
    cat << EOF
Usage: ${0##*/} [-h] [-in FASTA_FILE] [-genome GENOME_FILE] [-out OUTPUT_FILE] 
                [-sim SIMILARITY_THRESHOLD] [-afterblast] [-keepblast]

This script performs a BLAST search on a given FASTA file against a specified genome 
and processes the results based on similarity thresholds.

    -h, --help              Display this help and exit.
    -in FASTA_FILE          Specify the input FASTA file.
    -out OUTPUT_FILE        Specify the output file for results.
    -genome GENOME_FILE     Specify the genome file to create BLAST database.    
    -sim SIMILARITY_THRESHOLD 
                            Set the similarity threshold for BLAST (default: 85).
    -afterblast             Use this flag to process existing BLAST results.
    -keepblast              Use this flag to keep the BLAST temporary files.

Examples:
    ${0##*/} -in input.fasta -genome genome.fasta -out out_85.txt

    ${0##*/} -in input.fasta -genome genome.fasta -out out.txt -sim 90 -keepblast
    mv out.txt out_90.txt

    ${0##*/} -in input.fasta -genome genome.fasta -out out.txt -sim 95 -afterblast 
    mv out.txt out_95.txt

EOF
}

# ----------------------------------------------------------------------------
#            PARAMETERS
# ----------------------------------------------------------------------------

#!/bin/bash

after_blast_flag=0
keep_blast_flag=0

# Read arguments
while [ "$1" != "" ]; do
    case $1 in
      -h | --help )  show_help
                       exit
                       ;;
        -in )    shift
                 fasta_file=$1
                 ;;
        -out )   shift
                 output_file=$1
                 ;;
        -sim )   shift
                 sim_threshold=$1
                 ;;
        -genome ) shift
                  genome_file=$1
                  ;;
        -afterblast ) 
                  after_blast_flag=1  # Устанавливаем флаг в 1, если флаг -afterblast присутствует
                  ;;
        -keepblast ) 
                  keep_blast_flag=1  # Устанавливаем флаг в 1, если флаг -afterblast присутствует
                  ;;
        * )      echo "Invalid parameter: $1"
                 exit 1
    esac
    shift
done

# Check if FASTA file parameter is provided
if [ -z "$fasta_file" ]; then
    echo "FASTA file not specified"
    exit 1
fi

# Check if output file parameter is provided
if [ -z "$output_file" ]; then
    echo "Output file not specified"
    exit 1
fi

# Check if genome file parameter is provided
if [ -z "$genome_file" ]; then
    echo "Genome file not specified"
    exit 1
fi

# Check if the FASTA file exists
if [ ! -f "$fasta_file" ]; then
    echo "Input FASTA file not found: $fasta_file"
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


# Check if BLAST database exists
if [ ! -f "${genome_file}.nhr" ]; then
    echo "BLAST database for $fasta_file not found. Creating database..."
    makeblastdb -in "$genome_file" -dbtype nucl > /dev/null
fi


blast_res="${output_file}.blast.tmp"

if [ "$after_blast_flag" -eq 1 ]; then
    if [ ! -f "${blast_res}" ]; then
        echo "Blast results file not found: ${fasta_file}"
        exit 1
    fi
else
    pokaz_stage "BLAST search..."
    blastn -db ${genome_file} -query ${fasta_file} -out ${blast_res} -outfmt "6 qseqid qstart qend sstart send pident length sseqid" -perc_identity ${sim_threshold}
fi


pokaz_stage "Similarity search..."
Rscript sim/sim_search.R --in_file ${fasta_file} --res ${blast_res} --out ${output_file} --sim ${sim_threshold}

# Remove BLAST temporary file if not needed
if [ "$keep_blast_flag" -ne 1 ] && [ "$after_blast_flag" -ne 1 ]; then
    rm ${blast_res}
fi

pokaz_message "Done!"