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
#            PARAMETERS
# ----------------------------------------------------------------------------

#!/bin/bash


# Read arguments
while [ "$1" != "" ]; do
    case $1 in
        -in )    shift
                 fasta_file=$1
                 ;;
        -out )   shift
                 output_file=$1
                 ;;
        -sim )   shift
                 similarity_threshold=$1
                 ;;
        -genome ) shift
                  genome_file=$1
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
if [ -z "$similarity_threshold" ]; then
    similarity_threshold=90
    echo "Similarity threshold not specified, default: ${similarity_threshold}"
fi

# Your script code goes here

# ----------------------------------------------------------------------------
#            MAIN
# ----------------------------------------------------------------------------


# Check if BLAST database exists
db_name=$(basename "$genome_file" .fasta)
if [ ! -f "${db_name}.nhr" ]; then
    echo "BLAST database for $fasta_file not found. Creating database..."
    makeblastdb -in "$genome_file" -dbtype nucl
fi

blast_res="${output_file}.blast.tmp"
blastn -db ${genome_file} -query ${fasta_file} -out ${blast_res} -outfmt "6 qseqid qstart qend sstart send pident length sseqid" 

Rscript sim_search.R -in ${fasta_file} -res ${blast_res} -out ${output_file} -sim ${similarity_threshold}

