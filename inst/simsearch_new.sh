#!/bin/bash

INSTALLED_PATH=$(Rscript -e "cat(system.file(package = 'pannagram'))")

echo $INSTALLED_PATH

source "$INSTALLED_PATH/utils/chunk_error_control.sh"
source "$INSTALLED_PATH/utils/utils_bash.sh"
source "$INSTALLED_PATH/utils/help_simsearch.sh"
source "$INSTALLED_PATH/utils/argparse_simsearch.sh"

# ----------------------------------------------------------------------------
#            MAIN
# ----------------------------------------------------------------------------

# Fix the output file result
output_pref=$(add_symbol_if_missing "$output_pref" "/")
if [ ! -d "${output_pref}" ]; then
    mkdir -p "${output_pref}"
fi
output_pref="${output_pref}simsearch"
pokaz_message "Prefex for the output file was changed to ${output_pref}"

# ---------------------------------------------
# Files for the blast

# Add all FASTA files from path_genome to db_files if path_genome is not empty
if [ ! -z "$path_genome" ]; then

    path_genome=$(add_symbol_if_missing "$path_genome" "/")
    db_files=()

    fasta_extensions=("fa" "fasta" "fna" "fas" "ffn" "frn")

    for ext in "${fasta_extensions[@]}"; do
        for genome_file in "$path_genome"/*.$ext; do
            if [ -e "$genome_file" ]; then
                db_file=$(basename "$genome_file")
                db_files+=("$db_file")
            fi
        done
    done

fi

# Add file_seq to db_files if it's not empty
if [ ! -z "$file_seq" ]; then
    path_genome="$(dirname "$file_seq")/"
    db_files=("$(basename "$file_seq")")
fi

# Add file_genome to db_files if it's not empty
if [ ! -z "$file_genome" ]; then
    path_genome="$(dirname "$file_genome")/"
    db_files=("$(basename "$file_genome")")
fi

# ---------------------------------------------
# Run the pipeline

for db_file in "${db_files[@]}"; do

    db_pref=$(basename "$db_file")
    db_pref=${db_pref%%.*}
    file_out_cnt="${output_pref}.${db_pref}_${sim_threshold}_${coverage}.cnt"
    echo "File with counts ${file_out_cnt}"
    if [ -f "$file_out_cnt" ]; then
       echo "Counts for ${db_name} estimated."
       continue
    fi

    # ---------------------------------------------
    # Check if the BLAST database exists for the current file
    db_file_full="${path_genome}${db_file}"
    echo "Database ${db_file_full}"
    if [ ! -f "${db_file_full}.nhr" ]; then
        pokaz_stage "Creating database for $db_file..."
        makeblastdb -in "$db_file_full" -dbtype nucl > /dev/null
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
        blastn \
            -db "$db_file_full" \
            -query "$file_input" \
            -out "$blast_res" \
            -outfmt "6 qseqid qstart qend sstart send pident length sseqid qlen slen" \
            -perc_identity "$((sim_threshold - 1))"
    fi

    # Check if the BLAST results file is empty
    if [ ! -s "${blast_res}" ]; then
        pokaz_message "Blast result is empty for ${db_file}"
        continue
    fi

    # ---------------------------------------------
    # Proceed to similarity search
    pokaz_stage "Search in ${db_file}: similarity ${sim_threshold}, coverage ${coverage}..."

    # Determine if the search is on a set of sequences or a genome
    if [ -n "$file_seq" ]; then
        # On a set of sequences
        Rscript "$INSTALLED_PATH/sim/sim_in_seqs.R" \
            --in_file "$file_input" \
            --res "$blast_res" \
            --out "${output_pref}.${db_name}.rds" \
            --sim "$sim_threshold" \
            --use_strand "$use_strand" \
            --db_file "$db_file_full" \
            --coverage "${coverage}"
    else
        # On a genome
        Rscript "$INSTALLED_PATH/sim/sim_in_genome.R" \
            --in_file "$file_input" \
            --res "$blast_res" \
            --out "${output_pref}.${db_name}" \
            --sim "$sim_threshold" \
            --coverage "$coverage"
    fi

    # Remove the BLAST temporary file if not needed
    if [ "$keep_blast_flag" -ne 1 ] && [ "$after_blast_flag" -ne 1 ]; then
        rm "$blast_res"
    fi

done

# Combine all files to the total count file
if [[ -z "$file_seq" && -z "$file_genome" ]]; then
    Rscript "$INSTALLED_PATH/sim/sim_in_genome_combine.R" \
        --out "$output_pref" \
        --sim "$sim_threshold" \
        --coverage "$coverage"
fi

pokaz_message "Done!"