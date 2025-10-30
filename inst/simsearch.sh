#!/bin/bash

INSTALLED_PATH=$(Rscript -e "cat(system.file(package = 'pannagram'))")

FASTA_SUFFIX=("fa" "fasta" "fas" "fna" "fn" "ffn" "faa")

source "$INSTALLED_PATH/utils/chunk_error_control.sh"
source "$INSTALLED_PATH/utils/utils_bash.sh"
source "$INSTALLED_PATH/utils/help_simsearch.sh"
source "$INSTALLED_PATH/utils/argparse_simsearch.sh"

# ----------------------------------------------------------------------------
#            MAIN
# ----------------------------------------------------------------------------

output_pref=$(add_symbol_if_missing "$output_pref" "/")
# pokaz_message "Prefex for the output file was changed to ${output_pref}"

# Create intermediate directory
intermediate_dir="${output_pref}.intermediate"
mkdir -p "$intermediate_dir"
# pokaz_message "Intermediate files will be stored in $intermediate_dir"


# Add all FASTA files from path_genome to db_files if path_genome is not empty
if [ ! -z "$path_genome" ]; then

    path_genome=$(add_symbol_if_missing "$path_genome" "/")
    db_files=()

    for ext in "${FASTA_SUFFIX[@]}"; do
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
    file_out_cnt="${output_pref}${db_pref}_${sim_threshold}_${coverage}.cnt"
    if [ -f "$file_out_cnt" ]; then
       pokaz_message "Counts for ${db_pref} have been estimated."
       continue
    fi

    # ---------------------------------------------
    # Check if the BLAST database exists for the current file
    db_file_full="${path_genome}${db_file}"
    db_name="${db_pref}"  # Use consistent name
    
    # Create DB in intermediate_dir
    db_path_intermediate="${intermediate_dir}/${db_name}"
    if [ ! -f "${db_path_intermediate}.phr" ] && [ ! -f "${db_path_intermediate}.nhr" ]; then
        pokaz_stage "Creating database for $db_file..."
        makeblastdb -in "$db_file_full" -dbtype "$dbtype" -out "$db_path_intermediate" > /dev/null
    fi

    blast_res="${intermediate_dir}/${db_name}.blast.tmp"

    # Check if BLAST results should be used from an existing file
    if [ "$after_blast_flag" -eq 1 ]; then
        if [ ! -f "${blast_res}" ]; then
            pokaz_error "Blast results file not found: ${file_input}"
            exit 1
        fi
    else
        # Perform BLAST search
        pokaz_stage "BLAST search $(basename "$file_input") in $(basename "$db_file")..."
        if [ "$use_aa" -eq 1 ]; then
            blast_res_pre="${blast_res}.pre"

            $blast_cmd \
                -db "$db_path_intermediate" \
                -query "$file_input" \
                -out "$blast_res_pre" \
                -outfmt "6 qseqid qstart qend sstart send pident length sseqid qlen slen" \
                -num_threads "$cores"
            
            # Remove previous file with aa reults
            rm -f "$blast_res"

            Rscript "$INSTALLED_PATH/sim/sim_modify_blast_results.R" \
                --file.init "$blast_res_pre" \
                --file.mod "$blast_res"
        else
            $blast_cmd \
                -db "$db_path_intermediate" \
                -query "$file_input" \
                -out "$blast_res" \
                -outfmt "6 qseqid qstart qend sstart send pident length sseqid qlen slen" \
                -perc_identity "$((sim_threshold - 1))" \
                -num_threads "$cores"
        fi
    fi

    # Check if the BLAST results file is empty
    if [ ! -s "${blast_res}" ]; then
        pokaz_message "Blast result is empty for ${db_file}"
        continue
    fi

    if [ "$stop_after_blast_flag" -eq 1 ]; then
        pokaz_message "Simsearch was stopped after BLAST, but before the output file was created."
        exit 0
    fi

    # ---------------------------------------------
    # Proceed to similarity search
    pokaz_stage "Processing BLAST-results..."

    # Determine if the search is on a set of sequences or a genome
    if [ -n "$file_seq" ]; then
        # On a set of sequences
        Rscript "$INSTALLED_PATH/sim/sim_in_seqs.R" \
            --in_file "$file_input" \
            --res "$blast_res" \
            --out "${output_pref}${db_name}.rds" \
            --sim "$sim_threshold" \
            --use_strand "$use_strand" \
            --db_file "$db_file_full" \
            --coverage "${coverage}"
    else
        # On a genome
        Rscript "$INSTALLED_PATH/sim/sim_in_genome.R" \
            --in_file "$file_input" \
            --res "$blast_res" \
            --out "${output_pref}${db_name}" \
            --sim "$sim_threshold" \
            --coverage "$coverage"

        # pokaz_message "File with counts was generated: ${file_out_cnt}"
    fi

    # Remove the BLAST temporary file if not needed
    if [ "$keep_blast_flag" -ne 1 ] && [ "$after_blast_flag" -ne 1 ]; then
        rm -f "$blast_res"
        rm -f "$blast_res_pre"
    fi

done

# Combine all files to the total count file
if [[ -z "$file_seq" && -z "$file_genome" ]]; then
    Rscript "$INSTALLED_PATH/sim/sim_in_genome_combine.R" \
        --out_dir "$output_pref" \
        --sim "$sim_threshold" \
        --coverage "$coverage"
fi

pokaz_message "Done!"