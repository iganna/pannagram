# Check the number of chromosomes in the accessions genome
if [[ ${nchr} -ne 0 ]]; then
    for acc_name in "${acc_set[@]}"; do
        for ext in "${genome_extensions[@]}"; do
            # Find the file with the current acc and extension
            file="${path_in}/${acc_name}.${ext}"
            
            # Check if the file exists
            if [[ -f "$file" ]]; then
                count=$(grep -c '^>' "$file")
                if [[ ${count} -lt ${nchr} ]]; then
                    pokaz_error "Error: Number of chromosomes in the reference genome ${acc_name} is ${count}. Must be at least ${nchr}."
                    exit 1
                fi
            fi
        done
    done
fi