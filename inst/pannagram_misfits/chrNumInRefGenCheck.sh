# Check the number of chromosomes in the reference genome
if [[ ${nchr_ref} -ne 0 ]]; then
    for ref_name in "${refs_all[@]}"; do
        for ext in "${genome_extensions[@]}"; do
            # Find the file with the current acc and extension
            file="${path_ref}/${ref_name}.${ext}"
            # Check if the file exists
            if [[ -f "$file" ]]; then
                count=$(grep -c '^>' "$file")
                if [[ ${count} -lt ${nchr_ref} ]]; then
                    pokaz_error "Error: Number of chromosomes in the reference genome ${ref_name} is ${count}. Must be at least ${nchr_ref}."
                    exit 1
                fi
            fi
        done
    done
fi