#!/bin/bash
INSTALLED_PATH=$(Rscript -e "cat(system.file(package = 'pannagram'))")
# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

# source $INSTALLED_PATH/utils/chunk_error_control.sh

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source $INSTALLED_PATH/utils/utils_bash.sh

# ----------------------------------------------------------------------------
#             PARAMETERS
# ----------------------------------------------------------------------------

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -cores) cores="$2"; shift 2;;
        -path_mafft_in) path_mafft_in="$2"; shift 2;;
        -path_mafft_out) path_mafft_out="$2"; shift 2;;
        -log_path) log_path=$2; shift 2;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
done



if [ -z "$cores" ]; then
    cores=1
fi
# echo "Number of cores ${cores}"

check_missing_variable "path_mafft_in"
check_missing_variable "path_mafft_out"

path_mafft_in=$(add_symbol_if_missing "$path_mafft_in" "/")
path_mafft_out=$(add_symbol_if_missing "$path_mafft_out" "/")



if [ ! -d "${path_mafft_out}" ]; then
    mkdir -p "${path_mafft_out}"
fi


# ----------------------------------------------------------------------------
#                MAIN
# ----------------------------------------------------------------------------

trap "echo 'Script interrupted'; exit" INT

# Defining the function to be executed in parallel
mafft_task() {
    input_file=$1
    log_path=$2

    base_name=$(basename "${input_file}" .fasta)
    output_file="${path_mafft_out}/${base_name}_aligned.fasta"

    # if [ -e "$output_file" ]; then
    #     return  # Skip if output file exists
    # fi

    # Create a log file
    if [ -d "$log_path" ]; then
        file_log="${log_path}${base_name}.log"
        > "$file_log"
    else
        file_log="/dev/null"
    fi

    # Check if timeout exists
    if command -v timeout >/dev/null 2>&1; then
        TIMEOUT_CMD="timeout --foreground 60"
    elif command -v gtimeout >/dev/null 2>&1; then
        TIMEOUT_CMD="gtimeout --foreground 60"
    else
        TIMEOUT_CMD=""  # нет timeout — запускаем напрямую
    fi

    # Run MAFFT
    {
        $TIMEOUT_CMD mafft --quiet --op 3 --ep 0.1 --treeout "$input_file" > "$output_file" 2>/dev/null
    } || {
        rm -f "$output_file"
        echo "MAFFT command failed for $input_file, but continuing with other files."
    }

    if [ -d "$log_path" ]; then
        echo "Done." >> "$file_log"
    fi
}

export -f mafft_task
export path_mafft_out

# Run tasks in parallel, using a number of parallel jobs equal to the number of cores
# find "${path_mafft_in}" -maxdepth 1 -name "*.fasta" | parallel --will-cite -j $cores mafft_task {}  ${log_path}

find "${path_mafft_in}" -maxdepth 1 -name "*.fasta" ! -name "*_extra.fasta" | parallel --will-cite -j $cores mafft_task {} ${log_path}

