#!/bin/bash

# Выделение текста цветом
echo -e "\e[38;2;52;252;252m* BLAST of gaps between syn blocks\e[0m"

# Функция показа инструкции по использованию
print_usage() {
  echo "-path_gaps"
}

# Парсинг аргументов командной строки
while [ $# -gt 0 ]; do
    case $1 in
        -path_gaps) 
            path_gaps=$2
            shift
            ;;
        *) 
            print_usage
            echo "$0: error - unrecognized option $1" 1>&2
            exit 1
            ;;
    esac
    shift
done


cores=30
pids=""

# echo  ${path_gaps}

# Path to databases
path_db=${path_gaps}/db
if [ ! -d "${path_db}" ]; then
  mkdir -p "${path_db}"
fi


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


# ======================================================================================================
# Create blast databases
# Declare an array to store the process IDs (PIDs) of the background processes.
declare -a pids=()

# Loop through all files with a pattern '*query*.fasta' in the directory specified by 'path_gaps'.
for query_file_path in ${path_gaps}*query*.fasta; do
    # Extract the file name from the full file path.
    query_file=$(basename "$query_file_path")
    # Replace 'query' with 'base' in the file name.
    base_file="${query_file/query/base}"

    # Check if BLAST database files do not exist.
    if [ ! -f "${base_file}.nhr" ] && [ ! -f "${base_file}.nin" ] && [ ! -f "${base_file}.nsq" ]; then
        # Create BLAST database and run the process in the background.
        makeblastdb -in ${path_gaps}${base_file} -dbtype nucl -out ${path_db}${base_file}  &> /dev/null &
        # Add the PID of the last background process to the array.
        pids+=($!)
    fi

    # # If the number of background processes equals or exceeds the number specified by '$cores'.
    # if (( ${#pids[@]} >= $cores )); then
    #     # Wait for any one of the background processes to complete.
    #     wait -n
    #     # Check and remove completed processes from the array.
    #     for pid in "${pids[@]}"; do
    #         # If the process is no longer running, remove its PID from the array.
    #         if ! kill -0 $pid 2>/dev/null; then
    #             pids=("${pids[@]/$pid}")
    #         fi
    #     done
    # fi
done

# Wait for all remaining background processes to complete.
for pid in "${pids[@]}"; do
    wait $pid
done


# ======================================================================================================
#  BLAST
# Declare an array to store the process IDs (PIDs) of the background processes.
declare -a pids=()

# Loop through all files matching the pattern '*query*.fasta' in the directory specified by 'path_gaps'.
for query_file_path in ${path_gaps}*query*.fasta; do
    query_file=$(basename "$query_file_path")
    base_file="${query_file/query/base}"

    out_file="${query_file/query/out}"
    out_file="${out_file%.fasta}.txt"
    
    # Execute BLAST search
    if [[ ! -e ${path_gaps}${out_file} ]]; then
        blastn -db ${path_db}${base_file} \
               -query ${path_gaps}${query_file}  \
               -out ${path_gaps}${out_file} \
               -outfmt "7 qseqid qstart qend sstart send pident length qseq sseq sseqid" \
               -max_hsps 10 > /dev/null 2>> log_err.txt &
        # Store the PID of the last background process.
        pids+=($!)
    fi
    # If the number of background processes equals or exceeds the number specified by '$cores'.
    if (( ${#pids[@]} >= $cores )); then
        # Wait for any one of the background processes to complete.
        wait -n
        # Check and remove completed processes from the array.
        for pid in "${pids[@]}"; do
            # If the process is no longer running, remove its PID from the array.
            if ! kill -0 $pid 2>/dev/null; then
                pids=("${pids[@]/$pid}")
            fi
        done
    fi


    # ==============================================================
    # BLAST search in "cross" mode
    if [[ $query_file != *residual* ]]; then
        base_file="${query_file/query/residual_base}"
        out_file="${query_file/query/out_on_residual}"
        out_file="${out_file%.fasta}.txt"
    else
        base_file="${query_file/residual_query/base}"
        out_file="${query_file/query/out_on_core}"
        out_file="${out_file%.fasta}.txt"
    fi

    if [[ ! -e ${path_gaps}${out_file} ]]; then
        blastn -db ${path_db}${base_file} \
               -query ${path_gaps}${query_file}  \
               -out ${path_gaps}${out_file} \
               -outfmt "7 qseqid qstart qend sstart send pident length qseq sseq sseqid" \
               -max_hsps 5 > /dev/null 2>> log_err.txt &
        # Store the PID of the last background process.
        pids+=($!)
    fi

    # If the number of background processes equals or exceeds the number specified by '$cores'.
    if (( ${#pids[@]} >= $cores )); then
        # Wait for any one of the background processes to complete.
        wait -n
        # Check and remove completed processes from the array.
        for pid in "${pids[@]}"; do
            # If the process is no longer running, remove its PID from the array.
            if ! kill -0 $pid 2>/dev/null; then
                pids=("${pids[@]/$pid}")
            fi
        done
    fi
done

# Wait for all remaining background processes to complete.
for pid in "${pids[@]}"; do
    wait $pid
done


echo "  Done!"