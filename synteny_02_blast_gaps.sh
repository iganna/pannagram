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
        -path_gaps) path_gaps=$2 ;shift ;;
        -cores) cores=$2; shift ;;
        *) 
            print_usage
            echo "$0: error - unrecognized option $1" 1>&2
            exit 1
            ;;
    esac
    shift
done

# echo  ${path_gaps}

# Path to databases
path_db=${path_gaps}db/
if [ ! -d "${path_db}" ]; then
  mkdir -p "${path_db}"
fi

#echo ${path_db}


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


function process_db {
    query_file_path="$1"
    path_gaps="$2"
    path_db="$3"

    # Extract the file name from the full file path.
    query_file=$(basename "$query_file_path")
    # Replace 'query' with 'base' in the file name.
    base_file="${query_file/query/base}"

    # Check if BLAST database files do not exist.
    if [ -f "${path_gaps}${base_file}" ] && [ ! -f "${base_file}.nhr" ] && [ ! -f "${base_file}.nin" ] && [ ! -f "${base_file}.nsq" ]; then
        # Create BLAST database
        makeblastdb -in ${path_gaps}${base_file} -dbtype nucl -out ${path_db}${base_file}  &> /dev/null
    fi
}

export -f process_db

find ${path_gaps} -name '*query*.fasta' | parallel -j ${cores} process_db {} $path_gaps $path_db



# # ======================================================================================================
# #  BLAST


#!/bin/bash

function process_blast {
    query_file_path="$1"
    path_gaps="$2"
    path_db="$3"

    query_file=$(basename "$query_file_path")
    base_file="${query_file/query/base}"

    out_file="${query_file/query/out}"
    out_file="${out_file%.fasta}.txt"


    # Execute BLAST search
    if [[ ! -e ${path_gaps}${out_file} ]] && [[ -e ${path_db}${base_file}.nhr  ]] && [[ -e ${path_db}${base_file}.nin ]] && [[ -e ${path_db}${base_file}.nsq ]] && [[ -e ${path_gaps}${query_file} ]]; then
        blastn -db ${path_db}${base_file} \
               -query ${path_gaps}${query_file}  \
               -out ${path_gaps}${out_file} \
               -outfmt "7 qseqid qstart qend sstart send pident length qseq sseq sseqid" \
               -max_hsps 10 > /dev/null 2>> log_err.txt 
        echo "works!"
    fi

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

    if [[ ! -e ${path_gaps}${out_file} ]] && [[ -e ${path_db}${base_file} ]] && [[ -e ${path_gaps}${query_file} ]]; then
        blastn -db ${path_db}${base_file} \
               -query ${path_gaps}${query_file}  \
               -out ${path_gaps}${out_file} \
               -outfmt "7 qseqid qstart qend sstart send pident length qseq sseq sseqid" \
               -max_hsps 5 > /dev/null 2>> log_err.txt 
    fi
}

export -f process_blast

find ${path_gaps} -name '*query*.fasta' | parallel -j ${cores} process_blast {} $path_gaps $path_db


echo "  Done!"


