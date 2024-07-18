#!/bin/bash

# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

source utils/chunk_error_control.sh

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source utils/utils_bash.sh

print_usage() {
  echo "-path_gaps"
  echo "-log_path"
}

# ----------------------------------------------------------------------------
#                 PARAMETERS
# ----------------------------------------------------------------------------

while [ $# -gt 0 ]; do
    case $1 in
        -path_gaps) path_gaps=$2 ;shift 2;;
        -log_path) log_path=$2; shift 2;;
        -penalty) penalty=$2; shift 2;;
        -gapopen) gapopen=$2; shift 2;;
        -gapextend) gapextend=$2; shift 2;;
        -max_hsps) max_hsps=$2; shift 2;;
        -cores) cores=$2; shift 2;;
        *) 
            print_usage
            echo "$0: error - unrecognized option $1" 1>&2
            exit 1
            ;;
    esac
done

# echo  ${path_gaps}

penalty="${penalty:--2}"
gapopen="${gapopen:-10}"
gapextend="${gapextend:-2}"
max_hsps="${max_hsps:-1}"
cores="${cores:-30}"


# Path to databases
path_db=${path_gaps}db/
if [ ! -d "${path_db}" ]; then
  mkdir -p "${path_db}"
fi


# ----------------------------------------------------------------------------
#                 CREATE BLAST DATABASES
# ----------------------------------------------------------------------------

# pokaz_stage "Step 5. BLAST of gaps between syntenic matches"


function process_db {
    query_file_path="$1"
    path_gaps="$2"
    path_db="$3"
    log_path=${4}

    # Extract the file name from the full file path.
    query_file=$(basename "$query_file_path")
    # Replace 'query' with 'base' in the file name.
    base_file="${query_file/query/base}"

    # Create a log file
    if [ -d "$log_path" ]; then
        file_log="${log_path}${p_filename}_${ref_chr}_db.log"
        > "$file_log"
    else
        file_log="/dev/null"
    fi

    # Check if BLAST database files do not exist.
    if [ -f "${path_gaps}${base_file}" ] && [ ! -f "${base_file}.nhr" ] && [ ! -f "${base_file}.nin" ] && [ ! -f "${base_file}.nsq" ]; then
        # Create BLAST database
        makeblastdb -in ${path_gaps}${base_file} -dbtype nucl -out ${path_db}${base_file}  &> /dev/null
    fi

    if [ -d "$log_path" ]; then
        echo "Done." >> "$file_log"
    fi
}

export -f process_db

find ${path_gaps} -name '*query*.fasta' | parallel -j ${cores} process_db {} $path_gaps $path_db ${log_path}


# ----------------------------------------------------------------------------
#                 RUN BLAST
# ----------------------------------------------------------------------------

function process_blast {
    query_file_path="$1"
    path_gaps="$2"
    path_db="$3"
    log_path="$4"

    query_file=$(basename "$query_file_path")
    base_file="${query_file/query/base}"

    out_file="${query_file/query/out}"
    out_file="${out_file%.fasta}.txt"

    # ------------------------------------------
    # Create a log file
    if [ -d "$log_path" ]; then
        file_log="${log_path}${query_file}_normal.log"
        > "$file_log"
    else
        file_log="/dev/null"
    fi

    # Execute BLAST search
    if [[ ! -e ${path_gaps}${out_file} ]] && \
       [[ -e ${path_db}${base_file}.nhr ]] && \
       [[ -e ${path_db}${base_file}.nin ]] && \
       [[ -e ${path_db}${base_file}.nsq ]] && \
       [[ -e ${path_gaps}${query_file} ]]; then

        blastn -db ${path_db}${base_file} \
               -query ${path_gaps}${query_file}  \
               -out ${path_gaps}${out_file} \
               -outfmt "6 qseqid qstart qend sstart send pident length qseq sseq sseqid" \
               -max_hsps 10  >> "$file_log" 2>&1
    fi

    if [ -d "$log_path" ]; then
        echo "Done." >> "$file_log"
    fi

    # ------------------------------------------
    # BLAST search in "cross" mode

    # Create a log file
    if [ -d "$log_path" ]; then
        file_log="${log_path}${query_file}_cross.log"
        > "$file_log"
    else
        file_log="/dev/null"
    fi

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
               -outfmt "6 qseqid qstart qend sstart send pident length qseq sseq sseqid" \
               -max_hsps 5 >> "$file_log" 2>&1
    fi

    if [ -d "$log_path" ]; then
        echo "Done." >> "$file_log"
    fi
}

export -f process_blast

find ${path_gaps} -name '*query*.fasta' | parallel -j ${cores} process_blast {} $path_gaps $path_db ${log_path}


# pokaz_message "Done!"
