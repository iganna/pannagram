#!/bin/bash
INSTALLED_PATH=$(Rscript -e "cat(system.file(package = 'pannagram'))")
# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

source $INSTALLED_PATH/utils/chunk_error_control.sh

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source $INSTALLED_PATH/utils/utils_bash.sh

print_usage() {
  echo "-path_gaps"
  echo "-log_path"
  echo "-p_ident"
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
        -p_ident) p_ident=$2; shift 2;;
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
p_ident="${p_ident:-85}"


# Path to databases
path_db=${path_gaps}db/
if [ ! -d "${path_db}" ]; then
  mkdir -p "${path_db}"
fi

# ----------------------------------------------------------------------------
#                 FUNC: CREATE DATABASE
# ----------------------------------------------------------------------------

export path_gaps
export path_db
export log_path
export p_ident

function process_db {
    query_file_path="$1"

    # Extract the file name from the full file path.
    query_file=$(basename "$query_file_path")
    # Replace 'query' with 'base' in the file name.
    base_file="${query_file/query/base}"


    # If log file exists and has the word "Done" - then don't run the blast again
    file_log="${log_path}${query_file}_db.log"
    if [ -f "$file_log" ]; then
        if grep -q "Done" "$file_log"; then
            echo "Over." >> "$file_log"
            return 0
        fi
    fi
    echo "New attempt:" > "$file_log"  # Create or empty the log file

    # Check if BLAST database files do not exist.
    if [ -f "${path_gaps}${base_file}" ] && [ ! -f "${base_file}.nhr" ] && [ ! -f "${base_file}.nin" ] && [ ! -f "${base_file}.nsq" ]; then
        # Create BLAST database
        makeblastdb -in ${path_gaps}${base_file} -dbtype nucl -out ${path_db}${base_file}  &> /dev/null
    fi

    if [ -d "$log_path" ]; then
        echo "Done." >> "$file_log"
    fi
}


# ----------------------------------------------------------------------------
#                 FUNC: RUN BLAST
# ----------------------------------------------------------------------------

function process_blast_normal {
    query_file_path="$1"

    query_file=$(basename "$query_file_path")
    base_file="${query_file/query/base}"

    out_file="${query_file/query/out}"
    out_file="${out_file%.fasta}.txt"

    # ------------------------------------------
    # If log file exists and has the word "Done" - then don't run the blast again
    file_log="${log_path}${query_file}_normal.log"
    if [ -f "$file_log" ]; then
        if grep -q "Done" "$file_log"; then
            echo "Over." >> "$file_log"
            return 0
        fi
    fi
    echo "New attempt:" > "$file_log"  # Create or empty the log file

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
               -perc_identity "${p_ident}" \
               -max_hsps 10  >> "$file_log" 2>&1
    fi

    if [ -d "$log_path" ]; then
        echo "Done." >> "$file_log"
    fi
}


# BLAST search in "cross" mode
function process_blast_cross {
    query_file_path="$1"

    query_file=$(basename "$query_file_path")
    base_file="${query_file/query/base}"

    out_file="${query_file/query/out}"
    out_file="${out_file%.fasta}.txt"

    # ------------------------------------------
    # If log file exists and has the word "Done" - then don't run the blast again
    file_log="${log_path}${query_file}_cross.log"
    if [ -f "$file_log" ]; then
        if grep -q "Done" "$file_log"; then
            echo "Over." >> "$file_log"
            return 0
        fi
    fi
    echo "New attempt:" > "$file_log"  # Create or empty the log file

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
               -perc_identity "${p_ident}" \
               -max_hsps 5 >> "$file_log" 2>&1
    fi

    if [ -d "$log_path" ]; then
        echo "Done." >> "$file_log"
    fi
}

export -f process_blast_normal
export -f process_blast_cross
export -f process_db

# ----------------------------------------------------------------------------
#                 MAIN
# ----------------------------------------------------------------------------

files_acc=($(find ${path_gaps} -name '*query*.fasta'))

echo "${files_acc[@]}"

parallel -j ${cores} process_db ::: ${files_acc[@]}
parallel -j ${cores}  process_blast_normal ::: "${files_acc[@]}" 
parallel -j ${cores}  process_blast_cross ::: "${files_acc[@]}" 


# find ${path_gaps} -name '*query*.fasta' | parallel -j ${cores} process_db {} $path_gaps $path_db ${log_path}
# find ${path_gaps} -name '*query*.fasta' | parallel -j ${cores} process_blast_normal {} $path_gaps $path_db ${log_path} ${p_ident}
# find ${path_gaps} -name '*query*.fasta' | parallel -j ${cores} process_blast_cross {} $path_gaps $path_db ${log_path} ${p_ident}


# pokaz_message "Done!"
