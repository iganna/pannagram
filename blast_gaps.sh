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

echo  ${path_gaps}

#  BLAST
for query_file_path in ${path_gaps}*queryseq*.fasta; do

    query_file=$(basename "$query_file_path")
    
    # Extract base file prefix from query file name
    base_pref=$(echo "$query_file" | awk -F'_' '{print $1 "_" $4 "_" $5 "_" $8}')
    
    # Construct base file name
    base_file="${base_pref}_base.fasta"
    
    # Construct output file name
    out_file="${query_file%.fasta}.txt"
    
    # --- Get the exact name of base-sequence
    base_seq=$(echo "$query_file" | sed -e 's/queryseq_//' -e 's/\.fasta$//')
    # Write base sequence to a temporary file
    tmp_file="tmp.txt"
    echo "$base_seq" > ${path_gaps}${tmp_file}

    # Create BLAST database
    if [ ! -f "${base_file}.nhr" ] && [ ! -f "${base_file}.nin" ] && [ ! -f "${base_file}.nsq" ]; then
        makeblastdb -in ${path_gaps}${base_file} -dbtype nucl -parse_seqids -out ${base_file}  # &> /dev/null
    fi
    
    # Execute BLAST search
    blastn -db ${base_file} \
           -query ${path_gaps}${query_file}  \
           -out ${path_gaps}${out_file} \
           -outfmt "7 qseqid qstart qend sstart send pident length qseq sseq sseqid" \
           -seqidlist ${path_gaps}${tmp_file} \
           -max_hsps 5 &

    # Process tracking for parallel tasks
    pids="$pids $!"
    blast_number=$(pgrep -c blastn)

    if (( ${blast_number} > $cores )); then
        wait -n
    fi

done


wait $pids

echo "  Done!"