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
for query_file_path in ${path_gaps}*query*.fasta; do

    query_file=$(basename "$query_file_path")

    base_file="${query_file/query/base}"

    out_file="${query_file/query/out}"
    out_file="${out_file%.fasta}.txt"

    # Create BLAST database
    if [ ! -f "${base_file}.nhr" ] && [ ! -f "${base_file}.nin" ] && [ ! -f "${base_file}.nsq" ]; then
        makeblastdb -in ${path_gaps}${base_file} -dbtype nucl -out ${base_file}  &> /dev/null
    fi
    
    # Execute BLAST search
    if [[ ! -e ${path_gaps}${out_file} ]]; then
        blastn -db ${base_file} \
               -query ${path_gaps}${query_file}  \
               -out ${path_gaps}${out_file} \
               -outfmt "7 qseqid qstart qend sstart send pident length qseq sseq sseqid" \
               -max_hsps 10 &
    fi


    # BLAST seqarch in "cross" mode
    if [[ $query_file != *residual* ]]; then
        base_file="${query_file/query/residual_base}"
        out_file="${query_file/query/out_on_residual}"
        out_file="${out_file%.fasta}.txt"
    else
        base_file="${query_file/residual_query/residual}"
        out_file="${query_file/query/out_on_core}"
        out_file="${out_file%.fasta}.txt"
    fi

    if [[ ! -e ${path_gaps}${out_file} ]]; then
        blastn -db ${base_file} \
               -query ${path_gaps}${query_file}  \
               -out ${path_gaps}${out_file} \
               -outfmt "7 qseqid qstart qend sstart send pident length qseq sseq sseqid" \
               -max_hsps 5 &
    fi

    # Process tracking for parallel tasks
    pids="$pids $!"
    blast_number=$(pgrep -c blastn)

    if (( ${blast_number} > $cores )); then
        echo ${blast_number}
        wait -n
    fi

done


wait $pids

echo "  Done!"