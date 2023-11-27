# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

# Exit immediately if any command returns a non-zero status
set -e

# Keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# Define a trap for the EXIT signal
trap 'catch $?' EXIT

# Function to handle the exit signal
catch() {
    # Check if the exit code is non-zero
    if [ $1 -ne 0 ]; then
        echo "\"${last_command}\" command failed with exit code $1."
    fi
}

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source utils_bash.sh

# ----------------------------------------------------------------------------
#             PARAMETERS
# ----------------------------------------------------------------------------

# Инициализация переменных


# Разбор аргументов командной строки
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --cores) cores="$2"; shift ;;
        --path.mafft.in) path_mafft_in="$2"; shift ;;
        --path.mafft.out) path_mafft_out="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done


# Установка значения по умолчанию для cores, если оно не было задано
if [ -z "$cores" ]; then
    cores=1
fi

check_missing_variable "path_mafft_in"
check_missing_variable "path_mafft_out"

path_mafft_in=$(add_symbol_if_missing "$path_mafft_in" "/")
path_mafft_out=$(add_symbol_if_missing "$path_mafft_out" "/")


# Проверяем, существует ли выходная директория, если нет - создаем
if [ ! -d "${path_mafft_out}" ]; then
    mkdir -p "${path_mafft_out}"
fi


# ----------------------------------------------------------------------------
#                MAIN
# ----------------------------------------------------------------------------

trap "echo 'Script interrupted'; exit" INT

# Перебираем все .fasta файлы в входной директории

# for input_file in "${path_mafft_in}"/*.fasta; do
find "${path_mafft_in}" -maxdepth 1 -name "*.fasta" | head -n 100 | while read input_file; do
    # Проверяем, есть ли файлы .fasta
    if [ ! -e "$input_file" ]; then
        echo "No .fasta files found in input directory for MAFFT."
        break
    fi

    # Извлекаем имя файла без расширения для использования в выходном файле
    base_name=$(basename "${input_file}" .fasta)
    output_file="${path_mafft_out}/${base_name}_aligned.fasta"


    if [ -e "$output_file" ]; then
        continue
    fi
    
    echo ${output_file}

    # Запускаем mafft с ограничением времени выполнения
    # timeout 10 mafft --op 5 --quiet --maxiterate 100 "${input_file}" > "${output_file}"
    # exit_status=$?

    # if [ $exit_status -eq 124 ]; then
    #     echo "Command 'mafft' on file ${input_file} took too long and was terminated."
    # elif [ $exit_status -ne 0 ]; then
    #     echo "Command 'mafft' failed on file ${input_file} with exit status $exit_status."
    # fi
done
