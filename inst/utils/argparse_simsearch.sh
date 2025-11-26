# ----------------------------------------------------------------------------
#            PARAMETERS
# ----------------------------------------------------------------------------

# Function to check if a file has a valid FASTA extension
has_extension() {
    local filename="$1"
    shift
    local extensions=("$@")
    local ext="${filename##*.}"

    for allowed_ext in "${extensions[@]}"; do
        if [[ "$ext" == "$allowed_ext" ]]; then
            return 0
        fi
    done
    return 1
}

if [ $# -eq 0 ]; then
    pokaz_error "No arguments provided!"
    help_in_box
    exit 0
fi

if command -v nproc &>/dev/null; then
    max_cores=$(nproc)
else
    max_cores=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)
fi

after_blast_flag=0
stop_after_blast_flag=0
keep_blast_flag=0
use_strand=T
use_aa=0
cores=1

# Read arguments
cmdline="$(basename "$0") $@"
while [ "$1" != "" ]; do
    case $1 in
        -h | --help )         show_help; exit ;;
        -in_seq | -in-seq | -query_seqs )   file_input=$2;    shift 2 ;;
        -out )                             output_pref=$2;   shift 2 ;;
        -sim | -similarity )               similarity=$2;    shift 2 ;;
        -cov | -coverage )                 coverage=$2;      shift 2 ;;

        -on_seq    | -on-seq    | -target_seqs )    file_seq=$2;    shift 2 ;;
        -on_genome | -on-genome | -target_genome ) file_genome=$2; shift 2 ;;
        -on_path   | -on-path   | -target_path )   path_genome=$2; shift 2 ;;

        -afterblast ) after_blast_flag=1;         shift ;;
        -keepblast )  keep_blast_flag=1;          shift ;;
        -stopafterblast) stop_after_blast_flag=1; shift ;;
        -aa|-prot )   use_aa=1;                   shift ;;

        -strandfree ) use_strand=F;               shift ;;

        -cores)
        if [[ -n "${2-}" ]] && [[ "$2" =~ ^[0-9]+$ ]]; then
            cores="$2"
            shift 2
        else
            pokaz_error "Error: -cores requires a positive integer."
            exit 1
        fi
        ;;

        * ) pokaz_error "Unknown parameter: $1"; 
            pokaz_attention "!!! Please check your command: ${cmdline}";
            help_in_box; exit 1;;
    esac
done

if (( cores < 1 )); then
  cores=1
elif (( cores > max_cores )); then
  cores=$max_cores
fi

# Ensure only one of -on_seq, -on_genome, -on_path is set
count=0
[ ! -z "$file_seq" ] && count=$((count + 1))
[ ! -z "$file_genome" ] && count=$((count + 1))
[ ! -z "$path_genome" ] && count=$((count + 1))

if [ $count -ne 1 ]; then
    pokaz_error "Error: One or mode of required paraments were not provided."
    help_in_box
    exit 1
fi

# Check if FASTA file parameter is provided
if [ -z "$file_input" ]; then
    pokaz_error "Error: FASTA file not specified"
    help_in_box
    exit 1
fi

# Check if the FASTA file exists
if [ ! -f "$file_input" ]; then
    pokaz_error "Error: Input FASTA file not found: $file_input"
    exit 1
fi

# Check if output file parameter is provided
if [ -z "$output_pref" ]; then
    pokaz_error "Error: Output file not specified"
    help_in_box
    exit 1
fi

# Validate -in_seq
if [ -n "$file_input" ]; then
    if ! has_extension "$file_input" "${FASTA_SUFFIX[@]}"; then
        pokaz_error "Error: Invalid FASTA file: '$file_input'. Acceptable suffixes are: ${FASTA_SUFFIX[*]}"
        exit 1
    fi
fi

# Validate -on_seq if set
if [ -n "$file_seq" ]; then
    if [ ! -f "$file_seq" ]; then
        pokaz_error "Error: File not found: $file_seq"
        exit 1
    fi
    if ! has_extension "$file_seq" "${FASTA_SUFFIX[@]}"; then
        pokaz_error "Error: Invalid FASTA file: '$file_seq'. Acceptable suffixes are: ${FASTA_SUFFIX[*]}"
        exit 1
    fi
fi

# Validate -on_genome if set
if [ -n "$file_genome" ]; then
    if [ ! -f "$file_genome" ]; then
        pokaz_error "Error: File not found: $file_genome"
        exit 1
    fi
    if ! has_extension "$file_genome" "${FASTA_SUFFIX[@]}"; then
        pokaz_error "Error: Invalid FASTA file: '$file_genome'. Acceptable suffixes are: ${FASTA_SUFFIX[*]}"
        exit 1
    fi
fi

# Validate -on_path if set
if [ -n "$path_genome" ]; then
    if [ ! -d "$path_genome" ]; then
        pokaz_error "Error: Path genome directory does not exist: $path_genome"
        exit 1
    fi

    found_fasta=0
    for genome_file in "$path_genome"/*; do
        if [ -f "$genome_file" ] && has_extension "$genome_file" "${FASTA_SUFFIX[@]}"; then
            found_fasta=1
            break
        fi
    done

    if [ "$found_fasta" -eq 0 ]; then
        pokaz_error "Error: No FASTA files found in directory: $path_genome"
        exit 1
    fi
fi

# Check if similarity threshold parameter is provided
if [ -z "$similarity" ]; then
    if [ "$use_aa" -eq 1 ]; then
        similarity=60
        pokaz_message "Similarity threshold not specified, default for '-aa': ${similarity}"
    else
        similarity=85
        pokaz_message "Similarity threshold not specified, default: ${similarity}"
    fi
else
    pokaz_message "Similarity threshold: ${similarity}"
fi

# Check if coverage parameter is provided. If not - set qeual to sim
if [ -z "$coverage" ]; then
    coverage=${similarity}
    pokaz_message "Coverage threshold not specified, default: ${coverage}"
else
    pokaz_message "Coverage threshold: ${coverage}"
fi


# Check if 'similarity' is empty
if [[ -z "$similarity" ]]; then
    echo "Error: similarity is empty"
    exit 1
fi

# Check if 'coverage' is empty
if [[ -z "$coverage" ]]; then
    echo "Error: coverage is empty"
    exit 1
fi


# Determine BLAST command and database type
if [ "$use_aa" -eq 1 ]; then
    # pokaz_message "Searching for proteins in nucleotide db"
    blast_cmd="tblastn"
    dbtype="nucl"
else
    # pokaz_message "Searching for nucleotides in nucleotide db"
    blast_cmd="blastn"
    dbtype="nucl"
fi