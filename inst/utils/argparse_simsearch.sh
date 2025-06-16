# ----------------------------------------------------------------------------
#            PARAMETERS
# ----------------------------------------------------------------------------
if [ $# -eq 0 ]; then
    pokaz_error "No arguments provided!"
    help_in_box
    exit 0
fi

after_blast_flag=0
keep_blast_flag=0
use_strand=T

# Read arguments
while [ "$1" != "" ]; do
    case $1 in
        -h | --help ) show_help; exit ;;
        -in_seq )    file_input=$2; shift 2 ;;
        -out )       output_pref=$2; shift 2 ;;
        -sim )       sim_threshold=$2; shift 2 ;;
        -cov )       coverage=$2; shift 2 ;;

        -on_seq )    file_seq=$2; shift 2 ;;
        -on_genome ) file_genome=$2; shift 2 ;;
        -on_path )   path_genome=$2; shift 2 ;;

        -afterblast ) after_blast_flag=1; shift ;;
        -keepblast )  keep_blast_flag=1; shift ;;

        -strandfree ) use_strand=F; shift ;;
        * ) pokaz_error "Unknown parameter: $1"; help_in_box; exit 1;;
    esac
done

# Ensure only one of -on_seq, -on_genome, -on_path is set
count=0
[ ! -z "$file_seq" ] && count=$((count + 1))
[ ! -z "$file_genome" ] && count=$((count + 1))
[ ! -z "$path_genome" ] && count=$((count + 1))

if [ $count -ne 1 ]; then
    pokaz_error "Error: You must specify exactly one of -on_seq, -on_genome, or -on_path."
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

# Check if similarity threshold parameter is provided
if [ -z "$sim_threshold" ]; then
    sim_threshold=85
    pokaz_message "Similarity threshold not specified, default: ${sim_threshold}"
fi

# Check if coverage parameter is provided. If not - set qeual to sim
if [ -z "$coverage" ]; then
    coverage=${sim_threshold}
    pokaz_message "Coverage not specified, default: ${sim_threshold}"
fi
