#!/bin/bash

# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------
INSTALLED_PATH=$(Rscript -e "cat(system.file(package = 'pannagram'))")

if [ -z "$INSTALLED_PATH" ]; then
    echo "Error: package 'pannagram' is not installed." >&2
    exit 1
fi

source $INSTALLED_PATH/utils/chunk_error_control.sh

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source $INSTALLED_PATH/utils/utils_bash.sh


print_usage() {
    cat << EOF
Usage: ${0##*/} -genomes_in PATH_GENOMES -genomes_out PATH_OUTPUT \\
                [-remove WORDS_REMOVE] [-remain WORDS_REMAIN] [-keep WORDS_KEEP] \\
                [-cores NUM_CORES] [-h]

   or: ${0##*/} -path_project PATH_PROJECT -genomes_out PATH_OUTPUT \\
                -rearrange [-cores NUM_CORES] [-h]

Description:
    This script processes genome files either by filtering based on specified keywords,
    or by rearranging chromosomes according to a Pannagram project alignment.

Options:
    -h, --help                    Display this help message and exit.
    -cores NUM_CORES              Number of cores for parallel processing. Default is 1.

Filtering mode (before alignment):
    -genomes_in PATH_GENOMES      Path to input genomes (FASTA format).
    -genomes_out PATH_OUTPUT      Path to output processed genomes.
    -remove WORDS_REMOVE          Comma-separated list of words to remove (e.g., "plasmid,contig").
    -remain WORDS_REMAIN          Comma-separated list of words to retain.
    -keep WORDS_KEEP              Comma-separated list of words to keep regardless of other filters.

Rearrangement mode (after alignment):
    -path_project PATH_PROJECT    Path to the Pannagram output directory
    -genomes_out PATH_OUTPUT      Path to output rearranged genomes.
    -rearrange                    Enable rearrangement of chromosomes based on the reference.

Examples:

    # Example 1: Filter genome files before alignment
    ${0##*/} -genomes_in /data/genomes -genomes_out /data/genomes_filtered \\
             -remove "contig,mitochondria" -remain "chromosome" -cores 4

    # Example 2: Rearrange genomes based on alignment results
    ${0##*/} -path_project /data/pannagram_project -genomes_out /data/genomes_rearranged \\
             -rearrange -cores 4

EOF
}



# ----------------------------------------------------------------------------
#             PARAMETERS
# ----------------------------------------------------------------------------
if [ $# -eq 0 ]; then
    pokaz_error "No arguments provided!"
    help_in_box
    exit 0
fi


# Default values
cores=1
path_project=""
genomes_in=""
genomes_out=""
words_remove=""
words_remain=""
words_keep=""
path_aln=""
rearrange=false

while [ $# -gt 0 ]; do
    case $1 in
        -h|--help)        print_usage; exit 0 ;;
        -cores)           cores=$2;         shift 2 ;;
        -path_project)    path_project=$2;  shift 2 ;;
        -genomes_in)      genomes_in=$2;    shift 2 ;;
        -genomes_out)     genomes_out=$2;   shift 2 ;;
        -remove)          words_remove=$2;  shift 2 ;;
        -remain)          words_remain=$2;  shift 2 ;;
        -keep)            words_keep=$2;    shift 2 ;;
        -rearrange)       rearrange=true;   shift ;;
        *)                pokaz_error "Unknown parameter: $1"; help_in_box; exit 1 ;;
    esac
done

# Determine mode
mode_filter=false
mode_rearrange=false

# Check: either path_project or genomes_in must be set, but not both
if [[ -n "$path_project" && -n "$genomes_in" ]]; then
    echo "Error: -path_project and -genomes_in cannot be used together" >&2
    exit 1
elif [[ -z "$path_project" && -z "$genomes_in" ]]; then
    echo "Error: either -path_project or -genomes_in must be provided" >&2
    exit 1
elif [[ -n "$path_project" ]]; then
    mode_rearrange=true

    path_project=$(add_symbol_if_missing "$path_project" "/")

    if [ "$rearrange" = false ]; then
        echo "Error: -rearrange must be set" >&2
        exit 1
    fi

    if [[ -n "$words_remove" || -n "$words_remain" || -n "$words_keep" ]]; then
        echo "Error: None of the following should be set: -words_remove, -words_remain, -words_keep" >&2
        exit 1
    fi

elif [[ -n "$genomes_in" ]]; then
    mode_filter=true

    genomes_in=$(add_symbol_if_missing "$genomes_in" "/")

    if [[ -z "$words_remove" && -z "$words_remain" && -z "$words_keep" ]]; then
        echo "Error: At least one of the following must be set: -words_remove, -words_remain, -words_keep" >&2
        exit 1
    fi

    if [ "$rearrange" = true ]; then
        echo "Error: -rearrange should not be set" >&2
        exit 1
    fi
fi


# Path with modified genomes
genomes_out=$(add_symbol_if_missing "$genomes_out" "/")
mkdir -p "$genomes_out"


pokaz_message "Number of cores: ${cores}"

# ----------------------------------------------------------------------------
#                                   MAIN
# ----------------------------------------------------------------------------

# -------------------------------------------------
if [ "$mode_filter" = true ]; then
    pokaz_stage "Filtering chromosomes by keywords."

    param_words_remove=""
    param_words_remain=""
    param_words_keep=""

    if [[ -n "${words_remove}" ]]; then
        param_words_remove="--words.remove ${words_remove}"
    fi

    if [[ -n "${words_remain}" ]]; then
        param_words_remain="--words.remain ${words_remain}"
    fi

    if [[ -n "${words_keep}" ]]; then
        param_words_keep="--words.keep ${words_keep}"
    fi

    Rscript $INSTALLED_PATH/chromotools/filter_01.R \
        --path.genomes ${genomes_in} \
        --path.filtered ${genomes_out} \
        --cores ${cores} \
        ${param_words_remove} \
        ${param_words_remain} \
        ${param_words_keep}

fi

# -------------------------------------------------
if [ "$mode_rearrange" = true ]; then

    # Path with alignments
    shopt -s nullglob
    matches=("${path_project}.intermediate/alignments/"*)
    shopt -u nullglob

    # Handle number of matches
    case ${#matches[@]} in
      1) path_aln="${matches[0]}" ;;
      0) echo "Error: No matching paths found." >&2; exit 1 ;;
      *) echo "Error: Multiple matching paths found:" >&2; printf ' - %s\n' "${matches[@]}" >&2; exit 1 ;;
    esac

    # Path with chromosomes
    path_chr="${path_project}.intermediate/chromosomes/"


    # Fix folders
    path_aln=$(add_symbol_if_missing "$path_aln" "/")
    path_chr=$(add_symbol_if_missing "$path_chr" "/")

    # Get the name of the reference genome
    ref_name=$(basename "$path_aln")
    ref_name=${ref_name#alignments_}

    pokaz_stage "Find the best rearrangements based on alignment..."
    Rscript $INSTALLED_PATH/chromotools/rearrange_01_positions.R --path.aln ${path_aln} --ref ${ref_name} \
                                                                 --path.processed ${genomes_out} --path.chr ${path_chr}

    pokaz_stage "Rearranging (splitting/merging) chromosomes..."
    Rscript $INSTALLED_PATH/chromotools/rearrange_02_genomes.R --path.aln ${path_aln} --ref ${ref_name} \
                                                               --path.processed ${genomes_out} --path.chr ${path_chr}
fi

pokaz_message "Script completed successfully!"









