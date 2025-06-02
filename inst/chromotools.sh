#!/bin/bash

# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------
INSTALLED_PATH=$(Rscript -e "cat(system.file(package = 'pannagram'))")
source $INSTALLED_PATH/utils/chunk_error_control.sh

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source $INSTALLED_PATH/utils/utils_bash.sh


print_usage() {
    cat << EOF
Usage: ${0##*/}  -path_in PATH_IN  [-path_out PATH_OUT]
                [-remove WORDS_REMOVE] [-remain WORDS_REMAIN]
                [-cores NUM_CORES] [-h]

OR

Usage: ${0##*/}  -path_in PATH_IN  [-path_out PATH_OUT]
                -path_aln PATH_ALN [-resort] [-rearrange]
                [-cores NUM_CORES] [-h]


This script processes genomes either by filtering based on keywords 
or by reordering and rearranging chromosomes according to a given alignment.

Options:
    -h, --help                  Display this help message and exit.
    -cores NUM_CORES            Number of cores for parallel processing (default: 1).

Input/Output:
    -path_in PATH_IN            Path to the input directory with genomes or chromosomes depends on the mode.
    -path_out PATH_OUT          Path to the output directory for processed genomes.
                                Default: PATH_IN/processed/

Filtering mode (before alignment):
    -remove WORDS_REMOVE        Comma-separated list of words to filter out.
    -remain WORDS_REMAIN        Comma-separated list of words to retain.

Alignment-based mode (after alignment):
    -path_aln PATH_ALN          Path to the alignment file used to reorder genomes.
    -resort                     Reorder chromosomes according to the alignment.
    -rearrange                  Additionally rearrange (split/merge) chromosomes 
                                corresponding to the reference genome.

Examples:

# Example 1: Filtering genomes based on keywords (before alignment)
${0##*/} -path_in /data/genomes -path_out /data/genomes/filtered \\
         -remove "contig,mitochondria" -remain "chromosome" -cores 4


# Example 2: Reordering genomes based on alignment (after alignment)
${0##*/} -path_in /data/genomes -path_out /data/genomes/reordered \\
         -path_aln /data/alignment/aln.fasta -resort -cores 4


# Example 3: Rearranging genomes based on alignment (after alignment)
${0##*/} -path_in /data/chromosomes -path_out /data/genomes/rearrange \\
         -path_aln /data/alignments_reference/aln.fasta -rearrange -cores 4

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
path_in=""
path_out=""
words_remove=""
words_remain=""
words_keep=""
path_aln=""
resort=false
rearrange=false

while [ $# -gt 0 ]; do
    case $1 in
        -h|--help)        print_usage; exit 0 ;;
        -cores)           require_arg $1 $2; cores=$2; shift 2 ;;
        -path_in)         require_arg $1 $2; path_in=$2; shift 2 ;;
        -path_out)        require_arg $1 $2; path_out=$2; shift 2 ;;
        -remove)          require_arg $1 $2; words_remove=$2; shift 2 ;;
        -remain)          require_arg $1 $2; words_remain=$2; shift 2 ;;
        -keep)            require_arg $1 $2; words_keep=$2; shift 2 ;;
        -path_aln)        require_arg $1 $2; path_aln=$2; shift 2 ;;
        -resort)          resort=true; shift ;;
        -rearrange)       rearrange=true; shift ;;
        *)                pokaz_error "Unknown parameter: $1"; help_in_box; exit 1;;
    esac
done


# Check required parameters
if [[ -z "$path_in" ]]; then
    pokaz_error "Error: -path_in is required"
    help_in_box
    exit 1
fi
path_in=$(add_symbol_if_missing "$path_in" "/")

if [[ -z "$path_out" ]]; then
    path_out="${path_in}filtered/"
fi
path_out=$(add_symbol_if_missing "$path_out" "/")
mkdir -p "$path_out"

# Determine mode
mode_filter=false
mode_reorder=false
mode_rearrange=false

# Filtering mode (before alignment)
if [[ -n "$words_remove" || -n "$words_remain" ]]; then
    mode_filter=true
fi

# Reordering mode (after alignment, reorder chromosomes)
if [[ -n "$path_aln" && "$resort" == true && "$rearrange" == false ]]; then
    mode_reorder=true
fi

# Rearrangement mode (after alignment, split/merge chromosomes)
if [[ -n "$path_aln" && "$rearrange" == true && "$resort" == false ]]; then
    mode_rearrange=true
fi

# Check invalid combinations
if [[ "$mode_filter" == true && ( "$mode_reorder" == true || "$mode_rearrange" == true ) ]]; then
    pokaz_error "Error: Filtering mode (-remove/-remain) cannot be combined with reorder/rearrange modes."
    help_in_box
    exit 1
fi

if [[ "$mode_filter" == false && "$mode_reorder" == false && "$mode_rearrange" == false ]]; then
    pokaz_error "Error: Specify parameters for one of the modes:
    - Filtering:        -remove/-remain
    - Reordering:       -path_aln -resort
    - Rearrangement:    -path_aln -rearrange"
    help_in_box
    exit 1
fi

pokaz_message "Number of cores: ${cores}"

# ----------------------------------------------------------------------------
#                                   MAIN
# ----------------------------------------------------------------------------

# -------------------------------------------------
if [ "$mode_filter" = true ]; then
    pokaz_stage "Filtering chromosomes by keywords."

    path_plots="${path_consensus}plot_synteny/"
    mkdir -p ${path_plots}

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
        --path.genomes ${path_in} \
        --path.filtered ${path_out} \
        --cores ${cores} \
        ${param_words_remove} \
        ${param_words_remain} \
        ${param_words_keep}

fi


# -------------------------------------------------
if [ "$mode_reorder" = true ]; then
    pokaz_stage "Reordering chromosomes based on the preliminary alignment."

    path_aln=$(add_symbol_if_missing "$path_aln" "/")
    path_in=$(add_symbol_if_missing "$path_in" "/")
    path_out=$(add_symbol_if_missing "$path_out" "/")

    # Get the name of the reference genome
    ref_name=$(basename "$path_aln")
    ref_name=${ref_name#alignments_}

    # Run reordering
    Rscript $INSTALLED_PATH/chromotools/resort_01_find_best.R --path.aln ${path_aln} --ref ${ref_name} --path.resort ${path_out}
    Rscript $INSTALLED_PATH/chromotools/resort_02_resort.R --path.genomes ${path_in} --ref ${ref_name} --path.resort ${path_out}
fi


# -------------------------------------------------
if [ "$mode_rearrange" = true ]; then

    pokaz_stage "Find the best rearrangements based on alignment..."

    path_aln=$(add_symbol_if_missing "$path_aln" "/")
    path_in=$(add_symbol_if_missing "$path_in" "/")
    path_out=$(add_symbol_if_missing "$path_out" "/")

    # Get the name of the reference genome
    ref_name=$(basename "$path_aln")
    ref_name=${ref_name#alignments_}

    Rscript $INSTALLED_PATH/chromotools/rearrange_01_positions.R --path.aln ${path_aln} --ref ${ref_name} --path.processed ${path_out} --path.chr ${path_in}

    pokaz_stage "Rearranging (splitting/merging) chromosomes..."
    Rscript $INSTALLED_PATH/chromotools/rearrange_02_genomes.R --path.aln ${path_aln} --ref ${ref_name} --path.processed ${path_out} --path.chr ${path_in}
fi

pokaz_message "Script completed successfully!"









