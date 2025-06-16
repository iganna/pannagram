#!/bin/bash

if [ $# -eq 0 ]; then
    pokaz_error "No arguments provided!"
    help_in_box
    exit 0
fi

# Initialize variables
mode_pre="F"
mode_ref="F"
mode_msa="F"
clean="F"
one_step="F"
purge_reps="T"
rm_inter="F"
extra_steps="F"
purge_contigs="F"
unrecognized_options=()
step_start=0
step_end=100

# Required parameters
required_params=()

# Process command-line arguments
while [ $# -gt 0 ]; do
    case "$1" in
        -h|-help )       print_usage_detailed; print_examples; exit 0 ;;
        
        # Required parameters
        -path_out|-path_project) path_project="$2"; shift 2; required_params+=("path_project") ;;
        -path_in)                path_in="$2"; shift 2; required_params+=("path_in") ;; # path with all genomes in fasta format

        -s|-stage|-step) step_start="$2"; shift 2 ;; # first stage to run from, when the stage is not provided - the last interrupted stage withh be re-run
        -e|-end)         step_end="$2"; shift 2 ;;   # last stage to run to
        -log)            log_level="$2"; shift 2 ;;
        -cores)          cores="$2"; shift 2 ;;
        -clean|-cleanup) clean="T"; shift 1 ;;
        -one_step | -1s) one_step="T"; shift 1 ;;
        -rm_inter)       rm_inter="T"; shift 1 ;;
        
        
        # REF-based
        -pre)            mode_pre="T"; shift 1 ;;                # shitch to the preliminary mode
        -ref)            ref_name="$2"; mode_ref="T"; shift 2 ;; # name of the reference genome
        -path_ref)       path_ref="$2"; shift 2 ;;               # dont provide if it's the same folder as the path with query genomes
        
        # MSA
        -refs)           ref_set="$2"; mode_msa="T"; shift 2 ;;
        -nref)           ref_num="$2"; mode_msa="T"; shift 2 ;;
        
        # Number of chromosomes
        -nchr)           nchr="$2"; shift 2 ;;        # in every genome
        -nchr_ref)       nchr_ref="$2"; shift 2 ;;    # in reference genome

        -part_len)       part_len="$2"; shift 2 ;;    # fragments to which each chromosome should be cut, has a default value 5000
        -p_ident)        p_ident="$2"; shift 2 ;;     # percent of identity
        -p_ident_gap)    p_ident_gap="$2"; shift 2 ;; # percent of identity of gaps
        -max_len_gap)    max_len_gap="$2"; shift 2 ;; # Max length that can be aligned with MAFFT

        -one2one)        one2one="T"; shift 1 ;;                 # do compare chroms one-to-one (default in REF and MSA modes)
        -all2all)        one2one="F"; shift 1 ;;                 # do compare chroms all-to-all (default in PRE mode)
        -incl_reps |-include_repeats) purge_reps="F"; shift 1 ;; # repeats filtration
        -rev )           flag_rev="T"; shift 1 ;;                # reverse parts
        -orf )           flag_orf="T"; shift 1 ;;                # run ORF finder
        -purge_contigs)  purge_contigs="T"; shift 1 ;;
        -accessions)     acc_file="$2"; shift 2 ;;               # accessions file to analyze
        -combinations)   comb_file="$2"; shift 2 ;;              # chromosomal combinations tsv file to analyze: first column - query, second column - reference(base)
        -extra)          extra_steps="T"; shift 1 ;;             # perform extra steps
        *)               unrecognized_options+=("$1"); shift 1 ;;
    esac
done

# Validate required parameters
for param in "${required_params[@]}"; do
    if [[ -z "${!param}" ]]; then
        pokaz_error "Error: Required parameter -$param is missing!"
        help_in_box
        exit 1
    fi
done

# Handle unrecognized options
if [[ ${#unrecognized_options[@]} -gt 0 ]]; then
    pokaz_error "Unrecognized options: ${unrecognized_options[*]}"
    help_in_box
    exit 1
fi

# Determine mode_pangen
name_mode_pre='PRE'
name_mode_ref='REF'
name_mode_msa='MSA'

if [[ "$mode_pre" = "F" && "$mode_ref" = "F" && "$mode_msa" = "F" && -z "$path_ref" ]]; then
    mode_pangen="$name_mode_msa"
elif [[ "$mode_pre" = "T" && "$mode_ref" = "T" && "$mode_msa" = "F" ]]; then
    mode_pangen="$name_mode_pre"
elif [[ "$mode_ref" = "T" && "$mode_pre" = "F" && "$mode_msa" = "F" ]]; then
    mode_pangen="$name_mode_ref"
else
    pokaz_error "Invalid parameter combination for mode_pangen"
    exit 1
fi

# Validate one2one/all2all conflict
if [[ -n "$one2one" && "$one2one" != "T" && "$one2one" != "F" ]]; then
    pokaz_error "Error: -all2all and -one2one cannot be used together"
    exit 1
fi

path_in=$(add_symbol_if_missing "$path_in" "/")
path_project=$(add_symbol_if_missing "$path_project" "/")

path_ref="${path_ref:-${path_in}}"
path_ref=$(add_symbol_if_missing "$path_ref" "/")
