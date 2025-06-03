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
        -h| -help )       print_usage_detailed; print_examples; exit 0 ;;
        -s|-stage|-step ) step_start="$2"; shift 2 ;;
        -e | -end )       step_end="$2"; shift 2 ;;
        -log)             log_level="$2"; shift 2 ;;
        -cores)           cores="$2"; shift 2 ;;
        -clean | -cleanup)clean="T"; shift 1 ;;
        -one_step | -1s)  one_step="T"; shift 1 ;;
        -rm_inter)        rm_inter="T"; shift 1 ;;
        
        # Required parameters
        -path_out)        path_out="$2"; shift 2; required_params+=("path_out") ;;
        -path_in)         path_in="$2"; shift 2; required_params+=("path_in") ;;
        
        # REF-based
        -pre)             mode_pre="T"; shift 1 ;;
        -ref)             ref_name="$2"; mode_ref="T"; shift 2 ;;
        -path_ref)        path_ref="$2"; shift 2 ;;
        
        # MSA
        -refs)            ref_set="$2"; mode_msa="T"; shift 2 ;;
        -nref)            ref_num="$2"; mode_msa="T"; shift 2 ;;
        
        # Number of chromosomes
        -nchr)            nchr="$2"; shift 2 ;;
        -nchr_ref)        nchr_ref="$2"; shift 2 ;;

        -part_len)        part_len="$2"; shift 2 ;;
        -p_ident)         p_ident="$2"; shift 2 ;;
        -p_ident_gap)     p_ident_gap="$2"; shift 2 ;;
        -max_len_gap)     max_len_gap="$2"; shift 2 ;;
        -one2one)         one2one="T"; shift 1 ;;
        -all2all)         one2one="F"; shift 1 ;;
        -incl_reps |-include_repeats ) purge_reps="F"; shift 1 ;;
        -rev )            flag_rev="T"; shift 1 ;;
        -orf )            flag_orf="T"; shift 1 ;;
        -purge_contigs)   purge_contigs="T"; shift 1 ;;
        -accessions)      acc_file="$2"; shift 2 ;;
        -combinations)    comb_file="$2"; shift 2 ;;
        -extra)           extra_steps="T"; shift 1 ;;
        *)                unrecognized_options+=("$1"); shift 1 ;;
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

# Export variables for main script
export mode_pre mode_ref mode_msa clean one_step purge_reps rm_inter extra_steps purge_contigs
export step_start step_end log_level cores path_out path_in path_ref ref_name ref_set ref_num
export nchr nchr_ref part_len p_ident p_ident_gap max_len_gap one2one flag_rev flag_orf acc_file comb_file
export name_mode_pre name_mode_ref name_mode_msa mode_pangen