#!/bin/bash

if [ $# -eq 0 ]; then
    pokaz_error "No arguments provided!"
    help_in_box
    exit 0
fi

aln_type='msa_'
run_blocks=false
run_seq=false
run_aln=false
run_snp=false
run_snp_pi=false
run_sv_call=false
run_sv_sim=false
run_sv_sim_prot=false
run_sv_graph=false
run_annogroup=false

required_params=()

while [ $# -gt 0 ]; do
    case "$1" in
        -h|-help) print_usage;                                          exit 0  ;;
        -cores)   cores="$2";                                           shift 2 ;;
        -path_in) path_in="$2";                                         shift 2 required_params+=("path_in") ;;
        -ref)         ref_pref="$2";                                    shift 2 ;;
        -blocks)      run_blocks=true;                                  shift 1 ;;
        -seq)         run_seq=true;                                     shift 1 ;;
        -aln)         run_aln=true;                                     shift 1 ;;
        -snp)         run_snp=true;                                     shift 1 ;;
        -snp_pi)      run_snp_pi=true;                                  shift 1 ;;
        -sv_call|-sv) run_sv_call=true;                                 shift 1 ;;
        -sv_sim)      run_sv_sim=true;      set_file="$2";              shift 2 ;;
        -sv_sim_prot) run_sv_sim_prot=true; set_file_prot="$2";         shift 2 ;;
        -sv_graph)    run_sv_graph=true;                                shift 1 ;;
        -sim)         similarity_value="$2";                            shift 2 ;;
        -sv_acc)      acc_anal="$2";                                    shift 2 ;;
        -annogroup)   run_annogroup=true;   path_annot="$2";            shift 2 ;;
        -aln_type)    aln_type="$2";                                    shift 2 ;;
        *)            pokaz_error "Unknown parameter: $1"; help_in_box; exit 1  ;;
    esac
done

# Validate required parameters
for param in "${required_params[@]}"; do
    if [[ -z "${!param}" ]]; then
        pokaz_error "Error: -$param is required"
        help_in_box
        exit 1
    fi
done

cores="${cores:-1}"
acc_anal="${acc_anal:-NULL}"
ref_pref="${ref_pref:-NULL}"

path_in=$(add_symbol_if_missing "$path_in" "/")

export aln_type run_blocks run_seq run_aln run_snp run_snp_pi run_sv_call run_sv_sim run_sv_sim_prot run_sv_graph run_annogroup
export cores acc_anal ref_pref path_in set_file set_file_prot similarity_value path_annot