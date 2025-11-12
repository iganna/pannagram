#!/bin/bash

if [ $# -eq 0 ]; then
    pokaz_error "No arguments provided!"
    help_in_box
    exit 0
fi

aln_type=''
ref_pref=''
run_blocks=false
run_seq=false
run_snp=false
run_snp_pi=false
run_sv_call=false
run_sv_sim=false
run_sv_sim_prot=false
run_sv_graph=false
run_annogroup=false
run_sv_orf=false
plot_families="F"

required_params=()

cmdline="$(basename "$0") $@"
while [ $# -gt 0 ]; do
    case "$1" in
        -h|-help) print_usage;                                          exit 0  ;;
        -cores)   cores="$2";                                           shift 2 ;;
        -path_in|-path_proj|-path_project) path_project="$2"; required_params+=("path_project");shift 2 ;;
        -ref)             ref_pref="$2";                                    shift 2 ;;
        -blocks| -synteny)          run_blocks=true;                                  shift 1 ;;
        -seq|-consensus)  run_seq=true;                                     shift 1 ;;
        -snp)             run_snp=true;                                     shift 1 ;;
        -snp_pi)          run_snp_pi=true;                                  shift 1 ;;
        -sv)              run_sv_call=true;                                 shift 1 ;;
        -sv_sim)          run_sv_sim=true;      set_file="$2";              shift 2 ;;
        -sv_sim_prot)     run_sv_sim_prot=true; set_file_prot="$2";         shift 2 ;;
        -sv_orf)          run_sv_orf=true;                                  shift 1 ;;
        -sv_graph | \
        -sv_families)     run_sv_graph=true;                                shift 1 ;;
        -plot_families)   plot_families="T";                                shift 1 ;;
        -sim)             similarity_value="$2";                            shift 2 ;;
        -sv_acc)          acc_anal="$2";                                    shift 2 ;;
        -annogroup)       run_annogroup=true;   path_annot="$2";            shift 2 ;;
        -aln_type)        aln_type="$2";                                    shift 2 ;;
        * ) pokaz_error "Unknown parameter: $1"; 
            pokaz_attention "!!! Please check your command: ${cmdline}";
            help_in_box; exit 1;;
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

if ! $run_blocks && \
   ! $run_seq && \
   ! $run_snp && \
   ! $run_snp_pi && \
   ! $run_sv_call && \
   ! $run_sv_sim && \
   ! $run_sv_sim_prot && \
   ! $run_sv_graph && \
   ! $run_annogroup && \
   ! $run_sv_orf; then
    pokaz_error "Warning: no features were provided for analysis." >&2
    exit 0
fi

# if [[ "$run_snp_pi" == true && "$run_snp" == false ]]; then
#     pokaz_error "Error: -snp_pi won't run unless -snp flag is used"
#     exit 1
# fi

cores="${cores:-1}"
acc_anal="${acc_anal:-NULL}"

path_project=$(add_symbol_if_missing "$path_project" "/")


# Setup the alignment type
if [ -z "$aln_type" ]; then
  aln_type="pan"  # Default
fi


if [ -z "$ref_pref" ]; then
  ref_pref="NULL"
else
    if [[ -z "$aln_type" || "$aln_type" == "ref" ]]; then
      aln_type="ref"
    else
      pokaz_error "Error: aln_type is set to '$aln_type', but should be ref"
      exit 1
    fi
fi




