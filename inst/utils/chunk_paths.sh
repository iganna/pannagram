path_log="${path_project}.logs/"

path_features="${path_project}features/"
path_features_msa="${path_features}msa/"
path_extra="$path_features/extra/"
path_snp="$path_features/snp/"
path_seq="$path_features/seq/"
path_sv="$path_features/sv/"
path_gff="$path_sv/gff/"

path_inter="${path_project}.intermediate/"
path_alignment="${path_inter}alignments/"
path_inter_msa="${path_inter}msa/"
path_blast="${path_inter}blast/"
path_mafft="${path_inter}mafft/"
path_chrom="${path_inter}chromosomes/"
path_parts="${path_chrom}parts/"
path_parts_mirror="${path_chrom}parts_mirror/"
path_orf="${path_inter}orf/"

path_plots="${path_project}plots/"
path_plots_snp="${path_plots}snp/"
path_plots_sv="${path_plots}sv/"
path_plots_synteny="${path_plots}synteny_pangenome/"
path_plots_pairwise="${path_plots}synteny_pairwise/"


function check_dir() {
  local dir="$1"

  if [[ -z "$dir" ]]; then
    pokaz_error "Missing check_dir only argument" >&2
    return 1
  fi

  if [[ ! -d "$dir" ]]; then
    pokaz_error "Expected '$dir' does not exist or is not a directory." >&2
    return 1
  fi

  if [[ ! -w "$dir" ]]; then
    pokaz_error "Alas! '$dir' is not writable." >&2
    return 1
  fi

  return 0
}