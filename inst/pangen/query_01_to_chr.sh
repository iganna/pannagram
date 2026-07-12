#!/bin/bash
# Fast bash replacement for query_01_to_chr.R (STEP 1: genomes -> chromosomes).
# For each accession it streams the genome FASTA in a single awk pass and writes
# the first <n.chr> sequences as <path.out>/<acc>_chr<i>.fasta (header
# ">${acc}_Chr<i>", sequence on ONE line) plus a length table
# <path.out>/<acc>_chr_len.txt with columns "acc chr len name".
#
# Sequences are NOT upper-cased: toupper is not used downstream and it was the
# whole cost of the old R step (~10x slower). To restore case-folding, change
# `s=$0` to `s=toupper($0)` in the awk sequence-line rule.
#
# Accepts the same options as query_01_to_chr.R (extra ones are ignored).
set -uo pipefail

path_in=""; path_out=""; cores=1; nchr=0; file_acc=""
while [ $# -gt 0 ]; do
  case "$1" in
    --path.in)       path_in="$2";       shift 2;;
    --path.out)      path_out="$2";       shift 2;;
    --cores)         cores="$2";          shift 2;;
    --n.chr)         nchr="$2";           shift 2;;
    --accessions)    file_acc="$2";       shift 2;;
    --path.log)      shift 2;;            # accepted, not needed here
    --log.level)     shift 2;;
    --f.chr.anal)    shift 2;;
    --purge.contigs) shift 2;;            # unused in the R too
    *)               shift;;
  esac
done

# ensure trailing slash on the directories
[ "${path_in: -1}" = "/" ]  || path_in="${path_in}/"
[ "${path_out: -1}" = "/" ] || path_out="${path_out}/"
mkdir -p "$path_out"

# Split one accession's genome into per-chromosome files + length table.
split_one() {
  local acc="$1" path_in="$2" path_out="$3" nchr="$4"
  local genome="" ext
  for ext in fasta fna fa fas; do
    if [ -f "${path_in}${acc}.${ext}" ]; then genome="${path_in}${acc}.${ext}"; break; fi
  done
  if [ -z "$genome" ]; then echo "query_01_to_chr: genome for '${acc}' not found in ${path_in}" >&2; return 0; fi
  local n="$nchr"; [ "$n" = "0" ] && n=2147483647
  awk -v acc="$acc" -v nchr="$n" -v outdir="$path_out" -v lenfile="${path_out}${acc}_chr_len.txt" '
  /^>/ {
    if (prevout != "") { printf "\n" > prevout; prevout="" }
    idx++
    if (idx <= nchr) {
      name = substr($0, 2); gsub(/ /, "_", name)
      names[idx] = name; len[idx] = 0
      cur = sprintf("%s%s_chr%d.fasta", outdir, acc, idx)
      printf ">%s_Chr%d\n", acc, idx > cur
      prevout = cur
    }
    next
  }
  { if (idx >= 1 && idx <= nchr) { printf "%s", $0 > cur; len[idx] += length($0) } }
  END {
    if (prevout != "") printf "\n" > prevout
    print "acc\tchr\tlen\tname" > lenfile
    n = (idx < nchr ? idx : nchr)
    for (i = 1; i <= n; i++) printf "%s\t%d\t%d\t%s\n", acc, i, len[i], names[i] >> lenfile
  }
  ' "$genome"
}
export -f split_one

# One accession per task, in parallel across cores.
if command -v parallel >/dev/null 2>&1; then
  grep -v '^[[:space:]]*$' "$file_acc" | parallel --will-cite -j "$cores" split_one {} "$path_in" "$path_out" "$nchr"
else
  while IFS= read -r acc; do
    [ -n "$acc" ] && split_one "$acc" "$path_in" "$path_out" "$nchr"
  done < "$file_acc"
fi
