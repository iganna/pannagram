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
  echo "-path_chrom"
  echo "-path_parts"
  echo "-path_result"
  echo "-ref"
  echo "-one2one"
  echo "-p_ident"
  echo "-cores"
  echo "-penalty"
  echo "-gapopen"
  echo "-gapextend"
  echo "-max_hsps"
  echo "-word_size"
  echo "-log_path"
}

# ----------------------------------------------------------------------------
#                 PARAMETERS
# ----------------------------------------------------------------------------

all_vs_all="T"

while [ $# -gt 0 ]
do
    case $1 in
    # for options with required arguments, an additional shift is required
    -path_chrom) path_chrom=$2; shift 2;;
    -path_parts) path_parts=$2; shift 2;;
    -path_result) path_blast=$2; shift 2;;
    -ref) ref_name=$2; shift 2;;
    -p_ident) p_ident=$2; shift 2;;
    -cores) cores=$2; shift 2;;
    -penalty) penalty=$2; shift 2;;
    -gapopen) gapopen=$2; shift 2;;
    -gapextend) gapextend=$2; shift 2;;
    -xdrop_gap) xdrop_gap=$2; shift 2;;
    -xdrop_gap_final) xdrop_gap_final=$2; shift 2;;
    -max_hsps) max_hsps=$2; shift 2;;
    -word_size) w_size=$2; shift 2;;
    -path_log) path_log=$2; shift 2;;
    -combinations) file_combinations=$2; shift 2;;
    -accessions) file_accessions=$2; shift 2;;
    *) print_usage
       echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    esac
done

# -penalty -2 -gapopen 10 -gapextend 2 -max_hsps 5

penalty="${penalty:--2}"
gapopen="${gapopen:-10}"
gapextend="${gapextend:-2}"
# X-dropoff for gapped extension: low values stop extension at divergent patches,
# preventing ragged "staircase" gap runs and speeding up (default blastn: 30 / 100).
xdrop_gap="${xdrop_gap:-15}"
xdrop_gap_final="${xdrop_gap_final:-30}"
max_hsps="${max_hsps:-1}"
cores="${cores:-30}"
p_ident="${p_ident:-85}"
w_size="${w_size:-28}"  # parts-blast word_size fallback (tuned 2026-07-10 to 28; coverage-insensitive here). Normally supplied by pannagram.sh as -word_size ${w_size}; see the sweep table in pannagram.sh.

mkdir -p $path_blast

# ----------------------------------------------------------------------------
#                 FUNCTION
# ----------------------------------------------------------------------------

export penalty
export gapopen
export gapextend
export xdrop_gap
export xdrop_gap_final
export max_hsps
export cores
export p_ident
export w_size
export blastres
export all_vs_all
export log_path

# BLAST-search function
run_blast() {
    file_acc=$1
    file_ref=$2
    file_out=$3
    file_log=$4

    # If log file exists and has the word "Done" - then don't run the blast again
    if [ -f "$file_log" ]; then
        if grep -q "Done" "$file_log"; then
            echo "Over." >> "$file_log"
            return 0
        fi
    fi
    echo "New attempt:" > "$file_log"  # Create or empty the log file

    echo ${file_acc} >> ${file_log}
    echo ${file_ref} >> ${file_log}
    echo ${file_out} >> ${file_log}

    if [ ! -f ${file_acc} ]; then
        echo "File does not exist ${file_acc}"
        exit 1
        # return 0
    fi

    if [ ! -f ${file_ref} ]; then
        echo "File does not exist ${file_ref}"
        exit 1
        # return 0
    fi

    # Run BLAST
    if [ -f "${file_out}" ]; then 
      rm "${file_out}"    # Clean up the file
    fi
    blastn -db "${file_ref}" -query "${file_acc}" -out "${file_out}" \
           -outfmt "6 qseqid qstart qend sstart send pident length qseq sseq sseqid" \
           -perc_identity "${p_ident}" -penalty "$penalty" -gapopen "$gapopen" -gapextend "$gapextend" \
           -xdrop_gap "$xdrop_gap" -xdrop_gap_final "$xdrop_gap_final" \
           -max_hsps "$max_hsps" -word_size "$w_size" \
           >> "$file_log" 2>&1 # classic blastn + tuned gap penalties + low X-dropoff (15/30): clean alignments (no ragged gap staircases) and faster via early termination at divergent patches. Pair with -part_len 1000 for coverage. -word_size 50

    echo "Done." >> "$file_log"
}

export -f run_blast

# ----------------------------------------------------------------------------
#                 MAIN
# ----------------------------------------------------------------------------


# Count the number of chromosomes for every accession and the reference genome
# Initialize arrays to store accessions and their counts
accessions=()
acc_counts=()

# Initialize variable to store the count for the reference name
ref_count=0

# Read accessions from the file and process each one
while IFS= read -r accession; do
    # Skip the reference accession
    if [[ "$accession" != "$ref_name" ]]; then
        # Count files in ${path_parts} that start with the accession name
        count=$(ls "${path_parts}" | grep "^${accession}.*_chr.*\.fasta$" | wc -l)
        accessions+=("$accession")
        acc_counts+=("$count")
    fi
done < "${file_accessions}"

# Count files in ${path_chrom} that start with ${ref_name}
ref_count=$(ls "${path_chrom}" | grep "^${ref_name}.*_chr.*\.fasta$" | wc -l)

# # Initialize arrays to store file paths
# files_ref=()
# files_acc=()
# files_out=()
# files_log=()

# # Check if combinations file exists and is not empty
# if [ ! -s "${file_combinations}" ]; then
#     # All-to-all
#     for i in "${!accessions[@]}"; do
#         acc_name="${accessions[$i]}"
#         count="${acc_counts[$i]}"
#         for ((i_acc=1; i_acc<=count; i_acc++)); do
#             for ((i_ref=1; i_ref<=ref_count; i_ref++)); do
#                 files_ref+=("${path_chrom}${ref_name}_chr${i_ref}.fasta")
#                 files_acc+=("${path_parts}${acc_name}_chr${i_acc}.fasta")
#                 files_out+=("${path_blast}${acc_name}_${i_acc}_${i_ref}.txt")
#                 files_log+=("${path_log}${acc_name}_${i_acc}_${i_ref}.txt")

#                 # echo "${acc_name} ${i_acc} ${ref_name} ${i_ref}"
#             done
#         done
#     done
# else
#     # Only combinations from file_combinations for all accessions
#     while IFS= read -r acc_name; do
#         if [[ "$acc_name" != "$ref_name" ]]; then
#             # Read combinations from the file
#             while IFS=$'\t' read -r i_acc i_ref; do
#                 files_ref+=("${path_chrom}${ref_name}_chr${i_ref}.fasta")
#                 files_acc+=("${path_parts}${acc_name}_chr${i_acc}.fasta")
#                 files_out+=("${path_blast}${acc_name}_${i_acc}_${i_ref}.txt")
#                 files_log+=("${path_log}${acc_name}_${i_acc}_${i_ref}.txt")

#                 # echo "${acc_name} ${i_acc} ${ref_name} ${i_ref}"
#             done < "${file_combinations}"
#         fi
#     done < "${file_accessions}"
# fi

# # Run BLAST in parallel
# parallel --will-cite -j $cores --link run_blast ::: "${files_acc[@]}" ::: "${files_ref[@]}" ::: "${files_out[@]}" ::: "${files_log[@]}"


# Initialize arrays to store file paths
files_ref=()
files_acc=()
files_out=()
files_log=()

# Temporary files for combinations
temp_acc="${path_blast}temp_acc.txt"
temp_ref="${path_blast}temp_ref.txt"
temp_out="${path_blast}temp_out.txt"
temp_log="${path_blast}temp_log.txt"

> "$temp_acc"  # Clear or create the temporary files
> "$temp_ref"
> "$temp_out"
> "$temp_log"

# Check if combinations file exists and is not empty
if [ ! -s "${file_combinations}" ]; then
    # All-to-all
    for i in "${!accessions[@]}"; do
        acc_name="${accessions[$i]}"
        count="${acc_counts[$i]}"
        for ((i_acc=1; i_acc<=count; i_acc++)); do
            for ((i_ref=1; i_ref<=ref_count; i_ref++)); do
                echo "${path_chrom}${ref_name}_chr${i_ref}.fasta" >> "$temp_ref"
                echo "${path_parts}${acc_name}_chr${i_acc}.fasta" >> "$temp_acc"
                echo "${path_blast}${acc_name}_${i_acc}_${i_ref}.txt" >> "$temp_out"
                echo "${path_log}${acc_name}_${i_acc}_${i_ref}.txt" >> "$temp_log"
            done
        done
    done
else
    # Only combinations from file_combinations for all accessions
    while IFS= read -r acc_name; do
        if [[ "$acc_name" != "$ref_name" ]]; then
            # Read combinations from the file
            while IFS=$'\t' read -r i_acc i_ref; do
                echo "${path_chrom}${ref_name}_chr${i_ref}.fasta" >> "$temp_ref"
                echo "${path_parts}${acc_name}_chr${i_acc}.fasta" >> "$temp_acc"
                echo "${path_blast}${acc_name}_${i_acc}_${i_ref}.txt" >> "$temp_out"
                echo "${path_log}${acc_name}_${i_acc}_${i_ref}.txt" >> "$temp_log"
            done < "${file_combinations}"
        fi
    done < "${file_accessions}"
fi

# ----------------------------------------------------------------------------
# Chunked parallel BLAST.
# Previously one BLAST task per (query_chr, ref_chr) pair -> at most nchr tasks,
# so only nchr cores were ever busy (e.g. 5 of 16). We split each query parts
# file into sub-chunks of `parts_per_chunk` sequences, make one task per chunk,
# and let GNU parallel fill all cores. Chunk outputs are concatenated back into
# the original per-combination file, so the result is identical to before.
# ----------------------------------------------------------------------------
parts_per_chunk="${parts_per_chunk:-5000}"
chunk_dir="${path_blast}chunks/"
mkdir -p "$chunk_dir"

ctask_acc="${path_blast}ctask_acc.txt"
ctask_ref="${path_blast}ctask_ref.txt"
ctask_out="${path_blast}ctask_out.txt"
ctask_log="${path_blast}ctask_log.txt"
> "$ctask_acc"; > "$ctask_ref"; > "$ctask_out"; > "$ctask_log"

# Build chunk-level task lists from the per-combination lists
paste -d$'\t' "$temp_acc" "$temp_ref" "$temp_out" "$temp_log" | while IFS=$'\t' read -r qf rf of lf; do
    base=$(basename "$of" .txt)
    # split query parts fasta into <chunk_dir>/<base>.NNNN.fasta (parts_per_chunk seqs each)
    awk -v pre="${chunk_dir}${base}." -v cs="$parts_per_chunk" '
        /^>/ { if(n % cs == 0){ if(fh) close(fh); fh = sprintf("%s%04d.fasta", pre, int(n/cs)) } n++ }
        { print > fh }' "$qf"
    for cf in "${chunk_dir}${base}."*.fasta; do
        [ -f "$cf" ] || continue
        cb=$(basename "$cf" .fasta)
        echo "$cf"                   >> "$ctask_acc"
        echo "$rf"                   >> "$ctask_ref"
        echo "${chunk_dir}${cb}.txt" >> "$ctask_out"
        echo "${chunk_dir}${cb}.log" >> "$ctask_log"
    done
done

# Run all chunks in parallel across all cores
parallel --will-cite -j $cores --link run_blast :::: "$ctask_acc" :::: "$ctask_ref" :::: "$ctask_out" :::: "$ctask_log"

# Concatenate chunk outputs back into the per-combination result file
while IFS= read -r of; do
    base=$(basename "$of" .txt)
    cat "${chunk_dir}${base}."*.txt > "$of" 2>/dev/null
done < "$temp_out"

# Clean up temporary files
rm -rf "$chunk_dir"
rm -f "$temp_acc" "$temp_ref" "$temp_out" "$temp_log" \
      "$ctask_acc" "$ctask_ref" "$ctask_out" "$ctask_log"











