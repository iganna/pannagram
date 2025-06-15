show_help() {
    cat << EOF

╔════════════════════════════════════════╗
║   S e a r c h   f o r   s i m i l a r  ║
║            s e q u e n c e s           ║
╚════════════════════════════════════════╝

This script performs a BLAST search on a given FASTA file against a specified genome 
and processes the results based on similarity thresholds.

Usage: ${0##*/}  -in_seq FASTA_FILE -out OUTPUT_FILE
                [mode_option] [options]

Mode (at least one option is required:
    -on_seq SEQUENCE_FILE   Fasta-file containing sequences for comparison.
    -on_genome GENOME_FILE  Fasta-file containing genomes for comparison.
    -on_path SENOME_FOLDER  Path to the folder containing fasta-files with genomes.

Options:
    -in_seq FASTA_FILE     Path to the input FASTA file containing DNA sequences to be processed.
    -out OUTPUT_FILE       Path to the output file where the processed results will be saved.
    
    -sim SIMILARITY_CUTOFF (Optional) Similarity cutoff for sequence comparison.
    -cov SEQUENCE_COVERAGE (Optional) Cutoff for the coverage.
    -afterblast            (Optional) Flag to process existing BLAST results.
    -keepblast             (Optional) Flag to keep intermediate BLAST results.
    -strandfree            (Optional) Use both strands for coverage. This option is used together with -on_seq.
    -h                     Show this help message and exit.

Examples:
    ${0##*/} -in_seq input.fasta -on_genome genome.fasta -out out/
    ${0##*/} -in_seq input.fasta -on_genome genome.fasta -out out_90.txt -sim 90 -keepblast
    ${0##*/} -in_seq input.fasta -on_genome genome.fasta -out out_95.txt -sim 95 -afterblast 

    ${0##*/} -in_seq input.fasta -on_seq sequences.fasta -out out/

    ${0##*/} -in_seq input.fasta -on_path folder_with_genomes -out out_90.txt -sim 90 -keepblast

EOF
}