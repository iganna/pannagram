show_help() {
    cat << EOF

╔════════════════════════════════════════╗
║   S e a r c h   f o r   s i m i l a r  ║
║            s e q u e n c e s           ║
╚════════════════════════════════════════╝

This script performs a BLAST search on a given FASTA file against a specified genome 
and processes the results based on similarity thresholds.

Usage: ${0##*/}  -in_seq FASTA_FILE -out OUTPUT_DIR [-aa]
                 [-on_seq [-strandfree] [-positions] |-on_genome|-on_path]
                 [-sim SIMILARITY_CUTOFF] [-cov COVERAGE_CUTOFF]
                 [-afterblast] [-keepblast]
                 [-cores CPU_THREADS]
                 [-h]

Mode (at least one is required):
    -on_seq SEQUENCE_FILE   Fasta-file containing sequences for comparison.
    -on_genome GENOME_FILE  Fasta-file containing genomes for comparison.
    -on_path SENOME_FOLDER  Path to the folder containing fasta-files with genomes.

Options:
    -in_seq FASTA_FILE     Path to the input FASTA file containing sequences to be processed.
    -out OUTPUT_DIR        Path to the output directory for output files and blast db.
    
    -aa AMINOACID_SEARCH   Switch to tblastn engine if your FASTA_FILE contains proteins (blastn otherwise);
    -sim SIMILARITY_CUTOFF Similarity % cutoff. (default: blastn: 85, tblastn: 60);
    -cov COVERAGE_CUTOFF   Coverage % cutoff. (default: same as -sim);
    -afterblast            Flag to process existing BLAST results;
    -keepblast             Flag to keep intermediate BLAST results;
    -strandfree            Use both strands for coverage. This option is used together with -on_seq;
    -positions             Also output the covered span (first..last covered position) in the query
                           and in the target (on its forward strand). Used together with -on_seq;
    -cores CPU_THREADS     Specify BLAST cores usege (default: 1);
    -h                     Show this help message and exit.

Examples:
    ${0##*/} -in_seq input.fasta -on_genome genome.fasta -out OUT_DIR
    ${0##*/} -in_seq input.fasta -on_genome genome.fasta -out OUT_DIR -sim 90 -keepblast
    ${0##*/} -in_seq input.fasta -on_genome genome.fasta -out OUT_DIR -sim 95 -afterblast 

    ${0##*/} -in_seq input.fasta -on_seq sequences.fasta -out OUT_DIR

    ${0##*/} -in_seq input.fasta  -on_path folder_with_genomes -out OUT_DIR -keepblast
    ${0##*/} -in_seq proteins.faa -on_path folder_with_genomes -out OUT_DIR -sim 50 -cov 90

EOF
}