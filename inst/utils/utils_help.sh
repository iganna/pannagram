# Function to display a short help message
print_usage_short() {
    pokaz_help

    cat << EOF

Show short help:
${0##*/} -h

Show detailed help:
${0##*/} -help

Usage for the PREliminary alignment:
${0##*/} -path_in INPUT_FOLDER -path_out OUTPUT_FOLDER -pre 
         -ref REF_NAME 
         [-nchr NUM_OF_CHR] 
         [-nchr_ref NUM_OF_CHR_IN_REF] 
         [-one2one]
         [OTHER PARAMETERS]

Usage for the REF-base alignment:
${0##*/} -path_in INPUT_FOLDER -path_out OUTPUT_FOLDER 
         -ref REF_NAME 
         -nchr N_CHR_QUERY 
         [-nchr_ref N_CHR_REF]
         [-all2all]
         [OTHER PARAMETERS]

Usage for the Multiple Genome Alignment (MGA):
${0##*/} -path_in INPUT_FOLDER -path_out OUTPUT_FOLDER
        [-refs REFS_NAMES]
        [-nref NUM_OF_REFS]
        [-nchr NUM_OF_CHRS]
        [-all2all]
        [OTHER PARAMETERS]

OTHER PARAMETERS:
    [-h] [-s STAGE] [-cores CORES] [-log LOG_LEVEL]
    [-sort_len] [-purge_repeats]
    [-path_ref PATH_CHR_REF]
    [-p_ident P_IDENT] [-part_len PART_LEN]
    [-accessions ACC_FILE] [-combinations COMB_FILE]
    
EOF
}

#[-path_chrom PATH_CHROM] [-path_parts PATH_PARTS]
# -path_chrom PATH_CHROM      Path to the folder with individual chromosomes in separate files. 
# -path_parts PATH_PARTS      Path to the folder with files of chromosomal parts.

# Function to display a long help message
print_usage_detailed() {
    pokaz_help

    cat << EOF
This script performs alignment of query genomes to the reference genome.

Options:
    -h, -help                   Display this help messages and exit.
    -s, -stage STAGE            Specify the stage from which to run: a number from 1 to 12.
                                If not provided, the last interrupted stage will be re-run.
    -cores CORES                Number of cores for parallel processing. Default is 1.
    -log LOG_LEVEL              The level of logging to be shown on the screen. 
                                Options: 
                                    0 - off
                                    1 - step names (default)
                                    2 - step progress
                                    3 - full logging (visible only with "-cores 1")
        
    * Required parameters:
        -path_in INPUT_FOLDER       Folder to the directory with query genomes. 
                                    Possible file types: fasta, fna, fa, fas.
        -path_out OUTPUT_FOLDER     Folder where all results and intermediate files will appear.


    * Optional parameters for REF-based alignments:
        -pre                        Flag for only the "PRElimilary" stages
        -ref REF_NAME               Name of the reference genome: REF_NAME.fasta.
                                    Names should NOT contain "_chr" as a substring.
        -nchr_ref N_CHR_REF         Number of chromosomes in the reference genome.
        -nchr_acc N_CHR_QUERY       Number of chromosomes in the query genome.

    * Optional parameters for MSA:
        -refs REFS_NAMES            Names of reference genomes for randomisation.
                                    Names should NOT contain "_chr" as a substring.
        -nref NUM_OF_REFS           Number of reference genomes for randomising,
                                    when -refs is not specified.
        -nchr NUM_OF_CHRS           Number of chromosomes in both reference and query genomes.

    * Optional paths:
        -path_ref PATH_CHR_REF      Path where the reference genome is stored. 
                                    Do not provide if it's the same folder as the path with query genomes.

    * Input Design Handling:
        -sort_len                   Flag to sort chromosomes by length.
        -one2one                    Flag for straightforward pairwise alignment of chromosomes,
                                    matching each chromosome sequentially with its corresponding one 
                                    This is default: for REF and MSA mods.
        -all2all                    Flag for aligning all chromosomes to all.
                                    This is default for PRE mode.
        -accessions ACC_FILE        File with accessions to analyze. Accessions should be in rows.
        -combinations COMB_FILE     File with combinations to analyze.

    * Tuning parameters: 
        -purge_repeats              Flag for filtering of repeats (default is no filtering).
        -p_ident P_IDENT            Percentage identity threshold (default: 85).
        -part_len PART_LEN          Fragments to which each chromosome should be cut (default value: 5000).
        

EOF
}

# Function to display examples
print_examples() {
    cat << EOF
The simplest examples:
* MSA:
    ${0##*/} -path_in 'input_genomes' -path_out 'output_folder'
* REF-base:
    ${0##*/} -path_in 'input_genomes' -path_out 'output_folder' -ref 'ref_name'
* REF-base preliminary run:
    ${0##*/} -path_in 'input_genomes' -path_out 'output_folder' -ref 'ref_name' -pre

EOF
}
