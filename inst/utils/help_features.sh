print_usage() {
    cat << EOF
Usage: ${0##*/} -path_in PROJECT_DIR
                [-ref REF]
                [-h] [-cores NUM_CORES]  
                [-blocks] [-seq] [-aln] [-snp] 
                [-aln_type ALN_TYPE]


This script manages various genomic analyses and alignments.

Options:
    -h, --help                      Display this help message and exit;
    -cores NUM_CORES                Number of cores for parallel processing (default is 1);

    -path_in | -path_project PROJECT_DIR
                                    Path to pannagram (project) output directory;

    -ref REF                        Prefix for the gaccession, which was used to sort the alignment;
    -blocks                         RGet positions of synteny blocks between accessions;
    -seq                            Obtain consensus sequence for the pangenome alignment;
    -snp                            Get VCF file with SNPs;
    
    -sv                             SV calling;
    -sv_graph                       Create the Graph of SVs;


    -aln_type ALN_TYPE              Set the type of alignment (default: 'msa_');

Examples:
    ${0##*/}  -path_in '<project_dir>' -ref '<reference_name>' -blocks -seq -snp -sv -sv_graph -cores 4

EOF
}