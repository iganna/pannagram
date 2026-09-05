print_usage() {
    cat << EOF
Usage: ${0##*/} -path_project PROJECT_DIR
                [-ref REF]
                [-h] [-cores NUM_CORES]
                [-synteny] [-consensus] [-aln] [-snp]
                [-aln_type ALN_TYPE]


This script manages various genomic analyses and alignments.

Options:
    -h, --help                      Display this help message and exit;
    -cores NUM_CORES                Number of cores for parallel processing (default is 1);

    -path_project PROJECT_DIR       Path to pannagram (project) output directory;

    -ref REF                        Prefix for the gaccession, which was used to sort the alignment;
    -synteny                        Get positions of synteny blocks between accessions;
    -consensus                      Obtain consensus sequence for the pangenome alignment;
    -snp                            Get VCF file with SNPs;

    -sv                             SV calling;
    -sv_families                    Create the Graph of SVs;
    -sv_te_order                    Structural order (LTR / TIR / poly-A) of SV families;


    -aln_type ALN_TYPE              Set the type of alignment (default: 'pan');

Examples:
    ${0##*/}  -path_project '<project_dir>' -ref '<reference_name>' -synteny -consensus -snp -sv -sv_families -cores 4

EOF
}
