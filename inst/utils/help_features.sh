print_usage() {
    cat << EOF
Usage: ${0##*/}  -path_msa PATH_MSA  -path_chr PATH_CHR 
                [-ref REF]
                [-h] [-cores NUM_CORES]  
                [-blocks] [-seq] [-aln] [-snp] 
                [-aln_type ALN_TYPE]


This script manages various genomic analyses and alignments.

Options:
    -h, --help                  Display this help message and exit.
    -cores NUM_CORES            Specify the number of cores for parallel processing (default is 1).

    -path_msa PATH_MSA          Specify the global prefix for multiple sequence alignment. (PATH_OUT/intermediate/consensus/ from pannagram)
    -path_chr PATH_CHR          Specify the path to chromosome files. (PATH_OUT/intermediate/chromosomes/ from pannagram)

    -ref REF                    Specify the prefix for the gaccession, which was used to sort the alignment.
    -blocks                     RGet positions of synteny blocks between accessions.
    -seq                        Obtain consensus sequence for the pangenome alignment.
    -aln                        Produce a FASTA file with the pangenome alignment.
    -snp                        Get VCF file with SNPs.
    
    -sv                         SV calling
    -sv_graph                   Create the Graph of SVs
    -sv_sim                     Compare SVs with a dataset of sequences

    -annogroup                  Create the consensus annotation groups

    -aln_type ALN_TYPE          Set the type of alignment (default: 'msa_').

Examples:
    ${0##*/}  -path_msa /data/genomes -ref genome_ref -path_chr /data/chromosomes  -blocks -seq -snp

EOF
}