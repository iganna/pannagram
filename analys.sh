#!/bin/bash

# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

# Exit immediately if any command returns a non-zero status
set -e

# Keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# Define a trap for the EXIT signal
trap 'catch $?' EXIT

# Function to handle the exit signal
catch() {
    # Check if the exit code is non-zero
    if [ $1 -ne 0 ]; then
        echo "\"${last_command}\" command failed with exit code $1."
    fi
}

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

source utils/utils_bash.sh


print_usage() {
    cat << EOF
Usage: ${0##*/}  -path_msa PATH_MSA  -ref REF -path_chr PATH_CHR 
                [-h] [-cores NUM_CORES]  
                [-blocks] [-seq] [-aln] [-snp] 
                [-aln_type ALN_TYPE] [-path_cons PATH_CONS]


This script manages various genomic analyses and alignments.

Options:
    -h, --help                  Display this help message and exit.
    -cores NUM_CORES            Specify the number of cores for parallel processing (default is 1).

    -path_msa PATH_MSA          Specify the global prefix for multiple sequence alignment. The same as -path_out in pangen.sh
    -ref REF                    Specify the prefix for the gaccession, which was used to sort the alignment.
    -path_chr PATH_CHR          Specify the path to chromosome files.

    -blocks                     RGet positions of synteny blocks between accessions.
    -seq                        Obtain consensus sequence for the pangenome alignment.
    -aln                        Produce a FASTA file with the pangenome alignment.
    -snp                        Get VCF file with SNPs.
    

    -sv_call                    SV calling
    -sv_sim                     Compare SVs with a dataset of sequences
    -sv_graph                   Create the Graph of SVs


    -aln_type ALN_TYPE          Set the type of alignment (default: 'msa_').
    -path_cons PATH_CONS        Specify the path to the consensus folder (has the default value).

Examples:
    ${0##*/}  -path_msa /data/genomes -ref genome_ref -path_chr /data/chromosomes  -blocks -seq -snp

EOF
}

# ----------------------------------------------------------------------------
#             PARAMETERS
# ----------------------------------------------------------------------------


# Initialize variables to determine which scripts to run
aln_type='msa_'

# Simple analysis
run_blocks=false
run_seq=false
run_aln=false
run_snp=false

# analysis of SVs
run_sv_call=false
run_sv_sim=false
run_sv_graph=false

# Parse command line arguments
while [ $# -gt 0 ]; do
    case $1 in
        -h|-help) print_usage; exit 0;;
        -cores) cores=$2; shift 2 ;;

        -path_msa) path_consensus=$2; shift 2;;
        -ref) ref_pref=$2; shift 2;;
        -path_chr) path_chromosomes=$2; shift 2 ;;
        
        -blocks) run_blocks=true; shift;;  # Get position sof synteny blocks between accessions
        -seq)    run_seq=true; shift;;  # Get consencuc seqeunce
        -aln)    run_aln=true; shift;;  # Produce fasta file with the pangenome alignment
        -snp)    run_snp=true; shift;;  # Get VSF file with SNPs
        # -sv)     run_sv=true; shift;;

        -sv_call)   run_sv_call=true; shift;;               # SV calling from the alignment
        -sv_sim) run_sv_sim=true; set_file="$2"; shift 2;;   # File to compare SVs againts set of seequences
        -sv_graph)  run_sv_graph=true; shift;;              # Construction of a graph on SVs
        -sim) similarity_value=$2; shift 2;;                # Similarity value

        -sv_acc) sv_acc=$2; shift ;;  # file with accessions to analyse

        -aln_type) aln_type=$2; shift 2;;
        -path_cons) path_consensus=$2; shift 2;;
        *) print_usage; exit 1;;
    esac
done

cores="${cores:-1}"  # Number of cores
acc_anal="${acc_anal:-NULL}"   # Set of accessions to analyse

pokaz_message "Number of cores: ${cores}"

check_missing_variable "ref_pref"

# check_missing_variable "pref_global"
# pref_global=$(add_symbol_if_missing "$pref_global" "/")

# path_consensus="${path_consensus:-${pref_global}consensus/}"
path_consensus=$(add_symbol_if_missing "$path_consensus" "/")

# path_chromosomes="${path_chromosomes:-${pref_global}chromosomes/}"
path_chromosomes=$(add_symbol_if_missing "$path_chromosomes" "/")

# ----------------------------------------------------------------------------
#                                   MAIN
# ----------------------------------------------------------------------------


# ******************  General work with the alignment   **********************

# -------------------------------------------------
if [ "$run_blocks" = true ]; then

    Rscript analys/analys_01_blocks.R --path.cons ${path_consensus} --ref.pref  ${ref_pref} --cores ${cores} --aln.type ${aln_type}
fi

if [ "$run_seq" = true ]; then

    Rscript analys/analys_02_seq_cons.R --path.cons ${path_consensus} --ref.pref  ${ref_pref} --path.chromosomes ${path_chromosomes}  --aln.type ${aln_type} --cores ${cores}
fi

if [ "$run_aln" = true ]; then

    Rscript analys/analys_03_seq_aln.R --path.cons ${path_consensus} --ref.pref  ${ref_pref} --path.chromosomes ${path_chromosomes} --aln.type ${aln_type} --cores ${cores}
fi


if [ "$run_snp" = true ]; then

    Rscript analys/analys_04_snp.R --path.cons ${path_consensus} --ref.pref  ${ref_pref} --path.chromosomes ${path_chromosomes}  --aln.type ${aln_type} --cores ${cores}
fi


# ********************************   SVs   ***********************************



# -------------------------------------------------
# Sv calling
if [ "$run_sv_call" = true ]; then

    # Philosophy: GFF does not make any sense without a pangenome consensus fasta. 
    # So, sonsensus should be run before GFF
    # Therefore, sequences of seSVs could also be produced together with GFFs.

    Rscript analys/sv_01_calling.R --path.cons ${path_consensus} --ref.pref  ${ref_pref} --aln.type ${aln_type}
fi



# -------------------------------------------------
# Compare SVs with TEs
if [ "$run_sv_sim" = true ]; then
    check_missing_variable "set_file"

    if [ -z "${similarity_value}" ]; then
        pokaz_message "Simirarity value is 85% (default)"
        similarity_value=85
    fi

    # Check if BLAST database exists
    # if [ ! -f "${set_file}.nhr" ]; then
        makeblastdb -in "$set_file" -dbtype nucl > /dev/null
    # fi
    
    file_sv_big=${path_consensus}sv/seq_sv_big.fasta
    file_sv_big_on_set=${file_sv_big%.fasta}_on_set_blast.txt

    # if [ ! -f "${file_sv_big_on_set}" ]; then
        blastn -db "${set_file}" -query "${file_sv_big}" -out "${file_sv_big_on_set}" \
           -outfmt "7 qseqid qstart qend sstart send pident length sseqid" \
           -perc_identity "${similarity_value}"
    # fi

    file_sv_big_on_set_cover=${file_sv_big%.fasta}_on_set_cover.rds
    Rscript sim/sim_in_seqs.R --in_file ${file_sv_big} --db_file ${set_file} --res ${file_sv_big_on_set} \
            --out ${file_sv_big_on_set_cover} --sim ${similarity_value} --use_strand F

    rm "${set_file}.nin" "${set_file}.nhr" "${set_file}.nsq"
fi

# -------------------------------------------------
# SV on SVs
if [ "$sv_graph" = true ]; then

    if [ -z "${similarity_value}" ]; then
        pokaz_message "Simirarity value is 85% (default)"
        similarity_value=85
    fi

    file_sv_big=${path_consensus}sv/seq_sv_big.fasta
    file_sv_big_on_sv=${file_sv_big%.fasta}_on_sv_blast.txt

    # Check if BLAST database exists
    makeblastdb -in "$file_sv_big" -dbtype nucl > /dev/null
    
    # if [ ! -f "${file_sv_big_on_sv}" ]; then
        blastn -db ${file_sv_big} -query ${file_sv_big} -out ${file_sv_big_on_sv} \
           -outfmt "7 qseqid qstart qend sstart send pident length sseqid" \
           -perc_identity ${similarity_value} 
    # fi

    file_sv_big_on_sv_cover=${file_sv_big%.fasta}_on_sv_cover.rds
    Rscript sim/sim_in_seqs.R --in_file ${file_sv_big} --db_file ${file_sv_big} --res ${file_sv_big_on_sv} \
            --out ${file_sv_big_on_sv_cover} --sim ${similarity_value} --use_strand T

    rm "$file_sv_big".nin
    rm "$file_sv_big".nhr
    rm "$file_sv_big".nsq

fi
