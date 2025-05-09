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
    cat << EOF
Usage: ${0##*/}  -path_msa PATH_MSA  -path_chr PATH_CHR 
                [-ref REF]
                [-h] [-cores NUM_CORES]  
                [-blocks] [-seq] [-aln] [-snp] 
                [-aln_type ALN_TYPE] [-path_cons PATH_CONS]


This script manages various genomic analyses and alignments.

Options:
    -h, --help                  Display this help message and exit.
    -cores NUM_CORES            Specify the number of cores for parallel processing (default is 1).

    -path_msa PATH_MSA          Specify the global prefix for multiple sequence alignment. The same as -path_out in pangen.sh
    -path_chr PATH_CHR          Specify the path to chromosome files.

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
run_snp_pi=false

# analysis of SVs
run_sv_call=false
run_sv_sim=false
run_sv_graph=false
run_annogroup=false
run_sv_sim_prot=false

# Parse command line arguments
while [ $# -gt 0 ]; do
    case $1 in
        -h|-help)        print_usage;            exit 0;;
        -cores)          cores=$2;               shift 2;;

        -path_msa)       path_consensus=$2;      shift 2;;
        -ref)            ref_pref=$2;            shift 2;;
        -path_chr)       path_chromosomes=$2;    shift 2;;
        
        -blocks)         run_blocks=true;        shift;;           # Get position of synteny blocks between accessions
    
        -seq)            run_seq=true;           shift;;           # Get consensus sequence
        -aln)            run_aln=true;           shift;;           # Produce fasta file with the pangenome alignment
        -snp)            run_snp=true;           shift;;           # Get VCF file with SNPs
        -snp_pi)         run_snp_pi=true;        shift;;           # Run pi-diversity with VCF-tools
        # -sv)           run_sv=true;            shift;;

        -sv_call|-sv)    run_sv_call=true;       shift;;           # SV calling from the alignment
        -sv_sim)         set_file="$2";                            # File to compare SVs against set of sequences
                         run_sv_sim=true;        shift 2;;         
        -sv_sim_prot)    set_file_prot="$2";                            # File to compare SVs against set of sequences
                         run_sv_sim_prot=true;        shift 2;;         
        -sv_graph)       run_sv_graph=true;      shift;;           # Construction of a graph on SVs
        -sim)            similarity_value=$2;    shift 2;;         # Similarity value

        -sv_acc)         acc_anal=$2;            shift 2;;         # File with accessions to analyse

        -annogroup)      path_annot="$2"                           # Path with annotation
                         run_annogroup=true;     shift 2;;

        -aln_type)       aln_type=$2;            shift 2;;
        -path_cons)      path_consensus=$2;      shift 2;;
        *)               print_usage; echo "Wrong parameter ${1}";            exit 1;;
    esac
done


# Check if path_chromosomes is empty while any of run_seq, run_aln, or run_snp are set to true
if [ -z "$path_chromosomes" ] && ([ "$run_seq" = true ] || [ "$run_aln" = true ] || [ "$run_snp" = true ]); then
    pokaz_error "Error: -path_chr must be specified when any of -seq, -aln, or -snp options are used."
    exit 1
fi

cores="${cores:-1}"  # Number of cores
acc_anal="${acc_anal:-NULL}"   # Set of accessions to analyse
ref_pref="${ref_pref:-NULL}"   # Set of accessions to analyse

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
    pokaz_stage "Get blocks."

    path_plots="${path_consensus}plot_synteny/"
    mkdir -p ${path_plots}

    Rscript $INSTALLED_PATH/analys/analys_01_blocks3.R \
        --path.cons ${path_consensus} \
        --cores ${cores} \
        --ref  ${ref_pref} \
        --aln.type ${aln_type} \
        --path.figures ${path_plots}
fi

if [ "$run_seq" = true ]; then
    pokaz_stage "Get consensus sequences."
    Rscript $INSTALLED_PATH/analys/analys_02_seq_cons.R \
        --path.cons ${path_consensus} \
        --ref.pref  ${ref_pref} \
        --path.chr ${path_chromosomes} \
        --aln.type ${aln_type} \
        --cores ${cores}
fi

if [ "$run_aln" = true ]; then
    pokaz_stage "Get sequences of alignment."
    Rscript $INSTALLED_PATH/analys/analys_03_seq_aln.R \
        --path.cons ${path_consensus} \
        --ref.pref  ${ref_pref} \
        --path.chromosomes ${path_chromosomes} \
        --aln.type ${aln_type} \
        --cores ${cores}
fi


if [ "$run_snp" = true ]; then
    pokaz_stage "Get SNPs."
    # Rscript $INSTALLED_PATH/analys/analys_04_snp.R \
    #     --path.cons ${path_consensus} \
    #     --ref.pref  ${ref_pref} \
    #     --path.chr ${path_chromosomes} \
    #     --aln.type ${aln_type} \
    #     --cores ${cores}

    # ---------------
    # Pi divirsity
    if [ "$run_snp_pi" = true ]; then

        pokaz_stage "Pi diversity."
        path_snp="${path_consensus}snps/"
        vcf_files=$(find "$path_snp" -type f -name "*.vcf")
        if [ -z "$vcf_files" ]; then
            echo "VCF-files are not found in $path_snp!"
            exit 1
        fi

        path_plots="${path_snp}plot_snp/"
        mkdir -p ${path_plots}

        # Run VCF-tools
        for vcf_file in $vcf_files; do
            base_name=$(basename "$vcf_file" .vcf)
            output_file="${path_snp}${base_name}_output"
            #vcftools --vcf "$vcf_file" --site-pi --out "$output_file" > /dev/null
            #vcftools --vcf "${vcf_file}" --extract-FORMAT-info ID --out "$output_file"

            Rscript $INSTALLED_PATH/analys/analys_04_snp_plot.R \
                --path.figures ${path_plots} \
                --file.pi "${output_file}.sites.pi"

            #plink --vcf "${vcf_file}" --distance  --out "${vcf_file}.dist" --allow-extra-chr

            #Rscript $INSTALLED_PATH/analys/analys_04_snp_dist.R \
            #    --path.figures ${path_plots} \
            #    --file.pi "${vcf_file}"

        done
    fi

    # rm ${path_snp}*log

fi


# ********************************   SVs   ***********************************

# -------------------------------------------------
# Sv calling
if [ "$run_sv_call" = true ]; then
    pokaz_stage "SV-calling"
    # Philosophy: GFF does not make any sense without a pangenome consensus fasta. 
    # So, consensus should be run before GFF
    # Therefore, sequences of seSVs could also be produced together with GFFs.

    Rscript $INSTALLED_PATH/analys/sv_01_calling.R \
        --path.cons ${path_consensus} \
        --ref.pref  ${ref_pref} \
        --aln.type ${aln_type} \
        --acc.anal ${acc_anal}

    # path_plots="${path_consensus}plot_svs/"
    # mkdir -p ${path_plots}

    Rscript $INSTALLED_PATH/analys/sv_02_plot_stat.R \
        --path.cons ${path_consensus} 

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
           -outfmt "6 qseqid qstart qend sstart send pident length sseqid" \
           -perc_identity "${similarity_value}"
    # fi

    file_sv_big_on_set_cover=${file_sv_big%.fasta}_on_set_cover.rds
    Rscript $INSTALLED_PATH/sim/sim_in_seqs.R --in_file ${file_sv_big} --db_file ${set_file} --res ${file_sv_big_on_set} \
            --out ${file_sv_big_on_set_cover} --sim ${similarity_value} --use_strand F

    rm "${set_file}.nin" "${set_file}.nhr" "${set_file}.nsq"
fi

# -------------------------------------------------
# SV on SVs
if [ "$run_sv_graph" = true ]; then

    # pokaz_stage "Graph on SVs"
    # if [ -z "${similarity_value}" ]; then
    #     pokaz_message "Simirarity value is 85% (default)"
    #     similarity_value=85
    # fi

    # file_sv_big=${path_consensus}sv/seq_sv_big.fasta
    # file_sv_big_on_sv=${file_sv_big%.fasta}_on_sv_blast.txt

    # # Check if BLAST database exists
    # makeblastdb -in "$file_sv_big" -dbtype nucl > /dev/null
    
    # # if [ ! -f "${file_sv_big_on_sv}" ]; then
    #     blastn -db ${file_sv_big} -query ${file_sv_big} -out ${file_sv_big_on_sv} \
    #        -outfmt "6 qseqid qstart qend sstart send pident length sseqid" \
    #        -perc_identity ${similarity_value} 
    # # fi

    # file_sv_big_on_sv_cover=${file_sv_big%.fasta}_on_sv_cover.rds
    # Rscript $INSTALLED_PATH/sim/sim_in_seqs.R --in_file ${file_sv_big} --db_file ${file_sv_big} --res ${file_sv_big_on_sv} \
    #         --out ${file_sv_big_on_sv_cover} --sim ${similarity_value} --use_strand T

    # rm "$file_sv_big".nin
    # rm "$file_sv_big".nhr
    # rm "$file_sv_big".nsq


    pokaz_stage "Plotting SV-Graph..."
    Rscript $INSTALLED_PATH/analys/sv_03_plot_graph.R \
        --path.cons ${path_consensus} 


    Rscript $INSTALLED_PATH/analys/sv_04_orfs_in_graph.R \
        --path.cons ${path_consensus} 

    if [ "$run_sv_sim_prot" = true ]; then
        makeblastdb -in ${set_file_prot} -dbtype prot
        blastp -db ${set_file_prot} -query ${path_consensus}sv/sv_in_graph_orfs.fasta \
                -out blast_sv_orfs_on_set.txt \
                -outfmt "7 qseqid qstart qend sstart send pident length  sseqid" -num_threads ${cores}

    fi

fi


# -------------------------------------------------
# Annotation groups
if [ "$run_annogroup" = true ]; then

    pokaz_stage "Annotation groups"

    path_annot_res=${path_consensus}annotation/
    mkdir -p ${path_annot_res}

    Rscript $INSTALLED_PATH/analys/analys_05_annogroups_easier.R \
            --path.msa ${path_consensus} \
            --path.annot ${path_annot} \
            --path.res ${path_annot_res} \
            --aln.type ${aln_type}

fi

echo "Script completed successfully"
