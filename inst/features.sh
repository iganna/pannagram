#!/bin/bash

INSTALLED_PATH=$(Rscript -e "cat(system.file(package = 'pannagram'))")

source "$INSTALLED_PATH/utils/chunk_error_control.sh"
source "$INSTALLED_PATH/utils/utils_bash.sh"
source "$INSTALLED_PATH/utils/help_features.sh"
source "$INSTALLED_PATH/utils/argparse_features.sh" "$@"
source "$INSTALLED_PATH/utils/chunk_paths.sh" # requires path_project variable


# General alignment processing
if [ "$run_blocks" = true ]; then # -blocks
    pokaz_stage "Get blocks."

    check_dir "$path_inter_msa"    || exit 1
    check_dir "$path_features_msa" || exit 1

    mkdir -p ${path_plots_synteny}

    Rscript $INSTALLED_PATH/analys/analys_01_blocks3.R \
        --path.inter.msa ${path_inter_msa} \
        --path.features.msa ${path_features_msa} \
        --path.figures ${path_plots_synteny} \
        --cores ${cores} \
        --ref  ${ref_pref} \
        --aln.type ${aln_type}
fi

if [ "$run_seq" = true ]; then # -seq
    pokaz_stage "Get consensus sequences."

    check_dir "$path_features_msa" || exit 1
    check_dir "$path_chrom"        || exit 1

    mkdir -p $path_seq
    Rscript $INSTALLED_PATH/analys/analys_02_seq_cons.R \
        --path.features.msa ${path_features_msa} \
        --ref.pref  ${ref_pref} \
        --path.chr ${path_chrom} \
        --aln.type ${aln_type} \
        --cores ${cores}
fi

if [ "$run_aln" = true ]; then # -aln
    pokaz_stage "Get sequences of alignment."
    Rscript $INSTALLED_PATH/analys/analys_03_seq_aln.R \
        --path.cons ${path_consensus} \
        --ref.pref  ${ref_pref} \
        --path.chromosomes ${path_chromosomes} \
        --aln.type ${aln_type} \
        --cores ${cores}
fi


if [ "$run_snp" = true ]; then # -snp
    pokaz_stage "Get SNPs."

    check_dir "$path_features_msa" || exit 1
    check_dir "$path_seq"          || exit 1

    mkdir -p $path_snp

    Rscript $INSTALLED_PATH/analys/analys_04_snp.R \
        --path.cons ${path_consensus} \
        --ref.pref  ${ref_pref} \
        --path.chr ${path_chromosomes} \
        --aln.type ${aln_type} \
        --cores ${cores}

    # ---------------
    # Pi diversity
    if [ "$run_snp_pi" = true ]; then # -snp_pi

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
            echo "VCF-tools.."
            base_name=$(basename "$vcf_file" .vcf)
            output_file="${path_snp}${base_name}_output"
            vcftools --vcf "$vcf_file" --site-pi --out "$output_file" > /dev/null
            vcftools --vcf "${vcf_file}" --extract-FORMAT-info ID --out "$output_file"

            # Rscript $INSTALLED_PATH/analys/analys_04_snp_plot.R \
            #     --path.figures ${path_plots} \
            #     --file.pi "${output_file}.sites.pi"

            echo "Plink.."
            plink --vcf "${vcf_file}" --distance  --out "${vcf_file}.dist" --allow-extra-chr

            # Rscript $INSTALLED_PATH/analys/analys_04_snp_dist.R \
            #    --path.figures ${path_plots} \
            #    --file.pi "${vcf_file}"

        done
    fi

    # rm ${path_snp}*log

fi


# ********************************   SVs   ***********************************

# -------------------------------------------------
# Sv calling
if [ "$run_sv_call" = true ]; then # -sv_call|-sv
    pokaz_stage "SV-calling"
    # Philosophy: GFF does not make any sense without a pangenome consensus fasta. 
    # So, consensus should be run before GFF
    # Therefore, sequences of seSVs could also be produced together with GFFs.

    check_dir "$path_features_msa" || exit 1
    check_dir "$path_seq"          || exit 1

    mkdir -p $path_sv
    mkdir -p $path_gff

    Rscript $INSTALLED_PATH/analys/sv_01_calling.R \
        --path.cons ${path_consensus} \
        --ref.pref  ${ref_pref} \
        --aln.type ${aln_type} \
        --acc.anal ${acc_anal}

    # path_plots="${path_consensus}plot_svs/"
    # mkdir -p ${path_plots}

    pokaz_stage "Plot SV stat"
    Rscript $INSTALLED_PATH/analys/sv_02_plot_stat.R \
        --path.cons ${path_consensus} 

fi



# -------------------------------------------------
# Compare SVs with TEs
if [ "$run_sv_sim" = true ]; then # -sv_sim
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
           -perc_identity "${similarity_value}" -num_threads "${cores}"
    # fi

    file_sv_big_on_set_cover=${file_sv_big%.fasta}_on_set_cover.rds
    Rscript $INSTALLED_PATH/sim/sim_in_seqs.R --in_file ${file_sv_big} --db_file ${set_file} --res ${file_sv_big_on_set} \
            --out ${file_sv_big_on_set_cover} --sim ${similarity_value} --use_strand F

    rm "${set_file}.nin" "${set_file}.nhr" "${set_file}.nsq"
fi

# -------------------------------------------------
# SV on SVs
if [ "$run_sv_graph" = true ]; then # -sv_graph

    pokaz_stage "Graph on SVs"
    if [ -z "${similarity_value}" ]; then
        pokaz_message "Simirarity value is 85% (default)"
        similarity_value=85
    fi

    check_dir "$path_features_msa" || exit 1
    check_dir "$path_sv"           || exit 1
    check_dir "$path_plots_sv"     || exit 1

    file_sv_big=${path_sv}seq_sv_big.fasta
    file_sv_big_on_sv=${file_sv_big%.fasta}_on_sv_blast.txt

    # Check if BLAST database exists
    makeblastdb -in "$file_sv_big" -dbtype nucl > /dev/null
    
    # if [ ! -f "${file_sv_big_on_sv}" ]; then
        blastn -db ${file_sv_big} -query ${file_sv_big} -out ${file_sv_big_on_sv} \
           -outfmt "6 qseqid qstart qend sstart send pident length sseqid" \
           -perc_identity ${similarity_value} -num_threads "${cores}"
    # fi
    pokaz_message "Blast is done."

    file_sv_big_on_sv_cover=${file_sv_big%.fasta}_on_sv_cover.rds
    Rscript $INSTALLED_PATH/sim/sim_in_seqs.R --in_file ${file_sv_big} --db_file ${file_sv_big} --res ${file_sv_big_on_sv} \
            --out ${file_sv_big_on_sv_cover} --sim ${similarity_value} --use_strand T

    rm "$file_sv_big".nin
    rm "$file_sv_big".nhr
    rm "$file_sv_big".nsq


    pokaz_stage "Plotting SV-Graph..."
    Rscript $INSTALLED_PATH/analys/sv_03_plot_graph.R \
        --path.cons ${path_consensus} 

    Rscript $INSTALLED_PATH/analys/sv_04_orfs_in_graph.R \
        --path.cons ${path_consensus} 

    if [ "$run_sv_sim_prot" = true ]; then # -sv_sim_prot

        if [ -f "${path_consensus}sv/sv_in_graph_orfs.fasta" ]; then
            pokaz_stage "BLAST on proteins..."

            # makeblastdb -in ${set_file_prot} -dbtype prot
            blastp -db "${set_file_prot}" \
                   -query "${path_consensus}sv/sv_in_graph_orfs.fasta" \
                   -out "${path_consensus}sv/blast_sv_orfs_on_set.txt" \
                   -outfmt "7 qseqid qstart qend sstart send pident length sseqid" \
                   -num_threads "${cores}"
        else
            pokaz_error "File with ORFs does not exist, BLAST against proteins was not performed."
        fi

    fi

fi


# -------------------------------------------------
# Annotation groups
if [ "$run_annogroup" = true ]; then # -annogroup

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
