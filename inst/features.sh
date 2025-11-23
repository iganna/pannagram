#!/bin/bash

INSTALLED_PATH=$(Rscript -e "cat(system.file(package = 'pannagram'))")

source "$INSTALLED_PATH/utils/chunk_error_control.sh"
source "$INSTALLED_PATH/utils/utils_bash.sh"
source "$INSTALLED_PATH/utils/help_features.sh"
source "$INSTALLED_PATH/utils/argparse_features.sh" "$@"
source "$INSTALLED_PATH/utils/chunk_paths.sh" # requires path_project variable

check_dir ${path_project}

# General alignment processing
if [ "$run_blocks" = true ]; then # -blocks
    pokaz_stage "Get synteny blocks."

    # check_dir "$path_inter_msa"    || exit 1
    check_dir "$path_features_msa" || exit 1

    mkdir -p "${path_inter_msa}"
    mkdir -p ${path_plots_synteny}

    Rscript $INSTALLED_PATH/analys/analys_01_blocks3.R \
        --path.inter.msa ${path_inter_msa} \
        --path.features.msa ${path_features_msa} \
        --path.figures ${path_plots_synteny} \
        --cores ${cores} \
        --ref  ${ref_pref} \
        --aln.type ${aln_type} 
        
    pokaz_message "Step -blocks is done!"
fi


if [ "$run_seq" = true ]; then # -seq
    pokaz_stage "Get consensus sequences."

    check_dir "$path_features_msa" || exit 1
    check_dir "$path_chrom"        || exit 1

    mkdir -p $path_seq
    Rscript $INSTALLED_PATH/analys/analys_02_seq_cons.R \
        --path.features.msa ${path_features_msa} \
        --ref  ${ref_pref} \
        --path.chr ${path_chrom} \
        --path.seq ${path_seq} \
        --aln.type ${aln_type} \
        --cores ${cores}
    
    pokaz_message "Step -seq is done!"
fi


# Sequences of alignment
# if [ "$run_aln" = true ]; then # -aln
#     pokaz_stage "Get sequences of alignment."
#     Rscript $INSTALLED_PATH/analys/analys_03_seq_aln.R \
#         --path.features.msa ${path_features_msa} \
#         --ref.pref  ${ref_pref} \
#         --aln.type ${aln_type} \
#         --cores ${cores}
# fi


if [ "$run_snp" = true ]; then # -snp
    pokaz_stage "Get SNPs."

    check_dir "$path_features_msa" || exit 1
    check_dir "$path_seq"          || exit 1

    mkdir -p $path_snp

    Rscript $INSTALLED_PATH/analys/analys_04_snp.R \
        --path.features.msa ${path_features_msa} \
        --path.snp ${path_snp} \
        --path.seq ${path_seq} \
        --ref  ${ref_pref} \
        --aln.type ${aln_type} \
        --cores ${cores}

    pokaz_message "Step -snp is done!"

    # rm ${path_snp}*log
fi


# Pi diversity
if [ "$run_snp_pi" = true ]; then # -snp_pi

    pokaz_stage "Pi diversity."
    vcf_files=$(find "$path_snp" -type f -name "*.vcf")
    if [ -z "$vcf_files" ]; then
        echo "VCF-files are not found in $path_snp!"
        exit 1
    fi

    mkdir -p ${path_plots_snp}

    # Run VCF-tools
    for vcf_file in $vcf_files; do
        echo "VCF-tools.."
        base_name=$(basename "$vcf_file" .vcf)
        output_file="${path_snp}${base_name}_output"
        vcftools --vcf "$vcf_file" --site-pi --out "$output_file" &>/dev/null
        vcftools --vcf "$vcf_file" --extract-FORMAT-info ID --out "$output_file" &>/dev/null

        Rscript $INSTALLED_PATH/analys/analys_04_snp_plot.R \
            --path.figures ${path_plots_snp} \
            --file.pi "${output_file}.sites.pi"

        echo "Plink.."
        plink --vcf "${vcf_file}" --distance --out "${vcf_file}.dist" --allow-extra-chr > /dev/null

        Rscript $INSTALLED_PATH/analys/analys_04_snp_dist.R \
            --path.figures ${path_plots_snp} \
            --file.pi ${vcf_file} \
            --path.snp ${path_snp}


        rm -f "${path_snp}"*FORMAT
        rm -f "${path_snp}"*nosex
        rm -f "${path_snp}"*log
        rm -f "${path_snp}"*dist
        rm -f "${path_snp}"*id
        
    done
    pokaz_message "Step -snp_pi is done!"
fi


# SV calling
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
        --path.features.msa ${path_features_msa} \
        --path.seq ${path_seq} \
        --path.sv ${path_sv} \
        --path.gff ${path_gff} \
        --ref  ${ref_pref} \
        --aln.type ${aln_type} \
        --acc.anal ${acc_anal}
    
    mkdir -p $path_plots_sv

    pokaz_stage "Plot SV stat"
    Rscript $INSTALLED_PATH/analys/sv_02_plot_stat.R \
        --path.features.msa ${path_features_msa} \
        --path.sv ${path_sv} \
        --path.figures ${path_plots_sv}

    pokaz_message "Step -sv is done!"
fi

# ORFs in SVs
if [ "$run_sv_orf" = true ]; then # -sv_orf
    pokaz_stage "Get ORFs from SVs"

    Rscript $INSTALLED_PATH/analys/sv_05_orfs.R \
        --path.sv ${path_sv}
fi

# Compare SVs with TEs
# if [ "$run_annogroup" = true ]; then # -annogroup
#     pokaz_stage "Annotation groups"
#     path_annot_res=${path_consensus}annotation/
#     mkdir -p ${path_annot_res}
#     Rscript $INSTALLED_PATH/analys/analys_05_annogroups_easier.R \
#             --path.msa ${path_consensus} \
#             --path.annot ${path_annot} \
#             --path.res ${path_annot_res} \
#             --aln.type ${aln_type}
# fi


# SV on SVs
if [ "$run_sv_graph" = true ]; then # -sv_graph

    pokaz_stage "Graph on SVs"
    if [ -z "${similarity_value}" ]; then
        similarity_value=85
    fi
    coverage_value=${similarity_value}
    pokaz_message "Similarity: ${similarity_value}. Coverage: ${coverage_value}."

    check_dir "$path_features_msa" || exit 1
    check_dir "$path_sv"           || exit 1
    check_dir "$path_plots_sv"     || exit 1

    file_sv_large=${path_sv}seq_sv_large.fasta
    file_sv_large_on_sv=${file_sv_large%.fasta}_on_sv_blast.txt
    file_sv_large_on_sv_mv="${path_sv}nestedness_sv_large_85_85.txt"

    if [[ ! -f "$file_sv_large_on_sv_mv" ]]; then
        pokaz_message "Run Simsearch.."

        simsearch -in_seq ${file_sv_large} \
                  -on_seq ${file_sv_large} \
                  -sim ${similarity_value} \
                  -cov ${coverage_value} \
                  -out ${path_sv_simsearch} \
                  -cores "${cores}"

        if [ -f "${path_sv_simsearch}seq_sv_large_85_85.txt" ]; then
            mv "${path_sv_simsearch}seq_sv_large_85_85.txt" "${file_sv_large_on_sv_mv}"
            rm -rf ${path_sv_simsearch}
        fi
    else
        pokaz_message "Simsearch was executed before."
    fi

    pokaz_stage "Building SV-Graph for Mobile Element Families..."
    Rscript $INSTALLED_PATH/analys/sv_03_graph_build.R \
        --path.sv ${path_sv} \
        --path.figures ${path_plots_sv} \
        --file.nestedness ${file_sv_large_on_sv_mv} \
        --flag.plot ${plot_families}

    pokaz_stage "Get ORFs from Families..."
    Rscript $INSTALLED_PATH/analys/sv_04_orfs_in_graph.R \
        --path.features.msa ${path_features_msa} \
        --path.sv ${path_sv} 

    pokaz_message "Step -sv_graph is done!"
fi

# BLAST ORFs against the database
if [ "$run_sv_sim_prot" = true ]; then # -sv_sim_prot

    if [ ! -f "${set_file_prot}" ]; then
        pokaz_error "File with proteins does not exist, provide an existing file."
    elif [ -f "${path_sv}/sv_large_orfs.fasta" ]; then

        # Define the output file
        file_sv_large_on_set="${path_sv}sv_large_orfs_on_set.txt"
        if [ -f "$file_sv_large_on_set" ]; then
            pokaz_attention "Warning: you have already run the search for ORFs of SVs against a database."
            timestamp=$(date +%Y%m%d_%H%M%S)
            file_sv_large_on_set="${path_sv}sv_large_orfs_on_set_${timestamp}.txt"
            echo "The result of the new search will be saved in: ${file_sv_large_on_set}"
        else
            echo "The result of the search for ORFs of SVs against the database will be saved in: ${file_sv_large_on_set}."
        fi

        path_simsearch_out="${path_sv}.simsearch/"
        # simsearch -in_seq "${path_sv}seq_sv_large_orfs.fasta" \
        #           -on_seq ${set_file_prot} \
        #           -out ${path_simsearch_out} \
        #           -cores "${cores}" \
        #           -prot \
        #           -sim 10 \
        #           -cov 80 

        set_file_base_prot="$(basename "$set_file_prot")"
        makeblastdb -in "$set_file_prot" -dbtype prot -out "${path_simsearch_out}${set_file_base_prot}" > /dev/null

        # Run blast
        pokaz_stage "BLAST on proteins..."

        blastp -db "${path_simsearch_out}${set_file_base_prot}" \
               -query "${path_sv}seq_sv_large_orfs.fasta" \
               -out ${file_sv_large_on_set} \
               -outfmt "6 qseqid qstart qend sstart send pident length sseqid  qlen slen" \
               -num_threads "${cores}"
    else
        pokaz_error "File with ORFs does not exist, BLAST against proteins was not performed."
    fi

fi

# Similarity with the set of sequences


# Run SV similarity search if requested
if [[ "${run_sv_sim}" == "true" ]]; then  # -sv_sim
    check_missing_variable "set_file"

    
    set_file_base="$(basename "${set_file%.*}")"
    file_sv_large_on_set="${path_sv}/nestedness_sv_large_on_${set_file_base}.txt"

    # If result already exists, write to a timestamped file instead of clobbering
    if [[ -f "$file_sv_large_on_set" ]]; then
        pokaz_attention "Warning: you have already run the search for ORFs of SVs against a database."
        timestamp="$(date +%Y%m%d_%H%M%S)"
        file_sv_large_on_set="${path_sv}/sv_large_on_${set_file_base}_${timestamp}.txt"
        echo "The result of the new search will be saved in: ${file_sv_large_on_set}"
    else
        echo "The result of the search for ORFs of SVs against the database will be saved in: ${file_sv_large_on_set}."
    fi

    # Temporary output directory for simsearch; keep it hidden (leading dot)
    path_simsearch_out="${path_sv}/.simsearch_on_${set_file_base}/"
    mkdir -p "${path_simsearch_out}"

    # --- Run simsearch ---
    simsearch \
        -in_seq  "${path_sv}/seq_sv_large.fasta" \
        -on_seq  "${set_file}" \
        -out     "${path_simsearch_out}" \
        -cores   "${cores}" \
        -sim     85 \
        -cov     85

    # Expected simsearch summary file (adjust if your tool uses a different name)
    expected_out="${path_simsearch_out}${set_file_base}_85_85.txt"

    if [[ -f "${expected_out}" ]]; then
        mv -f "${expected_out}" "${file_sv_large_on_set}"
        rm -rf "${path_simsearch_out}"
        echo "Search finished. Results saved to: ${file_sv_large_on_set}"
    else
        pokaz_attention "simsearch did not produce the expected file: ${expected_out}"
        echo "Keeping intermediate directory for inspection: ${path_simsearch_out}"
        exit 1
    fi
fi

pokaz_message "Script completed successfully"
