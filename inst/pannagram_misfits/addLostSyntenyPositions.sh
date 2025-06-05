Add synteny positions which were lost

with_level 1 pokaz_stage "Step ${step_num}. Add synteny positions which were lost."

# Paths
path_extra="${path_inter}extra_regions/"
if [ ! -d "$path_extra" ]; then
    mkdir -p "$path_extra"
fi

# Logs
step_name="step${step_num}_comb_08"
step_file="${path_log}${step_name}_done"
path_log_step="${path_log}${step_name}/"
mkdir -p ${path_log_step}

# Start
if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

    # ---- Clean up the output folders ----
    if   [ "$clean" == "T" ]; then 
        touch ${path_cons}aln_fake_h5
        touch ${path_log_step}fake.log

        rm -f ${path_cons}aln*h5
        rm -f ${path_log_step}*
    fi  

    Rscript $INSTALLED_PATH/pangen/comb_08_insert_back.R  \
            --cores ${cores} \
            --path.cons ${path_cons} \
            --path.log ${path_log_step} \
            --log.level ${log_level}

    # Done
    touch "${step_file}"
fi

source $INSTALLED_PATH/utils/chunk_step_done.sh