# Get synteny blocks

with_level 1 pokaz_stage "Step ${step_num}. Get synteny blocks."

# Logs
step_name="step${step_num}_analys_01_blocks"
step_file="${path_log}${step_name}_done"
path_log_step="${path_log}${step_name}/"
mkdir -p ${path_log_step}

# Start
if [ "${step_num}" -ge "${step_start}" ] || [ ! -f ${step_file} ]; then

    Rscript $INSTALLED_PATH/analys/analys_01_blocks.R \
            --path.cons ${path_cons} \
            --cores ${cores}
    # Done
    touch "${step_file}"
fi

source $INSTALLED_PATH/utils/chunk_step_done.sh

with_level 1 pokaz_message "* The pipeline is done."