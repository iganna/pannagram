# ----------------------------------------------------------------------------
#           STEP DONE
# ----------------------------------------------------------------------------

with_level 1 pokaz_message "Step is done."

((step_num = step_num + 1))

if [ "${step_num}" -gt "${step_end}" ]; then
    exit 0
fi
