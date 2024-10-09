# ----------------------------------------------------------------------------
#           STEP DONE
# ----------------------------------------------------------------------------

# If only one step
if   [ "$one_step" = "T" ]; then 
    exit 0
fi
    
with_level 1 pokaz_message "Step is done."

((step_num = step_num + 1))