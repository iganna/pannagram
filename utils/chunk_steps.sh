# Path with steps
path_flags="${path_out}.flags/"
path_flags=$(add_symbol_if_missing "$path_flags" "/")
if [ ! -d "$path_flags" ]; then
    mkdir -p "$path_flags"
fi


if [ -z "${start_step}" ]; then
    max_step_file=$(ls ${path_flags}step*_done 2> /dev/null | sort -V | tail -n 1)

    if [ -z "$max_step_file" ]; then
        start_step=1
    else
        start_step=${max_step_file##*step}
        start_step=${start_step%_done}
        start_step=$((start_step + 1))
    fi

fi

# Looping through and deleting files of stages, which are greater or equal to the current one
for file_step in "${path_flags}"step*_done*; do
    # echo "Processing file: $file_step"
    if [ -f "$file_step" ]; then
        # Extracting step number from the file name
        step_tmp=$(echo "$file_step" | sed -e 's/.*step\([0-9]*\)_done.*/\1/')

        # Check if step number is greater or equal to start_step
        if [ "$step_tmp" -ge "$start_step" ]; then
            # echo "remove ${file_step}"
            rm -f "$file_step"
        fi
    fi
done
