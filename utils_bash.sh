# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

# Function to add a symbol to the end of a string if it's missing
add_symbol_if_missing() {
    local input_string="$1"  # Receive the string as an argument
    local symbol="$2"       # Receive the symbol to add

    # Check if the string and symbol are not empty
    if [ -n "$input_string" ] && [ -n "$symbol" ]; then
        # Check if the last character of the string is not the symbol to add
        if [ "${input_string: -1}" != "$symbol" ]; then
            input_string="$input_string$symbol"
        fi
    fi

    echo "$input_string"
}


# Function to remove a file if it exists
remove_file_if_exists() {
    local file="$1"
    if [ -f "${file}" ]; then
        rm "${file}"
    fi
}

# Function to check if a variable is set
check_missing_variable() {
    local var_name="$1"  # Name of the variable to check

    # Using indirect variable reference to check if the variable is set
    if [ -z "${!var_name}" ]; then
        echo "Error: Variable '$var_name' is not set."
        exit 1
    fi
}