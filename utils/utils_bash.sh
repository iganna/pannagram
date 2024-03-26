# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

# Add a symbol to the end of a string if it's missing
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


# Remove a file if it exists
remove_file_if_exists() {
    local file="$1"
    if [ -f "${file}" ]; then
        rm "${file}"
    fi
}


# Check if a variable is set
check_missing_variable() {
    local var_name="$1"  # Name of the variable to check

    # Using indirect variable reference to check if the variable is set
    if [ -z "${!var_name}" ]; then
        echo "Error: Variable '$var_name' is not set."
        exit 1
    fi
}


# Display a stage message
pokaz_stage() {
    local text="$1"
    local color_code="38;2;52;252;252"  
    echo -e "\e[${color_code}m* ${text}\e[0m"
}


# Display an attention message
pokaz_attention() {
    local text="$1"
    local color_code="38;2;52;252;252"  
    echo -e "\e[${color_code}m* ${text}\e[0m"
}


# Display a stage message
pokaz_message() {
    local text="$1"
    local color_code="38;5;195"  # Very light blue color code
    echo -e "\e[${color_code}m  ${text}\e[0m"
}


# Display the help message
pokaz_help() {
    pokaz_message "< Welcome to Hellp >"
}

