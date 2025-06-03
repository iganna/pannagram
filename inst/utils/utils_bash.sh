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
    color_code="38;5;37"
    printf "\e[${color_code}m* %s\e[0m\n" "$text"
}


# Display an attention message
pokaz_attention() {
    local text="$1"
    local color_code="38;5;203"
    printf "\e[${color_code}m  %s\e[0m\n" "$text"
}

# Display an error message in red
pokaz_error() {
    local text="$1"
    local color_code="31"  # ANSI color code for red
    printf "\e[${color_code}m%s\e[0m\n" "$text"
}

# Display a stage message
pokaz_message() {
    local text="$1"
    local color_code="38;5;158"  # Very light blue color code
    printf "\e[${color_code}m  %s\e[0m\n" "$text"
    # echo "  ${text}"
}


# Display a the command and evaluate
show_run() {
    local cmd="$@"    
    printf "\033[1;35mRunning command: $cmd\033[0m\n"
    eval "$cmd"
}


# Display the help message
pokaz_help() {
    pokaz_message "< Welcome to Hellp >"
}

help_in_box() {
    print_fancy_frame "Get help by running: ${0##*/} -h"
}

# Logging messages either to the console or to a specified file based on the given log level.
# Logging into files - always
# Logging into the console - based on the level
# log_message() {
#     local log_level_command=$1
#     local log_level=$2
#     local file_log=$3
#     shift 3

#     local pokaz_command=$1
#     shift
#     local message="$*"

#     # echo '----'
#     # echo ${log_level_command} 
#     # echo ${log_level} 
#     # echo ${file_log} 
#     # echo ${pokaz_command} 
#     # echo ${message}
#     # echo '==='

#     # Print to the console
#     if [ "${log_level_command}" -le "${log_level}" ]; then
#         $pokaz_command "$message"
#     fi

#     $pokaz_command "$message" | sed 's/\x1b\[[0-9;]*m//g' >> "$file_log"

# }

# Echo a fancy frame around the messase
print_fancy_frame() {
  local message="$1"
  local len=${#message}
  
  # Верхняя граница
  echo -n "┌"
  printf -- '─%.0s' $(seq 1 $((len + 2)))
  echo "┐"

  # Текст
  echo "│ $message │"

  # Нижняя граница
  echo -n "└"
  printf -- '─%.0s' $(seq 1 $((len + 2)))
  echo "┘"
}

require_arg() {
    if [ $# -lt 2 ]; then
        echo "Error: $1 requires an argument" >&2
        exit 1
    fi
}
