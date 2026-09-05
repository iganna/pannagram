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


# Absolute path with '..' and symlinks resolved, so that the same folder given
# in different ways is recognised as the same
canonical_path() {
    local path="$1"
    local real
    real=$(readlink -f "$path" 2>/dev/null)
    [ -n "$real" ] || real="$path"
    add_symbol_if_missing "$real" "/"
}


# Remove a file if it exists
remove_file_if_exists() {
    local file="$1"
    if [ -f "${file}" ]; then
        rm "${file}"
    fi
}


# Remove everything inside a directory, subdirectories included, but keep the
# directory itself. `rm -f "${dir}"*` cannot be used for this: it exits with 1
# on every subdirectory it meets, and `set -e` then kills the whole run.
clean_dir_content() {
    local dir="${1%/}"
    [ -n "${dir}" ] || return 0
    [ -d "${dir}" ] || return 0
    find "${dir}" -mindepth 1 -maxdepth 1 -exec rm -rf -- {} +
}


# Remove files matching a name pattern from a directory, non-recursively.
# `find -exec rm` and not a shell glob: these folders hold tens of thousands
# of files and the expanded glob would overflow the command line.
clean_dir_files() {
    local dir="${1%/}"
    local pattern="$2"
    [ -n "${dir}" ] && [ -n "${pattern}" ] || return 0
    [ -d "${dir}" ] || return 0
    find "${dir}" -maxdepth 1 -type f -name "${pattern}" -exec rm -f -- {} +
}


# Remove a directory if nothing is left in it
remove_dir_if_empty() {
    local dir="${1%/}"
    [ -n "${dir}" ] || return 0
    [ -d "${dir}" ] || return 0
    rmdir "${dir}" 2>/dev/null || true
}


# Remove a directory together with its content
remove_dir_if_exists() {
    local dir="${1%/}"
    [ -n "${dir}" ] || return 0
    [ -d "${dir}" ] || return 0
    rm -rf -- "${dir}"
}


# ----------------------------------------------------------------------------
#           STEP CHECKPOINTS
# ----------------------------------------------------------------------------
#
# Checkpoints are keyed on the step ID - the script the step runs - and NOT on
# the step number: `comb_07_long_done` and `comb_07_long/`, not
# `step16_comb_07_long_done` and `step16_comb_07_long/`. Inserting or removing a
# stage renumbers every step after it, and numbered markers would then all stop
# matching, so a finished project would silently recompute from that point on.
# The number is display-only and is written INSIDE the marker file, so that
# `-s <n>` still knows which markers to drop.

# Set step_file / path_log_step for the step with the given id.
set_step_paths() {
    local id="$1"
    step_file="${path_log}${id}_done"
    path_log_step="${path_log}${id}/"
    mkdir -p "${path_log_step}"
}


# Write the "step done" marker, recording the current step number inside it.
mark_step_done() {
    echo "${step_num}" > "${step_file}"
}


# The step number recorded in a marker file; empty if it cannot be determined.
step_marker_number() {
    local file="$1"
    local n
    n=$(head -n 1 "${file}" 2>/dev/null | tr -d '[:space:]')
    case "${n}" in
        ''|*[!0-9]*)
            # Marker of a project written before the number moved inside the file.
            n=$(basename "${file}" | sed -n 's/^step\([0-9]\{1,\}\)_.*/\1/p')
            ;;
    esac
    printf '%s' "${n}"
}


# Rename one legacy checkpoint marker to its number-free name, keeping the step
# number as the file content. Never overwrites an already migrated marker.
move_legacy_marker() {
    local dir="$1"
    local old="$2"
    local new="$3"
    local num="$4"

    [ -f "${dir}/${old}_done" ] || return 0
    if [ -f "${dir}/${new}_done" ]; then
        rm -f "${dir}/${old}_done"
        return 0
    fi
    echo "${num}" > "${dir}/${new}_done"
    rm -f "${dir}/${old}_done"
    with_level 2 pokaz_message "checkpoint ${old} -> ${new}"
}


# Same for the per-step log folder, which is also the checkpoint ledger of the
# workers. Sibling folders sharing the prefix (e.g. <step>_merge) come along.
move_legacy_dir() {
    local dir="$1"
    local old="$2"
    local new="$3"

    [ -d "${dir}/${old}" ] || return 0
    if [ -d "${dir}/${new}" ]; then
        return 0
    fi
    mv "${dir}/${old}" "${dir}/${new}"
}


# Bring the checkpoints of a project created by an older version over to the
# number-free layout. Without this every existing project would look unfinished
# and recompute from the first step.
migrate_legacy_step_markers() {
    local dir="${1%/}"
    [ -d "${dir}" ] || return 0

    local had_nullglob="F"
    if shopt -q nullglob; then
        had_nullglob="T"
    fi
    shopt -s nullglob

    local f name num id base round
    local -a legacy sorted

    # comb_11/12/13 run twice (extra long fragments, rounds 1 and 2) and used to
    # share an id, told apart only by the step number. Map them by pipeline
    # order: the lower step number is round 1.
    for base in comb_11 comb_12 comb_13; do
        legacy=()
        for f in "${dir}"/step[0-9]*_"${base}"_done; do legacy+=("${f}"); done
        [ ${#legacy[@]} -gt 0 ] || continue

        sorted=()
        while IFS= read -r f; do sorted+=("${f}"); done < <(printf '%s\n' "${legacy[@]}" | sort -V)

        round=0
        for f in "${sorted[@]}"; do
            round=$(( round + 1 ))
            name=$(basename "${f}")
            name=${name%_done}
            num=${name#step}
            num=${num%%_*}
            move_legacy_marker "${dir}" "${name}" "${base}_extra${round}" "${num}"
            move_legacy_dir    "${dir}" "${name}" "${base}_extra${round}"
        done
    done

    # Everything else: step<N>_<id> -> <id>, markers first, then the log folders
    # (a folder has no marker of its own, e.g. <step>_merge).
    for f in "${dir}"/step[0-9]*_done; do
        name=$(basename "${f}")
        name=${name%_done}
        num=${name#step}
        num=${num%%_*}
        id=${name#step${num}_}
        [ -n "${id}" ] && [ "${id}" != "${name}" ] || continue
        move_legacy_marker "${dir}" "${name}" "${id}" "${num}"
    done

    for f in "${dir}"/step[0-9]*_*; do
        [ -d "${f}" ] || continue
        name=$(basename "${f}")
        # A round-1/round-2 folder without its marker cannot be told apart, so it
        # is left as it is rather than migrated under a misleading id.
        case "${name}" in
            *_comb_11|*_comb_12|*_comb_13) continue ;;
        esac
        num=${name#step}
        num=${num%%_*}
        id=${name#step${num}_}
        [ -n "${id}" ] && [ "${id}" != "${name}" ] || continue
        move_legacy_dir "${dir}" "${name}" "${id}"
    done

    [ "${had_nullglob}" == "T" ] || shopt -u nullglob
}


# Verify that a list file exists and that every path listed in it is present.
# Guard against the case where a previous run's checkpoints survived but its
# outputs did not: the per-step checkpoints then report "already done", the
# pipeline marches through several steps in seconds and only dies much later,
# at the first step that actually opens the files. Fail here instead, loudly.
check_listed_files() {
    local list_file="$1"
    local what="${2:-input}"

    if [ ! -f "${list_file}" ]; then
        pokaz_error "Missing ${what} list: ${list_file}"
        exit 1
    fi

    local missing=0
    local f
    while IFS= read -r f || [ -n "${f}" ]; do
        [ -n "${f}" ] || continue
        if [ ! -f "${f}" ]; then
            missing=$(( missing + 1 ))
            if [ "${missing}" -le 5 ]; then
                pokaz_error "Missing ${what} file: ${f}"
            fi
        fi
    done < "${list_file}"

    if [ "${missing}" -ne 0 ]; then
        pokaz_error "${missing} file(s) listed in ${list_file} do not exist."
        pokaz_error "Outputs of an earlier step are gone while its checkpoints survived."
        pokaz_error "Rerun that step with -cleanup, which now removes both together."
        exit 1
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
  
  echo -n "┌"
  printf -- '─%.0s' $(seq 1 $((len + 2)))
  echo "┐"

  echo "│ $message │"

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

joinpath() {
    local joined=""
    local arg
    local first=1

    for arg in "$@"; do
        if [[ $first -eq 1 ]]; then
            joined="${arg%/}"
            first=0
        else
            if [[ "$arg" == /* ]]; then
                pokaz_error "Error: argument '$arg' should not start with a '/'" >&2
                exit 1
            fi
            joined="${joined%/}/${arg#/}"
        fi
    done

    echo "$joined"
}