# ----------------------------------------------------------------------------
#           LOGS
# ----------------------------------------------------------------------------

log_level=${log_level:-1}  # Set the default value to 'steps'

if ! [[ "$log_level" =~ ^[0-3]$ ]]; then
    pokaz_error "Error: log_level must be a number between 0 and 3."
    exit 1
fi

# Hidden path with logs
path_logs="${path_out}logs/"
make_dir ${path_logs}

# Path for logs of this script
path_log_ref="${path_logs}pangen_ref_${ref_name}/"
make_dir ${path_log_ref}

# File with steps logs
file_log_ref="${path_log_ref}steps.log"
