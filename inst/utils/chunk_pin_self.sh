# ----------------------------------------------------------------------------
#            PIN THE RUNNING SCRIPT
# ----------------------------------------------------------------------------
#
# Bash reads a script incrementally, not into memory: it keeps a byte offset in
# the file and comes back for the next command. Editing the file while a job is
# running shifts those offsets, the interpreter resumes mid-line and the job dies
# with a syntax error that `bash -n` cannot reproduce, hours in.
#
# The launchers in $CONDA_PREFIX/bin are symlinks straight into the working tree
# (see pannagram_checks.sh), which is deliberate - a developer edit takes effect
# without a reinstall - but it means an editor save, a `git pull` or a reinstall
# during a multi-hour run rewrites the very file that run is reading.
#
# So the script hands over to a private copy of itself and unlinks it right away:
# the copy is reachable only through the open file descriptor of the interpreter,
# so no later edit of the original can reach a run that has already started.
# Steps started after the handover still source the freshly installed helpers -
# only the entry script is pinned.
#
# Set PANNAGRAM_NO_PIN=T to skip this (tracing, profiling, `bash -x` wrappers).

if [ -z "${PANNAGRAM_PINNED}" ] && [ "${PANNAGRAM_NO_PIN}" != "T" ]; then

    pannagram_self="${BASH_SOURCE[1]:-$0}"
    pannagram_self=$(readlink -f "${pannagram_self}" 2>/dev/null || echo "${pannagram_self}")

    # Keep the tool name in the path: it is what shows up in the shell's error
    # messages and in `ps` from here on.
    pannagram_pin_dir="${TMPDIR:-/tmp}/.pannagram_run_$$"
    pannagram_pin="${pannagram_pin_dir}/$(basename "${pannagram_self}")"

    if [ -r "${pannagram_self}" ] && mkdir -p "${pannagram_pin_dir}" 2>/dev/null \
       && cp "${pannagram_self}" "${pannagram_pin}" 2>/dev/null; then

        chmod u+rx "${pannagram_pin}" 2>/dev/null || true
        export PANNAGRAM_PINNED="${pannagram_pin}"
        # Skip the second `Rscript -e system.file()` in the pinned copy.
        export PANNAGRAM_PATH="${INSTALLED_PATH}"
        exec bash "${pannagram_pin}" "$@"
    fi

    # Copying failed (read-only or full $TMPDIR): run from the original rather
    # than refuse to start, and say why the protection is not there.
    echo "Warning: could not pin ${pannagram_self} in ${pannagram_pin_dir};" >&2
    echo "         do not edit or reinstall the script while this job runs." >&2
    rm -rf "${pannagram_pin_dir}" 2>/dev/null || true
    unset pannagram_self pannagram_pin pannagram_pin_dir

elif [ -n "${PANNAGRAM_PINNED}" ]; then

    # This is the pinned copy. Drop its directory entry: the open descriptor keeps
    # it readable until the run ends, and nothing is left behind afterwards.
    rm -f "${PANNAGRAM_PINNED}" 2>/dev/null || true
    rmdir "$(dirname "${PANNAGRAM_PINNED}")" 2>/dev/null || true
    unset PANNAGRAM_PINNED

fi
