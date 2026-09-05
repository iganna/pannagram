# Building Conda package is currently under rework

## Current workflow
1. Create conda env
    ```sh
    conda env create -f pannagram.yml
    ```
2. Activate conda env
    ```sh
    conda activate pannagram
    ```
3. 
    * Run `./user.sh` to build R package normally and create symlinks to `bash` scripts

        or
    * Run `./developer.sh` to build R package in a quick way (no documentation) and create symlinks to `bash` scripts

## Editing the scripts while a job is running

`$CONDA_PREFIX/bin/pannagram`, `features`, `simsearch` and the rest are **symlinks
into this working tree** (`pannagram_checks.sh`), so an edit takes effect without a
reinstall. Bash, however, reads a script incrementally, keeping a byte offset in the
file: rewriting it under a running job shifts those offsets and the interpreter
resumes mid-line, dying with a syntax error that `bash -n` cannot reproduce - or,
worse, silently skipping the rest of the run.

Each entry script therefore hands over to a private, immediately unlinked copy of
itself (`inst/utils/chunk_pin_self.sh`) before doing any work. A job that has already
started is unaffected by later edits, reinstalls or `git pull`s of that file. Set
`PANNAGRAM_NO_PIN=T` to disable the handover (tracing, profiling).

Two things the pin does **not** cover:

* only the entry script is pinned. The R/Python workers and the `inst/utils/*.sh`
  chunks are read from the installed package at the moment each step starts, so
  reinstalling mid-run still changes the behaviour of the steps that have not
  started yet.
* the pipeline's own state. Checkpoints are keyed on the step **id**, not the step
  number, so inserting or removing a stage no longer invalidates the markers of an
  in-flight project; projects made by older versions are migrated on the next run
  (`migrate_legacy_step_markers` in `inst/utils/utils_bash.sh`).
