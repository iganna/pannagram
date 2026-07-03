# ---- Logging ----

log.level <- opt$log.level

# Check if the log level is provided
if (!is.null(log.level)) {
  log.level <- suppressWarnings(as.numeric(log.level))
  if (is.na(log.level)) {  # If conversion fails, default to 0
    log.level <- 0
  }
} else {
  log.level <- 0  # If no log level is provided, default to 0
}

# Define logging levels for the main code and the code in the loop (parallel processes)
ll.main <- 2
ll.loop <- 3 

# Determine if the main code should be echoed
echo.main <- log.level >= ll.main

# Determine if the code in the loop should be echoed
echo.loop <- log.level >= ll.loop

path.log <- opt$path.log
if (!is.null(path.log) & !is.null(log.level)) {
  if (!dir.exists(path.log)) dir.create(path.log)
  
  # Remove old log files
  log_files <- list.files(path = path.log, pattern = "\\.log$", full.names = TRUE)
  # file.remove(log_files)
  
  file.log.main <- paste0(path.log, 'script.log')
  invisible(file.create(file.log.main))
} else {
  path.log <- NULL
  file.log.main <- NULL
}

# ---- Loop logging as checkpoints (minimal number of files) ----
#
# Idea: decouple "who writes" from "what is done".
#   * one log file per PARALLEL WORKER  -> the number of files is bounded by the
#     number of cores, not by the number of loop items;
#   * a "Done: <id>" marker per item    -> per-item checkpoint granularity is kept
#     inside those few files.
#
# On resume the set of completed items is rebuilt by scanning ALL worker logs, so a
# rerun with a DIFFERENT number of cores still sees every previously written marker.

# Prefix used for checkpoint markers inside worker logs
loop.done.prefix <- 'Done:'

# Return the per-worker log file, creating it if needed.
# Uses a stable worker id (.worker.id, assigned to each cluster worker below) so the
# file names are reused across runs; falls back to the PID for the sequential path.
initLoopLog <- function(path.log, worker.id = NULL) {
  if (is.null(path.log)) return(NULL)
  if (is.null(worker.id)) {
    worker.id <- if (exists('.worker.id', envir = .GlobalEnv)) {
      get('.worker.id', envir = .GlobalEnv)
    } else {
      Sys.getpid()
    }
  }
  f <- paste0(path.log, 'core_', worker.id, '.log')
  if (!file.exists(f)) invisible(file.create(f))
  return(f)
}

# Build the set of already-completed item ids by scanning every core_*.log.
# Independent of the number of cores (globs all worker logs).
getDoneSet <- function(path.log, prefix = loop.done.prefix) {
  if (is.null(path.log)) return(character(0))
  files <- list.files(path = path.log, pattern = '^core_.*\\.log$', full.names = TRUE)
  if (length(files) == 0) return(character(0))
  pat <- paste0('^\\s*', prefix, '\\s*')
  done <- unlist(lapply(files, function(f) {
    lines <- readLines(f, warn = FALSE)
    lines <- grep(pat, lines, value = TRUE)
    trimws(sub(pat, '', lines))
  }), use.names = FALSE)
  unique(done)
}

# Write the checkpoint marker that flags an item as fully processed.
markDone <- function(id, file, echo = FALSE) {
  pokaz(paste0(loop.done.prefix, ' ', id), file = file, echo = echo)
}

