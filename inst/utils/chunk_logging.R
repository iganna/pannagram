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

