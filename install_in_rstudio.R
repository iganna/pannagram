#!/usr/bin/env Rscript

if (!dir.exists("R")) {
  message("Creating R/ directory...")
  dir.create("R")
}

r_files <- list.files(
  path       = "inst",
  pattern    = "\\.R$",
  recursive  = TRUE,
  full.names = TRUE
)

has_export <- vapply(r_files, function(path) {
  lines <- readLines(path, warn = FALSE)
  any(grepl("@export", lines, fixed = TRUE))
}, logical(1))
exported_files <- r_files[has_export]

for (src in exported_files) {
  dest <- file.path("R", basename(src))
  file.symlink(normalizePath(src), dest)
}
message("Linked ", length(exported_files), " files into R/")

message("Documenting package...")
devtools::document()

message("Installing package...")
devtools::install()

message("All done!")
