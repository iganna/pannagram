#' Save a ggplot2 plot to a PDF file
#'
#' @param geom A ggplot2 plot object to be saved
#' @param path The directory where the PDF file will be saved (default is current working directory)
#' @param name The base name of the PDF file (default is "some")
#' @param width The width of the PDF file in inches (default is 10)
#' @param height The height of the PDF file in inches (default is 10)
#' @export
savePDF <- function(geom, path = '.', name='some', width=10, height=10){
  invisible(suppressMessages({
    pdf(file.path(path, paste0(name, ".pdf")), width=width, height=height)
    print(geom)
    dev.off()
  }))
}

#' Save a ggplot2 plot to a PNG file
#'
#' @param geom A ggplot2 plot object to be saved
#' @param path The directory where the PNG file will be saved (default is current working directory)
#' @param name The base name of the PNG file (default is "some")
#' @param width The width of the PNG file in inches (default is 10)
#' @param height The height of the PNG file in inches (default is 10)
#' @param res The resolution of the PNG file in dots per inch (default is 300)
#' @export
savePNG <- function(geom, path = '.', name='some', width=10, height=10, res=300){
  invisible(suppressMessages({
    png(file.path(path, paste0(name, ".png")), width=width, height=height, units='in', res=res)
    print(geom)
    dev.off()
  }))
}

#' Get unique file prefixes in a given path
#'
#' This function gets all the files in the given path, splits the file names by the `_` character,
#' and returns the unique prefixes. Files without a `_` character are ignored.
#'
#' @param path Character string specifying the path to the directory.
#'
#' @return A character vector of unique file prefixes.
#'
#' @examples
#' \dontrun{
#'   # Get unique prefixes in the current directory
#'   getPrefixes(".")
#' }
#' @export
getPrefixes <- function(path) {
  files <- list.files(path)
  files_with_underscore <- files[grepl("_", files)]
  split_files <- strsplit(files_with_underscore, "_")
  prefixes <- sapply(split_files, function(x) x[1])
  return(unique(prefixes))
}

#' In given directory find FASTA file with chosen basename (no matter what extension)
#'
#' Returns absolute path to fasta file, given its basename is unique in the directory
#' Throws descriptive error otherwise. No support for compressed FASTAs currently
#'
#' @param directory Character string specifying the path to the directory.
#' @param basename Character string with no fasta suffix
#' @return A character string with absolute file path
#'
#' @export
findGenomeFile <- function(directory, basename) {
  extensions <- c('.fasta', '.fna', '.fa', '.fas')
  extensions <- c(extensions, toupper(extensions))
  
  compressed.extensions <- c('.gz', '.zip', '.bz2', '.xz')
  compressed.extensions <- c(compressed.extensions, toupper(compressed.extensions))
  
  found.files <- character(0)
  
  for (ext in extensions) {
    filename <- file.path(directory, paste0(basename, ext))
    if (file.exists(filename)) {
      found.files <- c(found.files, filename)
    }
  }

  if (length(found.files) > 1) {
    stop(paste0("For basename '", basename, "' multiple FASTA files were found in '", directory, "':\n", 
               paste0(found.files, collapse = "\n")))
  }
  
  if (length(found.files) == 1) {
    return(found.files[1])
  }
  
  for (ext in extensions) {
    for (comp.ext in compressed.extensions) {
      compressed.filename <- file.path(directory, paste0(basename, ext, comp.ext))
      if (file.exists(compressed.filename)) {
        stop(paste0("Compressed file found:\n", compressed.filename, 
                      "\nCompressed FASTA files are not currently supported."))
      }
    }
  }
  
  return(NULL)
}
