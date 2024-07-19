#' Save a ggplot2 plot to a PDF file
#'
#' @param geom A ggplot2 plot object to be saved
#' @param path The directory where the PDF file will be saved (default is current working directory)
#' @param name The base name of the PDF file (default is "some")
#' @param width The width of the PDF file in inches (default is 10)
#' @param height The height of the PDF file in inches (default is 10)
#' @export
savePDF <- function(geom, path = '.', name='some', width=10, height=10){
  pdf(file.path(path, paste0(name, ".pdf")), width=width, height=width)
  print(geom)
  dev.off()
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
#'   get_prefixes(".")
#' }
#' @export
get_prefixes <- function(path) {
  files <- list.files(path)
  files_with_underscore <- files[grepl("_", files)]
  split_files <- strsplit(files_with_underscore, "_")
  prefixes <- sapply(split_files, function(x) x[1])
  return(unique(prefixes))
}


#' Find Genome Files in Folder
#'
#' This function searches for files in the specified directory that start with a given prefix 
#' and end with one of the specified extensions.
#'
#' @param genome.pref A character string representing the prefix of the file names.
#' @param path.ref A character string representing the directory where the files are searched.
#' @param ext A character vector of file extensions to search for. Default is c('fasta', 'fna', 'fa', 'fas').
#' @return A character vector with the names of the matching files. If no files are found, an error is raised.
#' @examples
#' ref <- "ref"
#' path.ref <- "."
#' findGenomeFile(ref, path.ref)
#' @export
findGenomeFile <- function(genome.pref, path.genome, ext = c('fasta', 'fna', 'fa', 'fas')) {
  # Create the pattern
  ext.pattern <- paste(ext, collapse = "|")
  pattern <- paste0("^", genome.pref, ".*\\.(", ext.pattern, ")$")
  
  # Search for files in the specified directory
  ref.files <- list.files(path = path.genome, pattern = pattern, full.names = TRUE)
  
  # Check if any files were found and return the result or raise an error
  if (length(ref.files) == 0) {
    return(NULL)
  } else {
    return(ref.files)
  }
}
