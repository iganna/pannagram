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

