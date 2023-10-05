library(crayon)


#' Display styled text in console.
#'
#' This function takes one or more arguments, combines them into a single string,
#' and then displays the string with a custom style in the console.
#' 
#' @param ... One or more arguments to combine into a single string.
#' 
#' @return Styled string printed to the console.
#' 
#' @export
pokazStage <- function(...) {
  arguments <- paste('*', paste(..., sep = " "), sep = ' ')
  
  text.color <- make_style("#34FCFC")
  bg <- make_style("grey5", bg = TRUE)
  fancy <- combine_styles(text.color)
  # message(arguments)
  cat(fancy(arguments))
  cat('\n')
}

pokazAttention <- function(...) {
  arguments <- paste('  Attention:', paste(..., sep = " "), sep = ' ')
  
  text.color <- make_style("#FC345C")
  bg <- make_style("grey5", bg = TRUE)
  fancy <- combine_styles(text.color)
  # message(arguments)
  cat(fancy(arguments))
  cat('\n')
}


pokaz <- function(...) {
  arguments_list <- list(...)
  # Check if any arguments are vectors
  for (i in seq_along(arguments_list)) {
    if (is.character(arguments_list[[i]]) && length(arguments_list[[i]]) > 1) {
      arguments_list[[i]] <- paste(arguments_list[[i]], collapse = " ")
    }
  }
  
  arguments <- paste('  ', paste(arguments_list, collapse = " "), sep = '')
  
  text.color <- make_style("#FDFDFD")
  bg <- make_style("grey5", bg = TRUE)
  fancy <- combine_styles(text.color)
  # message(arguments)
  cat(fancy(arguments))
  cat('\n')
}