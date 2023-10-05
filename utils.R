library(crayon)

#' ----------------------------------------------------------------------
#' Display stylized stage messages
#'
#' This function displays a stylized message indicating the stage or step of a process.
#'
#' @param ... Arguments to be concatenated into a message string.
#'
#' @return No return value, called for side effects.
#' 
pokazStage <- function(...) {
  
  arguments_list <- list(...)
  # Check if any arguments are vectors
  for (i in seq_along(arguments_list)) {
    if (is.character(arguments_list[[i]]) && length(arguments_list[[i]]) > 1) {
      arguments_list[[i]] <- paste(arguments_list[[i]], collapse = " ")
    }
  }
  arguments <- paste('*', paste(..., sep = " "), sep = ' ')
  
  text.color <- make_style("#34FCFC")
  # bg <- make_style("grey5", bg = TRUE)
  fancy <- combine_styles(text.color)
  # message(arguments)
  cat(fancy(arguments))
  cat('\n')
}

#' ----------------------------------------------------------------------
#' Display attention messages with stylization
#'
#' This function displays a stylized attention message to alert the user.
#'
#' @param ... Arguments to be concatenated into an attention message string.
#'
#' @return No return value, called for side effects.
#'
pokazAttention <- function(...) {
  
  arguments_list <- list(...)
  # Check if any arguments are vectors
  for (i in seq_along(arguments_list)) {
    if (is.character(arguments_list[[i]]) && length(arguments_list[[i]]) > 1) {
      arguments_list[[i]] <- paste(arguments_list[[i]], collapse = " ")
    }
  }
  
  arguments <- paste('  Attention:', paste(..., sep = " "), sep = ' ')
  
  text.color <- make_style("#FC345C")
  # bg <- make_style("grey5", bg = TRUE)
  fancy <- combine_styles(text.color)
  # message(arguments)
  cat(fancy(arguments))
  cat('\n')
}


#' ----------------------------------------------------------------------
#' #' Display general stylized messages
#'
#' This function displays a general stylized message to the user.
#'
#' @param ... Arguments to be concatenated into a message string.
#'
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
  # bg <- make_style("grey5", bg = TRUE)
  fancy <- combine_styles(text.color)
  # message(arguments)
  cat(fancy(arguments))
  cat('\n')
}