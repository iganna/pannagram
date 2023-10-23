suppressMessages({
if (!require(crayon)) {
  install.packages("crayon")
}

library(crayon)
  
})

#' ----------------------------------------------------------------------
#' Read FASTA File
#'
#' This function reads a FASTA file and returns a named vector of sequences.
#' 
#' 
#' This function 3 times faster than the standard reading:
#' start_time <- Sys.time()
#' q.fasta = read.fasta(file.fasta)
#' end_time <- Sys.time()
#' elapsed_time <- end_time - start_time
#' print(elapsed_time)
#' 
#' @param file.fasta A string, the path to the FASTA file to be read.
#' 
#' @return A named vector where the names are the sequence headers 
#' 
#' @author Anna A. Igolkina 
#' 
readFastaMy <- function(file.fasta){
  
  start_time <- Sys.time()
  file.content <- readLines(file.fasta)
  header.idx <- c(which(substr(file.content, 1, 1) == ">"), length(file.content) + 1)
  n.seq = length(header.idx) - 1
 
  sequences <- c()
  for(i.seq in 1:n.seq){
    sequences = c(sequences, 
                  paste0(file.content[(header.idx[i.seq]+1):(header.idx[i.seq+1]-1)], collapse = ''))
  }
  seq.names = file.content[header.idx[1:n.seq]]
  seq.names = substr(seq.names, 2, nchar(seq.names))
  names(sequences) = seq.names
  
  rm(file.content)
  
  return(sequences)
}

#' ----------------------------------------------------------------------
#' Write Sequences to a FASTA File
#'
#' This function writes a named vector of sequences to a FASTA file. Each sequence
#' in the vector will be written to the file with its corresponding name as a header.
#'
#' @param sequences A vector of sequences to be written to the FASTA file.
#' #' @param file A string specifying the path to the output FASTA file.
#' @param seq.pref A string used as a prefix for sequence names when names are not provided
#' in the sequences vector. Default is "seq". Names will be constructed as seq.pref_1, seq.pref_2, etc.
#' @param append Logical. If TRUE, sequences will be appended to the file (if it exists). 
#' If FALSE, the function will overwrite the file with the new sequences. Default is FALSE.
#'
#' @return Invisible NULL.
#' 
#' @author Anna A. Igolkina
#'  
writeFastaMy <- function(sequences, file, seq.names = NULL, pref = NULL, append = FALSE){
  
  if(is.null(seq.names)){
    if(!is.null(names(sequences))){
      seq.names = names(sequences)
    } else if(!is.null(pref)){
      seq.names = paste(pref, 1:length(sequences), sep = '')
    } else {
      pref = 'seq_'
      seq.names = paste(pref, 1:length(sequences), sep = '')
    }
  }
  
  
  # Open the file connection for writing or appending
  mode <- ifelse(append, "a", "w")
  con <- file(file, mode)
  
  # Write each sequence to the file in FASTA format
  for (i in 1:length(sequences)) {
    cat(sprintf(">%s\n", seq.names[i]), file = con)
    cat(sprintf("%s\n", sequences[i]), file = con)
  }
  
  # Close the file connection
  close(con)

}

#' ----------------------------------------------------------------------
#' Split a Sequence into Fixed-Length Chunks
#'
#' This function splits a given sequence into chunks of a specified length. If the sequence's length 
#' is not a multiple of the specified chunk length, the last chunk will be padded with empty characters 
#' to match the specified length.
#'
#' @param sequence A character string representing the sequence to be split.
#' @param n Numeric. The desired length for each chunk. Default is 5000.
#'
#' @return A vector of character strings, where each string represents a chunk of the sequence.
#'
#' @examples
#' \dontrun{
#' seq <- "ACGTAnna0123456789?"
#' splitSeq(seq, n=4)
#'}
#'
#' @author Anna A. Igolkina 
#' 
splitSeq <- function(sequence, n = 5000) {
  
  # Split the sequence into individual characters
  sst <- strsplit(sequence, '')[[1]]
  n.add <- ceiling(length(sst) / n) * n - length(sst)
  sst <- c(sst, rep('', n.add))
  
  # Convert the sequence to a matrix with 'n' rows
  # This facilitates chunking the sequence into pieces of length 'n'
  m <- matrix(sst, nrow = n)
  
  # Convert each column of the matrix back to a string
  # Each column represents a chunk of the sequence
  s.chunks <- apply(m, 2, paste0, collapse = '')
  
  return(s.chunks)
}

#' ----------------------------------------------------------------------
#' Convert a string to a vector of nucleotides
#'
#' @param s A single string representing a sequence.
#' 
#' @return A character vector where each element is a single nucleotide.
#' 
#' @author Anna A. Igolkina 
#' 
seq2nt <- function(s){
  if(length(s) != 1) stop('There should be only one sequence')
  return(strsplit(s, '')[[1]])
}


#' ----------------------------------------------------------------------
#' Convert a vector of nucleotides to a string
#'
#' @param s A character vector where each element is a single nucleotide.
#' @return A single string representing a sequence.
#' 
#' @author Anna A. Igolkina 
#' 
nt2seq <- function(s){
  return(paste0(s, collapse = ''))
}


#' ----------------------------------------------------------------------
#' Generate the reverse complement of a vectorized sequence
#'
#' @param s A character vector where each element is a single nucleotide.
#' @return A character vector representing the reverse complement of the input.
#' 
#' @examples
#' revCompl(c("A", "T", "G", "C"))
#' 
#' @author Anna A. Igolkina 
#' 
revCompl <- function(s){
  
  if(nchar(s[1]) != 1) stop('Sequence should be vectorised by nucleotides')
  complementary_nts <- c(
    A='T', T='A', C='G', G='C', 
    R='Y', Y='R', S='S', W='W', 
    K='M', M='K', B='V', D='H', 
    H='D', V='B', N='N', 
    a='T', t='A', c='G', g='C', 
    r='Y', y='R', s='S', w='W', 
    k='M', m='K', b='V', d='H', 
    h='D', v='B', n='N'
  )
  
  seqs.rc = rev(complementary_nts[s])
  if(sum(is.na(seqs.rc)) != 0) stop('Wrong nucleotides are provided')
  return(seqs.rc)
}


#' ----------------------------------------------------------------------
#' Display stylized stage messages
#'
#' This function displays a stylized message indicating the stage or step of a process.
#'
#' @param ... Arguments to be concatenated into a message string.
#'
#' @return No return value, called for side effects.
#' 
#' @author Anna A. Igolkina 
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
#' @author Anna A. Igolkina 
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
#' @author Anna A. Igolkina 
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