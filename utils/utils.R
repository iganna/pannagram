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
splitSeq <- function(sequence, n = 5000, step = 0) {
  
  # Split the sequence into individual characters
  sst <- seq2nt(sequence)
  if(step != 0){
    sst = sst[-(1:step)]
  }
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
#' Convert Aligned Sequences to Nucleotide Matrix
#'
#' This function takes a vector of aligned sequences and converts them into a matrix,
#' where each row represents one sequence converted into nucleotides.
#'
#' @param s.aln A list of aligned sequences.
#' @return A matrix where each row represents a sequence from `s.aln` converted into nucleotides.
#' Each sequence is represented as a row.
#' @examples
#' 
#' s.aln <- c("ACGTACTTACGGACGT", "ACGCACGTACGTACAT", "ACGTGCGTACGTTCGT")
#' seqs.mx <- aln2mx(s.aln)
#' msaplot(seqs.mx)
#' 
aln2mx <- function(s.aln) {
  aln.len = unique(nchar(s.aln))
  if(length(aln.len) != 1) stop('Aligned sequences have different lengths')
  if(is.null(names(s.aln))){
    seq.names = paste('s', 1:length(s.aln), sep = '')
  } else {
    seq.names = names(s.aln)
  }
  seqs.mx = matrix('-', length(s.aln), aln.len, dimnames = list(seq.names, NULL))
  for(i.s in 1:length(s.aln)){
    seqs.mx[i.s,] = seq2nt(s.aln[i.s])
  }
  return(seqs.mx)
}

#' ----------------------------------------------------------------------
#' Translate Nucleotide Sequence to Amino Acid Sequence
#'
#' This function translates a given nucleotide sequence into its corresponding
#' amino acid sequence based on the standard genetic code. Non-standard codons are
#' represented by a '-'.
#'
#' @param seq A character string representing the nucleotide sequence to be translated.
#' @return A character vector representing the translated amino acid sequence.
#' @examples
#' translateSeq('atgctctgccagtgccacggcggaagcgacaaagccBBB')
#' 
translateSeq <- function(seq) {
  # Define the genetic code as a named vector
  genetic.code <- c(
    "TTT" = "F", "TTC" = "F", "TTA" = "L", "TTG" = "L",
    "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
    "ATT" = "I", "ATC" = "I", "ATA" = "I", "ATG" = "M",
    "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
    "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
    "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
    "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
    "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
    "TAT" = "Y", "TAC" = "Y", "TAA" = "*", "TAG" = "*",
    "CAT" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
    "AAT" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
    "GAT" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
    "TGT" = "C", "TGC" = "C", "TGA" = "*", "TGG" = "W",
    "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
    "AGT" = "S", "AGC" = "S", "AGA" = "R", "AGG" = "R",
    "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"
  )
  
  # Convert sequence to uppercase to ensure matching
  seq <- toupper(seq)
  seq.len = nchar(seq) - 2
  
  # Extract codons from the sequence
  codons <- sapply(seq(1, seq.len, by = 3), function(i) {
    substr(seq, i, i+2)
  })
  
  # Translate codons to amino acids
  aa.seq = genetic.code[codons]
  
  # Replace NA values with '-', indicating unknown codons
  aa.seq[is.na(aa.seq)] = '-'
  # aa.seq = aa.seq[!is.na(aa.seq)]
  # if(length(aa.seq) == 0) return(NULL)
  
  # Remove names for a clean output
  names(aa.seq) = NULL
  
  return(aa.seq)
}

#' ----------------------------------------------------------------------
#' Identify and Extract Open Reading Frames (ORFs) from a Nucleotide Sequence
#'
#' This function identifies open reading frames (ORFs) within a nucleotide sequence by translating it
#' to an amino acid sequence and finding sequences between stop codons that meet a minimum length criterion.
#' The sequence name, if available, is preserved and appended with ORF identification and positional information.
#'
#' @param seq A character string or named character string representing the nucleotide sequence.
#' @param orf.min.len An integer value specifying the minimum length for an ORF to be considered valid. Default is 25.
#'
#' @return A list with two components:
#' \itemize{
#'   \item{\code{pos}:}{A matrix with two columns indicating the start and end positions of each ORF within the amino acid sequence.}
#'   \item{\code{orf}:}{A named character vector where each element is an ORF sequence, and names provide the original sequence name (if available), ORF, start and end positions, and ORF length.}
#' }
#'
#' @examples
#' seq <- "ATGGCCATGGCCCCCGCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTAG"
#' result <- seq2orf(seq)
#' print(result$pos)
#' print(result$orf)
#'
#' @export
#' @importFrom stringr str_detect
seq2orf <- function(seq, orf.min.len = 25){
  # Handle sequence naming
  if(is.null(names(seq))){
    seq.name = 'seq' 
  } else {
    seq.name = names(seq)[1]
  }
  
  # Translate the nucleotide sequence to amino acids
  aa.seq = translateSeq(seq)
  if(is.null(aa.seq)){
    return(list(pos = NULL, orf = NULL))  
  }
  
  seq.len = length(aa.seq)
  
  # Identify stop codon positions and calculate ORF lengths
  idx.stop = c(0, which(aa.seq == '*'), seq.len + 1)
  block.len = diff(idx.stop) - 1
  
  # Filter ORFs based on minimum length
  idx.remain = which(block.len >= orf.min.len)
  if(length(idx.remain) == 0){
    return(list(pos = NULL, orf = NULL))  
  }
  
  pos.orf = cbind(idx.stop[idx.remain]+1, idx.stop[idx.remain+1]-1)
  
  # Convert ORF positions from amino acids to nucleotides
  pos.orf.nt = data.frame(cbind((pos.orf[,1]-1)*3+1, pos.orf[,2]*3))
  colnames(pos.orf.nt) = c('beg', 'end')
  
  # Extract ORF sequences
  aa.orf = c()
  idx.remove = c()
  for(i in 1:nrow(pos.orf)){
    pos1 = pos.orf[i, 1]
    pos2 = pos.orf[i, 2]
    seq.aa = aa.seq[pos1:pos2]
    seq.aa = seq.aa[seq.aa != '-']
    if(length(seq.aa) < orf.min.len) {
      seq.tmp = ''
      idx.remove = c(idx.remove, i)
    } else {
      seq.tmp = paste0(seq.aa, collapse = '')
    }
    aa.orf[i] = seq.tmp
  }
  
  if(length(idx.remove) > 0){
    aa.orf = aa.orf[-idx.remove]
    pos.orf = pos.orf[-idx.remove,,drop=F]
    pos.orf.nt = pos.orf.nt[-idx.remove,,drop=F]
  }
  
  if(length(aa.orf) == 0){
    return(list(pos = NULL, orf = NULL))  
  }
  
  # Return ORF positions and sequences
  return(list(pos = pos.orf.nt, orf = aa.orf))  
}


#' ----------------------------------------------------------------------
#' Identify Open Reading Frames (ORFs) in Both DNA Strands
#'
#' This function searches for open reading frames (ORFs) in both the forward and reverse complement
#' of the given DNA sequence. It detects ORFs in all three possible reading frames for each strand,
#' calculates their positions, lengths, and on which strand they are located.
#'
#' @param seq.init A character string representing the DNA sequence.
#' @param seq A character string or named character string representing the nucleotide sequence.
#' @return A list with two elements:
#'   \item{pos}{A data frame of ORF positions including columns for start (`beg`), end (`end`), length in nucleotides (`len`), 
#'   length in amino acids (`aalen`), and strand (`strand`).}
#'   \item{orf}{A character vector of the amino acid sequences for each ORF, named with their positions, strand, and length.}
#'
#' @examples
#' seq <- "ATGGCCATGGCCCCCGCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTAG"
#' orfFinderResult <- orfFinder(seq)
#' print(orfFinderResult$pos)
#' print(orfFinderResult$orf)
#'
#' @export
orfFinder <- function(seq.init, orf.min.len = 25){
  seq.len = nchar(seq.init) # Calculate the length of the initial sequence
  seq = seq.init
  
  orfs = c() # Initialize vector to store ORFs
  pos = c() # Initialize vector to store positions
  
  # Iterate over both strands
  for(i.strand in 1:2){
    # Iterate over the three possible reading frames
    for(i.orf in 0:2){
      # Adjust sequence for reading frames 2 and 3
      if(i.orf != 0){
        seq <- substr(seq, 2, nchar(seq))  
      }
      
      # Find ORFs in the current frame and strand
      orf.res = seq2orf(seq, orf.min.len = orf.min.len)
      if(is.null(orf.res$pos)) next
      
      # Adjust positions based on the reading frame and strand
      pos.tmp = orf.res$pos + i.orf
      if(i.strand == 2){
        pos.tmp = seq.len - pos.tmp + 1
      }
      pos.tmp$shift = i.orf
      
      # Store ORFs and their positions
      orfs = c(orfs, orf.res$orf)
      pos = rbind(pos, pos.tmp)
      
      if(nrow(pos) != length(orfs)) stop('Something went wrong')
      
    }
    # Reverse and complement the sequence for the second strand
    seq = revComplSeq(seq.init)
  }
  
  if(is.null(pos)){
    return(list(pos = NULL, orf = NULL))
  }
  
  # Calculate additional position info and order by ORF length
  pos$len = abs(pos[,2] - pos[,1]) + 1
  pos$aalen = nchar(orfs)
  pos$strand = c('+', '-')[(pos[,1] > pos[,2]) + 1]
  idx.order = order(-pos$len)
  pos = pos[idx.order,]
  orfs = orfs[idx.order]
  
  # Name ORFs with position, strand, and length info
  names(orfs) = paste('ORF', pos[,1], pos[,2], pos$strand, 'aaLEN', pos$aalen, sep = '|')
  
  # Return positions and ORFs
  return(list(pos = pos, orf = orfs))
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
#' Convert a Nucleotide Sequence to its Reverse Complement
#'
#' This function takes a nucleotide sequence, converts it into its reverse complement,
#' and returns the result as a single string.
#'
#' @param seq A character string representing the nucleotide sequence.
#' @return A character string representing the reverse complement of the input sequence.
#'
#' @examples
#' revComplSeq("ATGC")
#'
#' @author Anna A. Igolkina
revComplSeq <- function(seq){
  
  seq = seq2nt(seq)
  seq.rc = revCompl(seq)
  seq.rc = nt2seq(seq.rc)
  
  return(seq.rc)
}


#' ----------------------------------------------------------------------
#' Generate the complement of a vectorized sequence
#'
#' @param s A character vector where each element is a single nucleotide.
#' @return A character vector representing the reverse complement of the input.
#' 
#' @examples
#' justCompl(c("A", "T", "G", "C"))
#' 
#' @author Anna A. Igolkina 
#' 
justCompl <- function(s){
  
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
  
  seqs.c = complementary_nts[s]
  if(sum(is.na(seqs.c)) != 0) stop('Wrong nucleotides are provided')
  return(seqs.c)
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
  
  arguments <- paste('* ', paste(arguments_list, collapse = " "), sep = '')
  
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
  
  arguments <- paste('  Attention: ', paste(arguments_list, collapse = " "), sep = '')
  
  text.color <- make_style("#FC345C")
  # bg <- make_style("grey5", bg = TRUE)
  fancy <- combine_styles(text.color)
  # message(arguments)
  cat(fancy(arguments))
  cat('\n')
}


#' ----------------------------------------------------------------------
#' Display general stylized messages
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


#' ----------------------------------------------------------------------
#' Safe Removal of Variable from Global Environment
#'
#' This function safely removes a specified variable from the global environment
#' if it exists. It avoids errors associated with trying to remove non-existent variables.
#' The function uses the variable itself as an argument and internally determines its name.
#'
#' @param var Variable to be removed (provided as the variable itself, not as a string).
#'
#' @examples
#' x <- 42
#' rmSafe(x)
#' rmSafe(xxx)  # no error
#' # x is now removed from the global environment
#'
rmSafe <- function(var) {
  var.name = deparse(substitute(var))
  if (exists(var.name, envir = globalenv())) {
    rm(list = var.name, envir = globalenv())
  }
}

#' ----------------------------------------------------------------------
#' Convert BLAST results to GFF format
#'
#' @description
#' `blastres2gff` converts a data frame containing BLAST results into a GFF formatted file.
#'
#' @param v.blast A data frame containing the BLAST results. Expected columns are V1, V4, V5, V7, V8, len1, and strand.
#' @param f.gff The file path where the GFF output will be saved.
#' @param to.sort Boolean value indicating whether the output should be sorted. If `TRUE` (default), 
#' the output is sorted by column 4 and then by column 1. If `FALSE`, the output is not sorted.
#'
#' @return This function does not return a value. It writes the GFF formatted data to a file specified by `f.gff`.
#'
#' @examples
#' # Example usage (assuming `blast_results` is your data frame with BLAST results):
#' blastres2gff(blast_results, "output.gff")
#'
blastres2gff <- function(v.blast, f.gff, to.sort = T){
  v.gff = data.frame(col1 = v.blast$V8,
                     col2 = 'blast2gff',
                     col3 = 'query',
                     col4 = v.blast$V4,
                     col5 = v.blast$V5,
                     col6 = '.',
                     col7 = v.blast$strand,
                     col8 = '.',
                     col9 = paste('ID=Q', 1:nrow(v.blast),
                                  ';query=',v.blast$V1,
                                  ';len=', v.blast$len1,
                                  ';coverage=', v.blast$V7, sep = ''))
  
  # Sorting
  if(to.sort){
    v.gff = v.gff[order(v.gff$col4),]
    v.gff = v.gff[order(v.gff$col1),]
  }
  
  write.table(v.gff, f.gff, sep = '\t', quote = F, row.names=F, col.names = F)
}


#' ----------------------------------------------------------------------
#' Calculate the Repeat Score of a String
#'
#' @description
#' `repeatScore` computes the repeat score for a given string based on the frequency 
#' of repeating substrings within the string. The repeat score is the ratio of the 
#' total count of repeating substrings (exceeding a specified cutoff) to the total 
#' number of substrings.
#'
#' @param s The input string for which the repeat score is to be calculated.
#' @param wsize Window size for substring generation (default is 11). This determines 
#' the length of each substring extracted from the input string for analysis.
#' @param dup.cutoff Cutoff for considering a substring as repeating (default is 2). 
#' Substrings that appear more times than this cutoff contribute to the repeat score.
#'
#' @return A single numeric value representing the repeat score of the input string.
#' This score is the ratio of the count of substrings appearing more than `dup.cutoff` 
#' times to the total number of substrings generated.
#'
#' @examples
#' library(stringi)  # Ensure the stringi package is loaded
#' s = 'TAAACCCTAAACCCTAAACCCTAAACCCTAACCCTAAACCCTAAACCCTAAACCCTAAACCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCC'
#' repeatScore(s)
#'
#' @export
#'
repeatScore <- function(s, wsize = 11, dup.cutoff = 2){
  
  substrings <- unlist(stringi::stri_sub_all(s, 1:(nchar(s)-wsize+1), wsize -1 + 1:(nchar(s)-wsize+1)))
  substrings = sort(substrings)
  cnt <- rle(substrings)
  cnt = cnt$lengths
  
  return(sum(cnt[cnt > dup.cutoff]) / length(substrings))
}
# SLOW VERSION:
# old_repeatScore <- function(s, wsize = 11, dup.cutoff = 2){
#   substrings <- sapply(1:(nchar(s) - wsize + 1), function(i) {
#     substring(s, i, i + wsize - 1)
#   })
#   cnt = table(substrings)
#   return(sum(cnt[cnt > dup.cutoff]) / length(substrings))
# }


#' ----------------------------------------------------------------------
#' Convert Combination String to Reference Chromosome Number
#'
#' This function parses a combination string of the format 'X_Y' where X and Y are numbers.
#' It extracts the second number (Y) from the combination and returns it as a numeric value.
#' The function is designed to work with chromosome combinations.
#'
#' @param s.comb A string representing a chromosome combination in the format 'X_Y'.
#'               Both X and Y must be numeric values.
#'
#' @return An integer representing the second part (Y) of the combination.
#'
#' @examples
#' comb2ref("3_5") # returns 5
#' 
comb2ref <- function(s.comb){
  pattern <- "^\\d+_\\d+$"
  if(!grepl(pattern, s.comb)) stop('Something is wrong with the bombinations of chromosomes')
  ref.chr = as.numeric(strsplit(s.comb, '_')[[1]][2])
  return(ref.chr)
}


#' ----------------------------------------------------------------------
#' Load or Install and Load an R Package
#'
#' This function attempts to load an R package. If the package is not already installed,
#' it tries to install the package from CRAN and then loads it. This is useful for ensuring
#' that scripts or analyses have access to required packages with minimal manual intervention.
#'
#' @param package The name of the package to load. This should be a character string.
#'
#' @return Invisible NULL. This function is called for its side effect of loading a package
#' or displaying a message if the package cannot be installed.
#'
#' @examples
#' load.library("ggplot2") # Attempts to load ggplot2, installing it first if necessary
#' 
#' @author Anna Igolkina
load.library <- function(package) {
  # Attempt to load the package
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    # If the package is not installed, attempt to install it
    install_result <- tryCatch({
      install.packages(package, quiet = TRUE)
      TRUE  # Return TRUE if the installation is successful
    }, error = function(e) {
      FALSE  # Return FALSE if an error occurs during installation
    })
    
    # Check if the installation was successful
    if (!install_result) {
      message(paste("Failed to install package '", package, 
                    "'. Check if it is available for your version of R and if the package name is correct.", sep=""))
    } else {
      # Attempt to load the package after installation
      library(package, character.only = TRUE)
    }
  }
}


#' ----------------------------------------------------------------------
#' Find Sequences of Ones in Binary Vector
#'
#' Identifies and returns the start and end indices of consecutive sequences of ones in a binary vector.
#'
#' @param g.bin A vector containing zeros and ones.
#' 
#' @return A data frame with columns `beg` and `end`.
#' 
#' @examples
#' g.bin <- c(0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0)
#' findOnes(g.bin)
#' 
#'     beg end
#'   1   2   3
#'   2   5   5
#'   3   8  10
#'   
#' @export
findOnes <- function(g.bin) {
  changes <- diff(c(0, g.bin, 0))
  beg <- which(changes == 1)
  end <- which(changes == -1) - 1
  return(data.frame(beg = beg, end = end))
}

#' ----------------------------------------------------------------------
#' Finds Duplicates and return unique values
#'
#' @param x Vector to search for duplicates.
#' @return Vector of unique duplicated values from `x`.
#' @examples
#' uniqueDuplicates(c(1, 2, 2, 3, 4, 4, 5))
#' @export
uniqueDuplicates <- function(x){
  return(unique(x[duplicated(x)]))
}


#' ----------------------------------------------------------------------
#' Finds indexes of all non-unique values
#'
#' @param x Vector to search for duplicates.
#' @return Indexes of duplicated values from `x`.
#' @export
idxDuplicates <- function(x){
  y = unique(x[duplicated(x)])
  return(which(x %in% y))
}

#' ----------------------------------------------------------------------
#' Write GFF Data to File
#'
#' Writes first 9 columns of the table to the gff file.
#'
#' @param gff A data frame representing GFF data.
#' @param file The path to the file where the GFF data should be written.
#' @export
#' @examples
writeGFF <- function(gff, file){
  write.table(gff[,1:9], file, quote = F, row.names = F, col.names = F, sep = '\t')
}


#' ----------------------------------------------------------------------
#' Calculate Moving Window Mean
#'
#' This function calculates the moving window mean for a given numeric vector.
#' The window size is adjustable, allowing for flexibility in how the mean is calculated over the vector.
#'
#' @param v A numeric vector
#' @param wnd.size The size of the window over which the mean is calculated. Default is 10000.
#'
#' @return A numeric vector containing the moving window mean values.
#' @examples
#' v <- 1:100000
#' wndMean(v, wnd.size = 10000)
#' @export
wndMean <- function(v, wnd.size = 10000){
  
  # Calculate the length of the input vector
  len.v <- length(v)
  # Determine the number of full parts within the window size
  num.parts <- len.v %/% wnd.size
  # Initialize a vector for storing the mean values
  wnd.mean <- numeric(num.parts + ifelse(len.v %% wnd.size > 0, 1, 0))
  
  # Calculate the mean for each full part
  for(i in 1:num.parts) {
    wnd.mean[i] <- mean(v[((i-1) * wnd.size + 1) : (i * wnd.size)])
  }
  
  # Calculate the mean for the last part if it is smaller than the window size
  if(len.v %% wnd.size > 0) {
    wnd.mean[num.parts + 1] <- mean(v[(num.parts*wnd.size + 1):len.v])
  }
  
  return(wnd.mean)
}


#' Read BLAST Output File
#'
#' Reads a BLAST output file and returns its contents as a data frame. If the file contains
#' no data (only comments or empty), `NULL` is returned. It assumes that the file does not
#' have a header and all strings should not be converted to factors.
#'
#' @param file The BLAST output file to be read.
#' @return A data frame containing the contents of the BLAST file if it contains data; otherwise, `NULL`.
readBlast <- function(file) {
  if (any(grepl("^[^#]", readLines(file)))) {
    return(read.table(file, stringsAsFactors = F,  header = F))
  } else {
    return(NULL)
  }
}

readTableMy <- function(file){
  pokaz('Replace function readTableMy to readBlast')
  return(readBlast(file))
}


showt <- function(t, irow=NULL, chr=F, add = c()){
  idx = c(1:5,7, add)
  idx = intersect(idx, 1:ncol(t))
  irow = irow[irow <= nrow(t)]
  
  if(chr) idx = c(idx, 10)
  if(is.null(irow)){
    print(t[1:min(100, nrow(t)),idx])
  } else {
    print(t[irow,idx])
  }
}






