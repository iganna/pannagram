# Handy Functions

The Pannagram R library provides basic tools for working with biological sequences.  
It helps you analyze, transform, and visualize sequences in an easy and convenient way.  
It works with plain strings and character vectors, offering a lightweight alternative to specialized containers like Biostrings,   
which can be powerful but may feel overly complex for basic tasks.


With the successful installation of the Pannagram environment, the Pannagram R library becomes available.  
From within this environment, load the library in R/RStudio by running:

```r
library(pannagram)
```

Here we provide a quick start guide to the functions.
For detailed information about function parameters, run `<function_name> -h`.

## Read and Write FASTA

### `readFasta()`
Reading a FASTA file returns a named character vector:
```R
seqs <- readFasta("input_sequences.fasta")

# Example output:
# names(seqs)
# [1] "seq1" "seq2"
# seqs["seq1"]
# [1] "ATGCGTACGTA"
```

### `writeFasta()`
Writing works the same way: provide a named vector of sequences:
```
writeFasta(seqs, "output_sequences.fasta")
```


## Nucleotide Content

### `ntplot`() 

Plots nucleotide density across sliding windows:
```
ntplot(s)
```

<div style="width: 40%;">
<p align="left">
  <img src="images/ntplot.png" style="width:100%; object-fit:cover;"/>
</p>
</div>


You can also plot separate densities for each nucleotide:
```
ntplot(seq, nt.separate = T)
```

<div style="width: 40%;">
<p align="left">
  <img src="images/ntplot_separate.png" style="width:100%; object-fit:cover;"/>
</p>
</div>


### `gcContent()`

Compute GC-content for a vector of DNA sequences.

```R
# Example
seqs <- c(s1 = "ATGCGATCGATCG", s2 = "GCGCGCGCGCGCG", s3 = "ATATATATATAT")

gcContent(seqs)

# Output:
#        s1        s2        s3 
# 0.5384615 1.0000000 0.0000000
```

### `mx2cons()`
Generate a consensus sequence (or multiple consensus-like rows) from a character alignment matrix.
```R
# Example
mx
#     [,1] [,2] [,3] [,4] [,5]
# s1 "A"  "T"  "G"  "-"  "C"
# s2 "A"  "T"  "G"  "G"  "C"
# s3 "A"  "C"  "G"  "-"  "C"

cons <- mx2cons(mx)
cons
# Output:
# [1] "A" "T" "G" "G" "C"
```

### `mx2profile()`
Compute a nucleotide frequency profile for each alignment column.
```R
mx
#     [,1] [,2] [,3] [,4] [,5]
# s1 "A"  "T"  "G"  "-"  "C"
# s2 "A"  "T"  "G"  "G"  "C"
# s3 "A"  "C"  "G"  "-"  "C"

prof <- mx2profile(mx)
prof
# Output (no gaps counted):
#     [,1] [,2] [,3] [,4] [,5]
# A     3    0    0    0    0
# C     0    1    0    0    3
# G     0    0    3    1    0
# T     0    2    0    0    0
```

## Save Figures
All visualization functions in the Pannagram R package return a ggplot2 object,  
which can be exported using the following figure-saving helpers.

### `savePDF()`
```
library(ggplot2)

p <- ntplot(s)

savePDF(p, path = "folder_with_figures", name = "ntplot", width = 6, height = 4)
```

### `savePNG()`
```
library(ggplot2)

p <- dotself(s)

savePNG(p, path = "folder_with_figures", name = "dotplot", width = 6, height = 4)
```


## Sequence and Alignment Conversion Functions

### `revCompl()`
```
seq <- "ATGCRY"

rc <- revCompl(seq)
rc
# Output:
# [1] "RYGCAT"
```

### `seq2nt()`
Convert a string into a vector of individual characters (nucleotides).
```R
# Example
seq <- "ATGCGATC"

nts <- seq2nt(seq)
nts
# Output:
# [1] "A" "T" "G" "C" "G" "A" "T" "C"
```


### `nt2seq()`
Convert a vector of nucleotides back into a single sequence string.
```R
# Example
nts <- c("A", "T", "G", "C", "G", "A", "T", "C")

seq <- nt2seq(nts)
seq
# Output:
# [1] "ATGCGATC"
```

### `aln2mx()`
Convert an alignment (a named character vector where all strings have the same length) into a character matrix.  
Rows correspond to sequences, columns correspond to alignment positions.
```R
# Example
aln <- c(
  s1 = "ATG-C",
  s2 = "A-GGC",
  s3 = "ATGGC"
)

mx <- aln2mx(aln)
mx
# Output:
#     [,1] [,2] [,3] [,4] [,5]
# s1 "A"  "T"  "G"  "-"  "C"
# s2 "A"  "-"  "G"  "G"  "C"
# s3 "A"  "T"  "G"  "G"  "C"
```


### `mx2aln()`
Convert a character matrix back into a character-vector alignment.
```R
# Example
aln2 <- mx2aln(mx)
aln2
# Output:
#    s1    s2    s3 
# "ATG-C" "A-GGC" "ATGGC"
```

### `aln2pos()`

Convert an alignment (character vector of equal-length strings) into a matrix of ungapped positions.
```R
# Example
aln <- c(
  s1 = "ATG-C",
  s2 = "A-GGC"
)

pos <- aln2pos(aln)
pos
# Output:
#     [,1] [,2] [,3] [,4] [,5]
# s1    1    2    3    0    4
# s2    1    0    2    3    4
```

### `mx2pos()`
```R
# Example
mx
#     [,1] [,2] [,3] [,4] [,5]
# s1 "A"  "T"  "G"  "-"  "C"
# s2 "A"  "-"  "G"  "G"  "C"

pos <- mx2pos(mx)
pos
# Output:
#     [,1] [,2] [,3] [,4] [,5]
# s1    1    2    3    0    4
# s2    1    0    2    3    4
```






