# Dotplot functions

The Panngram R library includes a set of [dot-plot](https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics)) functions.  
Dot-plots help identify similarity, repeats, inversions, and other structural patterns between two sequences.  
All dot-plot functions use two core parameters:

- `wsize` – window size (length of the sliding window)
- `nmatch` – minimum number of matching positions within the window required to draw a dot  

A dot is drawn when at least nmatch characters match inside a wsize window.  
Both sequences start in the bottom-left of the plot:  
(**→**) moves along sequence 1,   
(**↑**) moves along sequence 2.


The Panngram R library provides 4 handy functions to visualise dotplots:
1. **`dotplot`** – pairwise dot-plot of two nucleotide sequences.
2. **`dotself`** – self-comparison of a nucleotide sequence.
3. **`dotprot`** – pairwise dot-plot for two protein sequences.
4. **`dotgrid`** – a grid of dot-plots across multiple nmatch thresholds.

## Input Sequences

All functions accept:
- Plain strings
```
seq1 <- "ACGTACGTACGT"
```
- Symbol vectors
```
seq2 <- c("A","C","G","T","A","C","G","T")
```

The sequences used in example figures can be found at the end of this page.

## 1. dotplot

Basic pairwise dot plot between two sequences using the default parameters (`wsize = 15`, `nmatch = 12`).
```
dotplot(seq1, seq2)
```

Example with 80% similarity, i.e. 80 matches out of 100:
```
dotplot(seq1, seq2, wsize=10, nmatch=8)
```

The results are shown below:
<div style="width: 40%; display: flex; justify-content: flex-start; gap: 10px;">
  <img src="images/dotplot.png" style="width: 50%; object-fit: cover;">
  <img src="images/dotplot80.png" style="width: 50%; object-fit: cover;">
</div>

### Interpretation

- **Dark dots** – matches on the **forward strand**  
- **Red dots** – matches on the **reverse complement**

In this dotplot, several structural features can be observed:

- Flanking-regions: duplications present in both sequences  
- Around 1 kb: a large inversion
- Around 2 kb: multiple duplications in seq2
- Around 3 kb: a deletion in seq2
- Around 4 kb: a duplication with inversion in seq2 (resulting in a palindromic region)

## 2. dotself

Displays a self-comparison of the sequence.
- Lower dark triangle: forward strand vs itself (duplications).
- Upper red triangle: forward strand vs reverse complement (inversions/palindromes).

Example:
```
dotself(seq2, wsize=15, nmatch=11)
```

The results are shown below:
<div style="width: 20%;">
<p align="left">
  <img src="images/dotplot_self.png" style="width:100%; object-fit:cover;"/>
</p>
</div>


## 3. dotprot

Pairwise dot-plot for protein sequences.
The function works similarly to dotplot, but Only the forward orientation of the sequences is compared and reverse-complement matches are ignored.

Default parameters: `wsize = 10`, `nmatch = 5` (~50% similarity).
```
dotprot(prot1, prot2)
```

Examples with 80% and 100% similarity:
```
dotprot(prot1, prot2, wsize = 10, nmatch = 8)
dotprot(prot1, prot2, wsize = 10, nmatch = 10)
```

The results are shown below:
<div style="width: 50%; display: flex; justify-content: flex-start; gap: 10px;">
  <img src="images/dotprot.png" style="width: 30%; object-fit: cover;">
  <img src="images/dotprot80.png" style="width: 30%; object-fit: cover;">
  <img src="images/dotprot100.png" style="width: 30%; object-fit: cover;">
</div>

## 4. dotgrid

This function generates multiple dot-plots at once using a vector of match thresholds.
This is useful when the optimal threshold is unclear — for noisy, diverged, or highly repetitive sequences.

```R
dotgrid(seq1, seq2, 15, c(15, 12, 10, 8), nrow = 2)
```

The output:
<div style="width: 40%;">
<p align="left">
  <img src="images/dotgrid.png" style="width:100%; object-fit:cover;"/>
</p>
</div>


## Examples used above

### Nucleotide sequences
**[View](../examples/seqs_dotplot.txt)** · **[Download](../examples/seqs_dotplot.fasta)**  
This example contains a *Arabidopsis thaliana* transposable element and a sequence with artificial modifications.  
To load these sequences in R, run:

```r
library(pannagram)

seqs <- readFasta("seqs_dotplot.fasta")
seq1 <- seqs[1]
seq2 <- seqs[2]
```

### Protein sequences
**[View](../examples/seqs_dotprot.txt)** · **[Download](../examples/seqs_dotprot.fasta)**  
This example compares a coding region from an *Arabidopsis thaliana* transposable element and a hypothetical protein from *Cucumis melo*.  
To load these sequences in R, run:

```r
library(pannagram)

prots <- readFasta("seqs_dotprot.fasta")
prot1 <- prots[1]
prot2 <- prots[2]
```







