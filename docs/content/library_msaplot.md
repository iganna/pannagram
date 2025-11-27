# Multiple Sequence Alignment Plots

The Pannagram R library provides two functions for visualising Multiple Sequence Alignments (MSAs):

1. `msaplot` – the alignment, colored by nucleotides or amino acids.
2. `msadiff` – the alignment, highlighting differences between sequence.

### Input anignment
Both functions expect as input either an array with aligned sequences or a matrix of aligned sequences.  
Examples:
```
# Aligned sequences 
aln <- c("ATAG", "CATG", "GCAT")

# Matrix of aligned sequences
aln <-  matrix(
  c("A","C","G","T",
    "A","C","A","T",
    "A","G","G","T"),
  nrow = 3
)
```

> Note:  
> Make sure all sequences are the same length and already aligned  
> (e.g. with an external MSA tool) before calling msaplot() or msadiff().


## 1. msaplot

The function shows the alignment, coloring each cell by nucleotide or amino acid or a custom colormap.  
- Rows correspond to sequences.
- Columns correspond to alignment positions.

```
msaplot(aln)
```

The result is shown below:
<div style="width: 40%;">
<p align="left">
  <img src="images/msaplot.png" style="width:100%; object-fit:cover;"/>
</p>
</div>

### Arguments
- `aln`: A multiple sequence alignment object.
- `seq.type`: Type of sequences in the alignment. Use "nt" for nucleotide sequences (default) or "aa" for amino acid sequences.
- `msa.cols`: Optional custom color palette for residues. Default: the colormap for "nt" or "aa" respectively.
- `show.legend`: Logical (default FALSE); whether to display a legend for the color scheme.


## 2. msadiff

This function is useful for visualizing differences in a multiple sequence alignment relative to a selected reference sequence.  
By default, the first sequence is used as the reference.

```
msadiff(aln)
```

The results are shown below:
<div style="width: 40%;">
<p align="left">
  <img src="images/msadiff.png" style="width:100%; object-fit:cover;"/>
</p>
</div>

```
msadiff(aln, i.ref = 3, show.legend=T)
```

The results are shown below:
<div style="width: 40%;">
<p align="left">
  <img src="images/msadiff_ref3.png" style="width:100%; object-fit:cover;"/>
</p>
</div>


### Arguments
- `aln`: A multiple sequence alignment object.
- `i.ref`: Index of the reference sequence used for comparison.
- `show.legend`: Logical (default FALSE); whether to display a legend for the color scheme.


## Alignment Used in Examples

**[View](../examples/seqs_msaplot.txt)** · **[Download](../examples/seqs_msaplot.fasta)**  
This example contains an alignment of *Arabidopsis thaliana* transposable elements from the *DNA/HAT* superfamily, *ATHATN2* family.
To load these sequences in R, run:
```r
library(pannagram)

aln <- readFasta("seqs_msaplot.fasta")
```

### Bonus Analysis for Inspiration

It is always interesting to examine the structural features of transposon sequences,  
and the current example is no exception:

```
dotgrid(aln[1], aln[1], 
	15, c(12, 10, 8))
```

<div style="width: 60%;">
<p align="left">
  <img src="images/msaplot_dot.png" style="width:100%; object-fit:cover;"/>
</p>
</div>


