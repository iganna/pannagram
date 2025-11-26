# Sequence Alignment plots

The Pannagram R library have two functions to visualise Multiple Sequence Alignments.
The input for the alignemnts should be provided as a matrix of

## 1. msaplot

```
msaplot(aln)
```

The results are shown below:
<div style="width: 40%;">
<p align="left">
  <img src="images/msaplot.png" style="width:100%; object-fit:cover;"/>
</p>
</div>


Other parameters:

msaplot <- function(seqs.mx, seq.type='nt', msa.cols = NULL, show.legend=F){

## 2. msadiff

Very useful when one wants to see the differences in the alignment rather than nucleotides.

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



## Alignments used above

**[View](../examples/seqs_msaplot.txt)** · **[Download](../examples/seqs_msaplot.fasta)**  
This example contains alignment of *Arabidopsis thaliana* transposable elements of DNA/HAT Superfamily, family ATHATN2.  
To load these sequences in R, run:

```r
library(pannagram)

aln <- readFasta("seqs_msaplot.fasta")
```


### Bonus: Additional analysis

It's always interesting to observe the structure of transposon sequence.  
И текущий пример - не исключение:

```
dotgrid(aln[1], aln[1], 
	15, c(12, 10, 8))
```

<div style="width: 60%;">
<p align="left">
  <img src="images/msaplot_dot.png" style="width:100%; object-fit:cover;"/>
</p>
</div>


