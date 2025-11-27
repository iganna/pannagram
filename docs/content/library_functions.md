# Functions with sequences


## Save Figures

savePDF
savePNG

## Read and Write Fasta

```
seqs = readFasta(<>)
writeFasta(<>)
```

## Nucleotide Content

### Plot nt density in windows
```
ntplot(s)
```

<div style="width: 40%;">
<p align="left">
  <img src="images/ntplot.png" style="width:100%; object-fit:cover;"/>
</p>
</div>


```
ntplot(seq, nt.separate = T)
```

<div style="width: 40%;">
<p align="left">
  <img src="images/ntplot_separate.png" style="width:100%; object-fit:cover;"/>
</p>
</div>


### Consensus
mx2cons
mx2profile


### GC
gcContent <- function(s, wnd = 1000)

## Converters between types

### String and a vector of characters
seq2nt()
nt2seq()


aln2mx
mx2aln


aln2pos
mx2pos




## Reverse Complement

degenerated nucleotides are suppotred
revCompl
revComplSeq