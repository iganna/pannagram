
## Cut a part of the alignment

if one want to check a special part of the alignment you should defibe the following variables:

```R
library(pannagram)
pos.beg = 10000
pos.end = 20000
acc <- 'name_genome'
```

Then one can cut the alignment in two modes: in sequence (default) or in matrix of posisiotns:
And then execute
```
aln.seq <- cutAln(path.proj=path_to_ref_free_alignment,
	          mode='seq')


aln.pos <- cutAln(path.proj=path_to_ref_free_alignment,
	          mode='pos')
```

The results are the matrix:
Rows -  genomes
Columns - postisions

For seqience matrxi one can visualise it:

```
msaplot(aln.seq)
```

The result will be as follows:







