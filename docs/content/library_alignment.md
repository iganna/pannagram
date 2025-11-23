# Pangenome Functions

There two main Pangenome functional abilities: liftover annotation from one genome for another and cut a part of the alignment.

## Visualisation of synteny: dotplot + pangenome plot ?

## Annotation liftover

If one have a gff file with annotation for one genome, and want to get the same annotation but in coordinates of another genome, `gff2gff` surves it:

```
# load the library
library(pannagram)

# Define variables
acc1 <- 'name_genome1'
acc2 <- 'name_genome2'
file.gff1 <- 'file_annotation_genome1.gff'
gff1 <- read.table(file.gff1)
```

If you want to liftover annotation through the reference-free alignment run this:
```
gff2.ref.free <- gff2gff(path.proj=path_to_ref_free_alignment,
	            acc1=acc1,
	            acc2=acc2,
	            gff1=gff1)

```

If you want to liftover annotation through the reference-based alignment run this:
```
gff2.ref.based <- gff2gff(path.proj=path_to_ref_based_alignment,
	            acc1=acc1,
	            acc2=acc2,
	            gff1=gff1,
	            ref='name_reference_genome',
	            aln.type = 'ref')

```

The same behaviour is provided by functions:
- bed2bed:
- pos2pos:

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







