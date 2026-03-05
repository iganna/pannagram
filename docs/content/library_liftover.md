# Coordinate Conversion and Liftover

The Pangenome alignment is a matrix of corresponding positions between accessions (genomes).  
This structure allows for converting coordinates between accessions.

The Pangenome R library provides three functions for performing coordinate conversion and liftover:

- **`gff2gff()`** – converts GFF annotations between accessions 
- **`bed2bed()`** – converts BED intervals between accessions  
- **`pos2pos()`** – converts individual genomic positions between accessions


## ⚠️ When Conversion/Liftover Is Not Possible

The pangenome alignment does not contain all genomic positions from every accession:

1. Homology in centromeric regions (and other highly repetative regions) is often impossible to find,   
and as a result, these regions are missing from the alignment.
2. Another type of missing region occurs near the flanking boundaries of inversions  
the so-called scars.

If a user attempts to convert **regions that overlap or span any of these non-aligned areas**,  
the conversion for those regions cannot be performed.  
Consequently, such regions are **skipped and do not appear in the output**.

## Parameters

All of the functions mentioned above share a common set of parameters:

- **path.proj** – the path to the project directory
- **acc1** – the name of the accession *from which* coordinates are converted  
  (or `"pangenome"` when using pangenome coordinates)
- **acc2** – the name of the accession *to which* coordinates are converted  
  (or `"pangenome"` when converting to pangenome coordinates)
- **s.chr** – chromosome name pattern  
  For example, with `s.chr = "_Chr"` (Default), the expected format is `*_ChrX`, where **X** is the chromosome number.
- **ref.acc** ⚠️ – the name of the reference genome, **only** if the alignment is reference-based
- **exact.match** ⚠️ – This parameter controls behavior when positions of region boundaries **do not exist** in the pangenome alignment.  
  If set to `FALSE` (**Default**), only exact boundary matches are returned.  
  If set to `TRUE`, overlapping aligned segments are returned so the region is not lost.


## 1. gff2gff

This function lifts over genomic annotations in [GFF format](https://en.wikipedia.org/wiki/General_feature_format) from one accession to another.


```
# Load the library
library(pannagram)

# Define variables
path.to.alignment <- 'path_to_alignment_project/'
acc1 <- 'name_genome1'
acc2 <- 'name_genome2'

# Example input files
file.gff1 <- 'file_annotation_genome1.gff'
file.bed1 <- 'file_regions_genome1.bed'

# Read annotation/regions
gff1 <- read.table(file.gff1)

gff2.ref.free <- gff2gff(
    path.proj = path_to_ref_free_alignment,
    acc1 = acc1,
    acc2 = acc2,
    gff1 = gff1
)

```

## 2. bed2bed

This function lifts over intervals in [BED format](https://en.wikipedia.org/wiki/BED_(file_format)) from one accession to another.

```
# Example BED file
file.bed1 <- 'file_regions_genome1.bed'

# Read BED intervals
bed1 <- read.table(file.bed1)

# Liftover BED regions
bed2.ref.free <- bed2bed(
    path.proj = path.to.alignment,
    acc1 = acc1,
    acc2 = acc2,
    bed1 = bed1
)

```

## 3. pos2pos

This function converts individual genomic positions between accessions.  
It accepts a data frame with the following three columns:

- **chrom** — chromosome name  
- **pos** — genomic position  
- **info** — optional annotation or identifier (returned in the output)

```
pos1 <- data.frame(
    chrom = c("Chr1", "Chr3"),
    pos   = c(12345, 98765),
    info  = c("geneA", "markerB")
)

# Convert positions between genomes
pos2 <- pos2pos(
    path.proj = path.to.alignment,
    acc1 = acc1,
    acc2 = acc2,
    pos1 = pos1
)
```

