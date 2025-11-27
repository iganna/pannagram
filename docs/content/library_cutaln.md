# Extracting a Specific Region from a Pangenome Alignment
The Pannagram R library provides a functionality to extract a specific region  
from a pangenome alignment by choosing a chromosome and accession, and then providing the region’s coordinates.

Before extraction, specify the chromosome index, accession name, region start and end positions,  
and the path to the alignment project (the same value used for the `-path_project` option in `pannagram`):

```R
library(pannagram)

# Chromosome index
i.chr <- 1

# Accession (genome) name
acc <- "name_genome"

# Start and end positions of the region
pos.beg <- 10000
pos.end <- 20000

# Path to the pangenome alignment project
path.project <- "path_to_alignment_project/"
```

## Get the region sequences
To extract the multiple-sequence alignment for a selected genomic region:
```
aln.seq <- getRegion(
  i.chr     = i.chr,
  acc       = acc,
  p.beg     = pos.beg,
  p.end     = pos.end,
  path.proj = path.project,
  mode      = "seq"          # ← type of output to return
)
```

The output is a character matrix where:
- Rows represent genomes (accessions)
- Columns represent alignment positions within the selected region

To visualize the extracted alignment:

```
msaplot(aln.seq)
```

## Get the region positions
If instead of sequences you need a matrix of corresponsing genomic positions, use `pos` mode:
```
aln.seq <- getRegion(
  i.chr     = i.chr,
  acc       = acc,
  p.beg     = pos.beg,
  p.end     = pos.end,
  path.proj = path.project,
  mode      = "pos"          # ← type of output to return
)

```

## ...for Reference-Based Alignments

For reference-based alignments, specify the reference name: `ref='name_of_reference'`.

```
aln.seq <- getRegion(
  i.chr     = i.chr,
  acc       = acc,
  p.beg     = pos.beg,
  p.end     = pos.end,
  path.proj = path.project,
  mode      = "seq",
  ref       = ref          # ← reference genome
)
```






