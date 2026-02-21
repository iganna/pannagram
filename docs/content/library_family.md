# Work with Mobile Element Families (Mobilome Families)

This section describes how to extract families of Mobile Elements from the graph of nestednedd built on Structural Variants (SVs).

## 0. Initial data
To start extracting Mobile Element Families, you need either:
– a FASTA file containing the sequences of structural variants, or
– a file with nested similarity information produced by comparing SVs to each other (for example, 85% similarity and coverage).

Both files are expected to be located in:
```
<path_project>/
└── features/
    └── sv/
        ├── seq_sv_large.fasta
        └── nestedness_sv_large_85_85.txt
```

If the file `nestedness_sv_large_85_85.txt` does not exist—because you haven’t run `features -sv_families`,  
or if you want finer control over similarity and coverage thresholds —  
you can generate a the similar file using the `simsearch` module from the Pannagram package.

Run the following commands:

```bash
cd ${PATH_PROJECT}/features/sv
PATH_SIMSEARCH="simsearch_sv_80_80"
simsearch -query_seqs seq_sv_large.fasta -target_seqs seq_sv_large.fasta -sim 80 -cov 80 -out ${PATH_SIMSEARCH}
```

After completion, the file you will use for downstream processing will be located at:

```
<path_project>/
└── features/
    └── sv/
        └── simsearch_sv_80_80/
           └── seq_sv_large_80_80.txt
```

> **Notes**  
> You do **not** need to use all large SVs to build families.  
> You can instead use a subset of SVs or any other sequences of interest.  
> If so, we sould:  
> 1. Preparing a FASTA file containing the sequences you want to cluster.  
> 2. Running `simsearch` with your chosen similarity and coverage thresholds.  
> 3. Using the resulting similarity file as the starting point for the downstream pipeline.

## 1. Read the Nestedness of sequences

```R
library(pannagram)
library(ggplot2)

file.nestedness <- 'nestedness_sv_large_85_85.txt'  # or the output of simsearch


res.sim = read.table(res.sim.file, header = T)

```

The table has the following columns:
```
           name.query         name.target strand length.query length.target coverage.query coverage.target
1 SVgr_1_id_0004|1349  SVgr_1_id_0169|324      -         1349           324      0.2105263       0.8765432
2 SVgr_1_id_0004|1349 SVgr_1_id_0594|1349      +         1349          1349      0.9933284       0.9933284
3 SVgr_1_id_0004|1349 SVgr_1_id_0764|1349      +         1349          1349      1.0000000       1.0000000
4  SVgr_1_id_0036|119   SVgr_1_id_0089|72      +          119            72      0.5126050       0.8472222

```

It guranteed that for the coverage in the table the similarity is not lower than was set up.
coverage.query  means which fraction name.query is covered by name.target 
coverage.target  means which fraction name.target is covered by name.query


## 2. Get edges and visualisation

```R
library(network)
library(igraph)
library(GGally)

# One can filter sequences by length:
min.len = 200
res.sim = filterCoverageMatrix(res.sim, min.len = min.len)

# Get edges for the graph of nestedness
coverage.cutoff = 80
edges = getGraphFromNestedness(res.sim, coverage.cutoff = coverage.cutoff)

# Visualise the Graph of nestedness
p <- ggnet2(g, label = F, edge.color = "black",
            node.size = 1,
            color = '#468B97',
            arrow.gap = 0.01, arrow.size = 2,
)
p

```

## 3. Clean up the graph to extract the Families

The resultant graph may not have separate Families of Mobile Elements. Therefore, some edges or nodes should be removed.



