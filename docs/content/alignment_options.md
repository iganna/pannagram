# Extra parameters for Alignment

### Number of chromosomes

If the maximum number of chromosomes is not specified, the program scans genomes in `PATH_DATA` to determine it automatically.  
If all query genomes have the same number of chromosomes, this value is used as `N`.  
Otherwise, an error will be raised.
If genomes contain different numbers of chromosomes, specify it manually using `-nchr N`.

In **reference-based** alignment modes, you can set a separate number of chromosomes for reference genomes using `-nchr_ref M`.  
In this case, output file names will include all combinations of the form `N_M`.
All reference genomes must have the same number of chromosomes (`M`), and all query genomes too (`N`).  

## Reference Genome Location

All query accessions to be aligned must be located in the `PATH_GENOMES` folder.  
However, reference genomes might be stored in a different folder, `PATH_REF`.  
In that case, the argument `-path_ref ${PATH_REF}` should be provided.

> ⚠️ **Warning:**  
> All reference genomes must be located **either** in `PATH_GENOMES` **or** in `PATH_REF`,  
> but **not in both** simultaneously.

## Select Genomes for Alignment

If you want to analyze only a subset of genomes, you can provide a text file containing the desired genome names (one per line): `-accessions <FILE_ACCESSIONS>`.

## Chromosome-to-Chromosome Combinations

There are two flags for comparing chromosomes between each other:

- `-one2one` — when all chromosomes in all genomes are sorted in the same order and should be aligned one-to-one.  
- `-all2all` — when all chromosomes should be compared with each other.

**Preliminary alignment mode** uses the `-all2all` flag by default, as this mode is needed to check the quality of genomes. The `-one2one` flag is also allowed.  
**Reference-based alignment mode** uses the `-one2one` flag by default. The `-all2all` flag is also allowed.  
**Reference-free alignment mode** uses the `-one2one` flag by default. The `-all2all` flag is **not allowed**.

## Provide Combinations of chromosomes


If you want to consider only specific combinations of chromosomes, you can provide them in a separate file with parameter `-combinations <FILE_COMB>`.

The file should contain **two columns**:

- First column: Chromosome number of the query accession  
- Second column: Chromosome number of the reference genome



