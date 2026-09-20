# SNP Calling

The function `getSNPs()` calls SNPs from the pangenome alignment, the same way as `features -snp`.  
It requires consensus sequences produced by `features -seq`.

```R
library(pannagram)

path.project <- "path_to_alignment_project/"

# SNPs in pangenome coordinates
getSNPs(path.proj = path.project)

# SNPs in coordinates of a genome
getSNPs(path.proj = path.project, acc = "name_genome")

# Input for SINGER, from the trustable syntenic positions only
getSNPs(path.proj = path.project, aln.type = "synteny", singer = TRUE)
```

Parameters:

- **path.proj** – the path to the project directory
- **acc** – the genome in whose coordinates positions are reported (Default: `"pangenome"`)
- **positions** – which positions to save: `"snp"` (Default), `"invariant"` or `"all"`
- **ref.acc** ⚠️ – the name of the reference genome, **only** if the alignment is reference-based
- **aln.type** – type of the alignment: `"pan"` (Default) or, e.g., `"synteny"` for the trustable syntenic positions. For non-default types the consensus files and the output VCF get the suffix `_<aln.type>`
- **cores** – number of cores (Default: 1)
- **singer** – if `TRUE`, the VCF is prepared for [SINGER](https://github.com/popgenmethods/SINGER): only biallelic positions without gaps, haploid genotypes `0`/`1`; run SINGER with `-ploidy 1` (Default: `FALSE`). Files: `singer_*.vcf`

Output `VCF` files, one per chromosome, are saved to `${PATH_PROJECT}/features/snp/`:  
`snps_*_pangen.vcf` or `snps_*_<acc>.vcf` (prefix `invariant_` or `all_` for other `positions`).
