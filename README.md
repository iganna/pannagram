# Pannagram

Pannagram takes genomes as input, some of which it identifies as reference genomes. It then performs multiple genome alignments in several stages. First, it aligns the genomes to all specified references, identifies blocks that are identical across all alignments, and then uses multiple alignments to fill in the gaps between blocks.

The pipeline for aligning genomes to the reference genome can be found in `pipeline.sh`.


## Parameters to run `pipeline.sh`


## Script Overview


- `-pref_global <value>`
   The global prefix for folders that will be created during the pipeline.

- `-ref_pref <value>`
   The basename of the reference genome file (without extension), that the pipeline will use as a reference for alignment.

- `-path_chr_ref <value>`
   The directory containing the reference genome.

- `-n_chr_ref <value>`
   The number of chromosomes in the reference genome.

- `-path_in <value>`
   The directory containing query genomes that must be aligned to the reference.

- `-n_chr_query <value>`
   The number of chromosomes in the query genomes.

- `-path_parts <value>`
   (Optional) The directory containing chunked chromosomes.

- `-path_chr_acc <value>`
   (Optional) The directory containing separate chromosomes.

- `-path_consensus <value>`
   The directory where the common consensus results will be stored.

- `-sort_chr_len <value>`
   (Optional) Flag to sort chromosomes in the genomes by length.

- `-part_len <value>`
   (Optional) The length of chunks used in the analysis.
   Default = 500bp

- `-all_cmp <value>`
   Flag to align all chromosomes to each other, not just the corresponding ones.

- `-p_ident <value>`
   The identity value used for the internal BLAST comparisons.

- `-acc_anal <value>`
   The file containing the set of accessions to analyze.

- `-cores <value>`
   The number of CPU cores to use for running the pipeline.


## Acknowledgements

To acknowledge the utilized process parallelization tool, reference:
O. Tange (2018): GNU Parallel 2018, Mar 2018, ISBN 9781387509881,
  DOI https://doi.org/10.5281/zenodo.1146014
