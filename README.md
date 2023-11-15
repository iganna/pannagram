# pannagram

Pannagram takes genomes as input, some of which it identifies as reference genomes. It then performs multiple genome alignments in several stages. First, it aligns the genomes to all specified references, identifies blocks that are identical across all alignments, and then uses multiple alignments to fill in the gaps between blocks.

The pipeline for aligning genomes to the reference genome can be found in `pipeline.sh`.


To acknowledge the utilized process parallelization tool, reference:
O. Tange (2018): GNU Parallel 2018, Mar 2018, ISBN 9781387509881,
  DOI https://doi.org/10.5281/zenodo.1146014
