<p align="left">
<img src="https://github.com/iganna/pannagram/blob/main/examples/pannagram_logo.png" width="30%" height="auto">
</p>

# Pannagram


Pannagram is a toolkit for analyzing pan-genomes, structural variants, and annotations. Additionally, Pannagram contains useful functions for visualizing sequences and alignments.

## Multiple genome alignment

Multiple genome alignment can be performed on a single reference genome or using multiple references to create a reference-free alignment. These modes are implemented in `pangen_ref.sh` and `pangen_consensus.sh`, respectively.

The quickstart toy-examples are:
```
./pangen_ref.sh -pref_global 'output_folder' -ref_pref '0' \
    -path_in 'thaliana_genomes' -n_chr_query 5 -n_chr_ref 5
./pangen_consensus.sh -pref_global 'output_folder' -ref_set '0,10024' \
    -path_in 'thaliana_genomes' -n_chr_query 5 -n_chr_ref 5
```

An extended description of the parameters can be read by executing:
```
./pangen_ref.sh -help
./pangen_consensus.sh -help
```

## Visualisation
Pannagram contains a number of useful methods for visualization.

* Dotplot
* MSA


## Additional useful tools
### Search for similar sequences...

There are two distinct approaches available for searching for similarities for a **set of sequences** (a fasta file with genes or structural variants). 
The first approach involves searching against entire genomes or individual chromosomes. 
The second approach, in contrast, is designed to search for similarities against another set of sequences. 

#### ...in the genome
The quickstart toy-example is:
```
./sim_in_genome.sh -in genes.fasta -genome genome.fasta -out out.txt
```
The result is a GFF file with hits matching the similarity threshold.


#### ...in another set
The quickstart toy-example is:
```
sim_in_seqs.sh -in genes.fasta -set genome.fasta -out out.txt
```
The result is an RDS (R Data Structure) table. 
This table shows the coverage of one sequence over another and 
includes a flag column that indicates whether the sequences meet the similarity threshold. 
Additionally, the second script takes into account the coverage strand, 
determining not just if a sequence is covered, but also if it's covered in a specific orientation.

<!--

## Parameters to run `pipeline_consensus.sh`


## Script Overview


- `-pref_global <value>`
   The global prefix for folders that will be created during the pipeline.

- `-ref_pref <value>`
   The basename of the reference genome file (without extension), that the pipeline will use as a reference for alignment.

- `-path_chr_ref <value>`
   (Optional) 
   The directory containing **chromosomes** of the reference genome.
   можно -path_chr_ref ваще не указывать, если референс находится там же, где и остальные

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


## Dependencies

install.packages("optparse")
install.packages("seqinr")
install.packages("foreach")
install.packages("doParallel")

stringi and stringr,



## Dependencies



BiocManager::muscle

foreach
doParallel
optparse
BiocManager::crayon
BiocManager::rhdf5
msa
dplyr
BiocManager::Biostrings
seqinr
foreach
stringr
ggplot2
utils.R сам устанавливает crayon.


## TODO:

следать хромосомные файлы для референса, если он не распарсен
прописать дефолты
все рабоает нормально. и работает без -path_chr_ref
-all_clp rкак-то переделать флаг

поменять фон


-->

## Acknowledgements

Thanks for the testing:
Anna Glushkevich
Elizaveta Grigoreva
Jilong Ma

Logo was generated with the help of DALL-E

To acknowledge the utilized process parallelization tool, reference:
O. Tange (2018): GNU Parallel 2018, Mar 2018, ISBN 9781387509881,
  DOI https://doi.org/10.5281/zenodo.1146014
