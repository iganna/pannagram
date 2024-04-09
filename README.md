<p align="left">
<img src="https://github.com/iganna/pannagram/blob/main/examples/pannagram_logo.png" width="30%" height="auto">
</p>

# Pannagram


Pannagram is a toolkit for analyzing pan-genomes, structural variants, and annotations. 
Additionally, Pannagram contains useful functions for visualizing sequences and alignments.


## Recreating working environment

Make sure you have [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) installed. To create and activate the package environment run:
```sh
conda env create -f pannagram_conda_env.yml
conda activate pannagram_conda_env
```


## Pangenome linear genome alignment

Pangenome alignment can be built in two modes:
1. **reference-free**:
```
./pangen.sh -path_in 'input_folder_with_all_genomes'  \
			-path_out 'output_folder' \
			-nchr_query 5 -nchr_ref 5 
```

2. **reference-based**:
```
./pangen_ref.sh  -ref 'tari10'  
                 -path_in 'input_folder_with_all_genomes'  \
			     -path_out 'output_folder' \
			     -nchr_query 5 -nchr_ref 5 
```


An extended description of the parameters can be read by executing:
```
./pangen.sh -help
./pangen_ref.sh -help
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


## Dependencies


BiocManager::muscle

foreach
doParallel
optparse
BiocManager::crayon
BiocManager::rhdf5
msa
dplyr
seqinr
foreach
stringr
ggplot2
utils.R сам устанавливает crayon.


-->

## Acknowledgements


Thanks for the testing:  
* Anna Glushkevich  
* Elizaveta Grigoreva  
* Jilong Ma  
* Alexander Bezlepsky

Logo was generated with the help of DALL-E

To acknowledge the utilized process parallelization tool, reference:
O. Tange (2018): GNU Parallel 2018, Mar 2018, ISBN 9781387509881,
  DOI https://doi.org/10.5281/zenodo.1146014
