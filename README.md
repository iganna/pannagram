<p align="left">
<img src="https://github.com/iganna/pannagram/blob/dev/images/pannagram_logo.png" width="30%" height="auto">
</p>

# Pannagram


Pannagram is a package for constructing pan-genome alignments, analyzing structural variants, and translating annotations between genomes.
Additionally, Pannagram contains useful functions for visualization.


### Recreating working environment

Make sure you have [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) installed. To create and activate the package environment run:
```sh
conda env create -f pannagram_conda_env.yml
conda activate pannagram_conda_env
```
The environment downloads required R libraries, [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/) and [MAFFT](https://mafft.cbrc.jp/alignment/software/manual/manual.html).

## Pangenome linear alignment

### Building the alignment
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


An extended description of the parameters can be read by executing scripts with the flag `-help`.

### Extract information from the pangenome alignment

```
./analys.sh -path_out 'output_folder' 
			-blocks  \	# Synteny block inforamtion for visualisation
			-seq  \		# Consensus sequence of the pangenome
			-snp		# SNP calling
```

### Calling structural variants

When the pangenome linear alignment is built, SVs can be called with the following:
```
./sv.sh -path_out 'output_folder' \
        -gff  \ 					# Output Gff files
        -te -te_file te.fasta  \ 	# Compare with known TE sequences
        -graph  					# Construct the graph of SVs
```


## Visualisation
Pannagram contains a number of useful methods for visualization.

### Visualisation of the pangenome alignment
All genomes together:
<p align="left">
<img src="https://github.com/iganna/pannagram/blob/dev/images/pangenome_alignment.png" width="50%" height="auto">
</p>

A dotplot for a pair of genomes:
<p align="left">
<img src="https://github.com/iganna/pannagram/blob/dev/images/syntenyplot.png" width="50%" height="auto">
</p>

### Seuqnce plot for a fragment of the alignment
In the ACTG-model:
<p align="left">
<img src="https://github.com/iganna/pannagram/blob/dev/images/msaplot.png" width="50%" height="auto">
</p>

In the Polymorphism mode:

<p align="left">
<img src="https://github.com/iganna/pannagram/blob/dev/images/msaplot_diff.png" width="50%" height="auto">
</p>

### Dotplots of sequences
Simultaneously in forward (dark color) and reverce comlement (pink color) strands:
<p align="left">
<img src="https://github.com/iganna/pannagram/blob/dev/images/dotplot.png" width="50%" height="auto">
</p>


### ORF-finder and visualisation


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
