# Auxiliary options

## Number of chromosomes
By default, if the maximum number of chromosomes is not specified (i.e., `-nchr N` is not provided), the genomes of the query accessions in `PATH_DATA` are scanned to determine the number of chromosomes. If all genomes have the same number of chromosomes, `N` is set to this number. Otherwise, an error will occur, and it's recommended to provide `-nchr N`.

It's also possible to specify different number of chromosomes for reference genomes `-nchr_ref M`,
then names of the output files contain all combinations combilations of type `N_M`.

> ⚠️ **Warning:**  
> All reference genomes are expected to have the same number of chromosomes, and all query accessions too. These numbers can differ between them.


## Combinations of chromosomes

If one wants to consider only specific combinations, they can be provided in a separate file in a two-column format:
- First column: chromosome number of the query accession
- Second column: chromosome number of the reference genome.
The falg to use specific combination is `-combinations <FILE_COMB>`.


## Reference genomes location

All query accessions that should be aligned must be located in the `PATH_DATA` folder. However, it might be the case that the reference genomes are located in another folder, `PATH_REF`. In that case, the flag `-path_ref ${PATH_REF}` should be provided. 

> ⚠️ **Warning:**  
> All reference genomes must be located either in `PATH_DATA` or in `PATH_REF`, but not in both simultaneously.


## Chromosome reordering to fit the reference genome


If after a preliminary analysis you observe that the chromosomes in the genomes do not match the expected order relative to the reference genome.


### Manual reodering
The following R functions can be useful for this task:
``` R
library(pannagram)

# Read the FASTA file
genome <- readFasta("genome.fasta")

# Get the reverse complement of the first sequence (chromosome)
chr_rev_compl <- revComplSeq(genome[1])

# Write the reverse complement to a new FASTA file
writeFasta(chr_rev_compl, "chromosome_out.fasta", seq.names = "Chr1_rev_compl")

```
- `readFasta`: Reads a FASTA file and stores sequences.
- `revComplSeq`: Generates the reverse complement of a sequence.
- `writeFasta`: Writes the modified sequence back to a FASTA file.

### Automatic reodering

You can also try the automatic procedure in bash.
First, set the necessary environment variables:
```bash
# Set the reference genome name, output folder, and input genomes path
REF_NAME="reference_genome_name"
PATH_OUT="path_to_output/"
PATH_DATA="path_to_genomes/"
```

Then, create the intermediate paths and run the following R scripts:
```
# Intermediate directories for alignments and resorted genomes
PATH_INTER="${PATH_OUT}intermediate/"
PATH_ALN="${PATH_INTER}alignments_${REF_NAME}/"
PATH_RESORT="${PATH_INTER}resort_${REF_NAME}/"

# Run the R scripts to find the best alignment and rearrange genomes
Rscript ${PATH_PAN}inst/pangen/resort_01_find_best.R --path.aln ${PATH_ALN} --ref ${REF_NAME} --path.resort ${PATH_RESORT}
Rscript ${PATH_PAN}inst/pangen/resort_02_rearrange.R --path.genomes ${PATH_DATA} --ref ${REF_NAME} --path.resort ${PATH_RESORT}
```





