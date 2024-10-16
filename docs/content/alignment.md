
# Construct the Multiple Genome Alignment

Before running the example, define the key working paths in the command line:
- `PATH_DATA` - the path to the data folder
- `PATH_PAN` - the path to the `pannagram` folder
- `PATH_OUT` - the path to the output directory

> ⚠️ **Warning:**  
> Ensure that `PATH_OUT` is set to a completely new folder. Existing files in this directory may be overwritten.

## Preliminary mode

This mode might be useful for quickly reviewing the data.
As a result, you will get a visualization of draft reference-based alignments.
To set up the reference, please define the `REF_NAME` variable (basename without extension).
And run Pannagram in preliminary mode:
```
conda activate pannagram
${PATH_PAN}pannagram.sh -path_in ${PATH_DATA} -path_out ${PATH_OUT} -cores 8 -ref ${REF_NAME} -pre 
```
<!-- pannagram.sh -path_in ${PATH_DATA} -path_out ${PATH_OUT} -cores 8 -nchr 1 -log 2 -ref ${REF_NAME} -accessions ${FILE_ACC} -->

The result is pairwise dot plots in PDF format:
```
cd "${PATH_OUT}plots/plots_${REF_NAME}/"
```

These plots are useful for deciding in which mode Pannagram should be used:
- `one2one`, when all chromosomes in all genomes are sorted in the same order and should be aligned one-to-one
- `all2all`, when all chromosomes should be compared with each other.

> 💡 **Tips:**  
> 1. If the chromosomes are not sorted but there is a one-to-one correspondence, it is recommended to reorder the original genome-files so that the chromosomes are arranged in the correct order for more convenient further analysis.
> 2. If some chromosomes are in the reverse complement orientation to the reference, it is recommended to reorient them for clearer analysis later on, but this is not strictly required.
> 3. If certain accessions appear suspicious based on visual inspection, and you decide not to analyze them, create a file called `FILE_ACC` that contains only the accessions you wish to analyze, with each accession listed on a separate line. For all further analyses, use the flag `-accessions ${FILE_ACC}` in every command.

> ⚠️ **Warning:**  
> If you have modified the genome-files after the preliminary mode, please **delete** the `PATH_OUT` folder completely before building the alignment.

## Alignment modes

In alignment modes, it is expected that the number of chromosomes in the query accessions is the same. 
If they differ — for instance, if the genome files contain not only chromosomes but also a varying number of scaffolds — you need to use the `-nchr N` flag, where `N` represents the number first sequences in the genome FASTA files.

### Reference-based mode

This mode produces the alignmnet of all genomes to the reference genome:
```
conda activate pannagram
${PATH_PAN}pannagram.sh -path_in ${PATH_DATA} -path_out ${PATH_OUT} -cores 8 -ref ${REF_NAME} 
```

The result files are located at and called `ref_*.h5`. The prefix `ref` is crucial for the `analysis.sh` script.
```
PATH_ALN="${PATH_OUT}intermediate/consensus/"
cd ${PATH_ALN}
ls -lrt ref*
```

### REFERENCE-FREE mode

This mode doesn't require specification of the reference genome.
```
conda activate pannagram
${PATH_PAN}pannagram.sh -path_in ${PATH_DATA} -path_out ${PATH_OUT} -cores 8
```

The result files are located at and called `msa_*.h5`. The prefix `msa` is crucial for the `analysis.sh` script.
```
PATH_ALN="${PATH_OUT}"
cd ${PATH_ALN}
ls -lrt msa*
```

> 💡 **Tips:**  
> In this mode, the randomisation of reference genomes occurs.
> By default, two random reference genomes are chosen, and the first genome is used to sort the positions in the final alignment.
> You can specify the number of references to use with the `-nref` flag or directly provide specific references using the `-refs` flag.

## Format of the output files

In the default case, the names of the output files contain combinations of `N_N`, where `N` represents the chromosome numbers that are aligned together. 
This implies that all chromosomes in genomes are sorted in the corresponding order, and all first chromosomes, all second chromosomes, and so on are aligned separately.