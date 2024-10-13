# Useful options

## Format of the output files

In the default case, the names of the output files contain combinations of `N_N`, where `N` represents the chromosome numbers that are aligned together. 
This implies that all chromosomes in genomes are sorted in the corresponding order, and all first chromosomes, all second chromosomes, and so on are aligned separately.


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


