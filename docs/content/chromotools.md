# Chromotools

Chromotools is a module designed to standardize the chromosomal structure across the input set of genomes.  
It rearranges genomes so that all of them follow the same chromosome order and orientation as a chosen reference genome.  


This step is important when the genomes in your collection do not have a one-to-one chromosomal correspondence,   
for example, when some chromosomes are split, merged, or inverted compared to the reference.  
Such inconsistencies are typically identified after visual inspection of the results obtained from the Preliminary mode of Pannagram alignment.  


Once these differences are detected, Chromotools can automatically rearrange the genomes to match the reference chromosomal organization.

Usage:
```bash
chromotools \
  -genomes_in <PATH_GENOMES> \
  -path_project <PATH_PROJECT> \
  -genomes_out <PATH_GENOMES_REARRANGED> \
  -rearrange
```

With the use of Chromotools the *Arabidopsis lyrata* chromosomal organization were converted to the *Arabidopsis thaliana* organization as demonstrated below:
