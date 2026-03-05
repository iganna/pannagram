# Building Alignments in Different Modes

The alignment can be built in three modes:

- **`I` Preliminary** – builds only a skeleton alignment to visually check the correspondence between aligned genomes.  
- **`II` Reference-based** – aligns all genomes to a chosen reference genome.  
- **`III` Reference-free** – creates a pangenome alignment independent of any reference genome.


Before running the alignment, define the key working paths in the command line:
- `PATH_GENOMES` - the path to the data folder containing files with whole genomes.
- `PATH_PROJECT` - the path to the output directory


> ⚠️ **Warning:**  
> Ensure that `PATH_PROJECT` is set to a completely new folder. Existing files in this directory may be overwritten.


## **`I`** Preliminary Mode

The **preliminary mode does not** perform a full pangenome reconstruction.  
It creates reference-based skeletion alignments for a quick visual check of your input data.  
This helps you verify chromosome correspondence before running a full analysis.

This mode **requires** the name of the genome that will be used as the **reference**.

Run Pannagram in preliminary mode:
```bash
# Define the reference genome name
REF_NAME="RefGenomeName"

# Run Pannagram in preliminary mode
pannagram -pre \
          -path_in ${PATH_DATA} \
          -path_out ${PATH_PROJECT} \
          -ref ${REF_NAME} \
          -cores 8

```

As a result, you will get a dot plot visualization of the reference-based alignments in PDF format.
It can be found at the following path:

```
${PATH_PROJECT}/
├── features/
└── plots/
    └── synteny_pairwise/
        └── ${REF_NAME}/   ← The resulting dot plots are here.
```

On the resulting plots, the chromosomes of the reference genome are shown on the Vertical axis, while those of the aligned accession are shown on the Horizontal axis.  
Chromosomes are separated by vertical and horizontal lines.  
Below are examples of such dot plots: one for the *Arabidopsis lyrata* genome aligned to *Arabidopsis thaliana*, and another showing a comparison of two *Euwallacea fornicatus* genomes.

<div style="width:50%; display:flex; justify-content:center; align-items:center;">
  <div style="width:50%; text-align:center;">
    <img src="images/synteny_pw_0-11B21.png" style="width:100%; object-fit:cover;"/>
    <div><sub><i>A. thaliana</i> vs <i>A. lyrata</i></sub></div>
  </div>
  <div style="width:50%; text-align:center;">
    <img src="images/synteny_pw_insect_euwallacea_fornicatus.png" style="width:100%; object-fit:cover;"/>
    <div><sub>Two genomes of <i>Euwallacea fornicatus</i></sub></div>
  </div>
</div>


> 💡 **Tips:**  
> 1. If the chromosomes are not sorted but there is a one-to-one correspondence, it is recommended to reorder the original genome-files so that the chromosomes are arranged in the same order. This can be done automatically using the `chromotools` module.
> 2. If some chromosomes are in the reverse complement orientation to the reference, it is recommended to reorient them for clearer analysis later on, but this is **not required**. This can be done automatically using the `chromotools` module.
> 3. If certain accessions appear suspicious based on visual inspection, and you decide not to analyze them, create a file `FILE_ACC` that contains only the accessions you wish to analyze, with each accession listed on a separate line. Then, for all further analyses, use the flag `-accessions ${FILE_ACC}` in every command.

> ⚠️ **Warning:**  
> If you have modified the genome files after the preliminary mode, please **change** the `PATH_PROJECT` variable before building the alignment on new modified genomes.

TODO:
<!-- In alignment modes, it is expected that the number of chromosomes in the query accessions is the same. 
If they differ — for instance, if the genome files contain not only chromosomes but also a varying number of scaffolds — you need to use the `-nchr N` flag, where `N` represents the number first sequences in the genome FASTA files.
 -->

## **`II`** Reference-Based mode

This mode produces the alignmnet of all genomes to the reference genome:
```shell
pannagram \
    -path_in '${PATH_DATA}' \
    -path_out '${PATH_PROJECT}' \
    -ref '<reference genome filename with no FASTA suffix>' \
    -cores 8
```

The result files are located here: 
```
${PATH_PROJECT}/
├── features/
|   └── alignments/
|       └── ref_*.h5      ← Resultant alignments are here
└── plots/
```

These HDF5-format files contain matrices of corresponding positions between genomes.  
If you want to extract alignment features (e.g., SNPs) after running the reference-based mode, please use the flag `-anl_type ref` when executing the `features` module.


## **`III`** Reference-Free mode

This mode does **not require a reference genome**, you can simply run it as follows:

```shell
pannagram \
    -path_in ${PATH_DATA}\
    -path_out ${PATH_PROJECT} \
    -cores 8
```



In Reference-Free Mode, Pannagram eliminates the need for a reference genome by performing internal randomization among several genomes used as temporary references.  
By default, it uses two randomly selected genomes as references. You can change this number with the parameter `-nref X`, where `X` specifies how many genomes to use for randomization.


If you already know the structure of your dataset, we recommend providing your own reference genomes for randomization.  
Choose several genomes that are relatively genetically distant within your dataset, and specify them using the `-refs` parameter without spaces between names, for example: `-refs 'Name1,Name2,Name3'`.


The output files are located here: 
```
${PATH_PROJECT}/
├── features/
|   └── alignments/
|       └── pan_*.h5      ← Resultant alignments are here
└── plots/
```

Each .h5 file (HDF5-format) contains matrices of corresponding positions between genomes.

If you want to extract alignment features after running the Reference-Free Mode, you **do not need to specify** the `-anl_type pan` flag when executing the **features** module.
However, if you wish, you can explicitly set it.


## Names of output HDF5-files

File names start with the prefix `ref` contains also two numbers in the format `_N_M_`, where:
- `N` is the chromosome number of the query genome 
- `M` is the chromosome number of the reference genome.

For `pan` files, the numbers are the same (`_N_N_`), indicating that all Nth chromosomes are aligned together in the reference-free manner.


