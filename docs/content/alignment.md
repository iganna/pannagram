
# Construct the Multiple Genome Alignment

Before running the example, define the key working paths in the command line:
- `PATH_DATA` - the path to the data folder
- `PATH_PAN` - the path to the `pannagram` folder
- `PATH_OUT` - the path to the output directory

> ⚠️ **Warning:**  
> Ensure that `PATH_OUT` is set to a completely new folder. Existing files in this directory may be overwritten.

## Test Dataset

The test dataset comprises bacterial genomes (Escherichia coli) from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=562). You can download some genomes manually and put fasta-files into the `PATH_DATA` folder or utilize a tool for automatic downloading.

### Automated Option (Recommended):

1. Create the list of genome IDs as follows:

```bash
echo -e "GCA_000005845.2\nGCA_000008865.2\nGCA_042692195.1\nGCA_042259155.1\nGCA_042189615.1" > "${PATH_DATA}ecoli.txt"
echo -e "GCA_042192135.1\nGCA_042192365.1\nGCA_042161425.1\nGCA_042138005.1\nGCA_964264615.1" >> "${PATH_DATA}ecoli.txt"
echo -e "GCA_042016495.1\nGCA_042017145.1\nGCA_042017895.1\nGCA_040964775.2\nGCA_041950585.1" >> "${PATH_DATA}ecoli.txt"
echo -e "GCA_041954705.1\nGCA_041927255.1\nGCA_041897545.1\nGCA_002853715.1\nGCA_000013265.1" >> "${PATH_DATA}ecoli.txt"
echo -e "GCA_003697165.2\nGCA_000210475.1\nGCA_000010385.1\nGCA_000019645.1\nGCA_003018455.1" >> "${PATH_DATA}ecoli.txt"
echo -e "GCA_024300685.1\nGCA_027925825.1\nGCA_040571355.1\nGCA_003018035.1\nGCA_027925745.1" >> "${PATH_DATA}ecoli.txt"
echo -e "GCA_027925805.1\nGCA_027925765.1\nGCA_027925785.1\nGCA_027925565.1\nGCA_027925845.1" >> "${PATH_DATA}ecoli.txt"
```

2. Specify the folder `PATH_TOOLS` where the download helper [poputils](https://github.com/iganna/poputils) will be cloned and clone the repo:
```bash
cd ${PATH_TOOLS}
git clone https://github.com/iganna/poputils.git
```

3. Download the genomes from NCBI to the folder `PATH_DATA`:
```bash
cd ${PATH_TOOLS}poputils/genomes
./genbank_download_list.sh -f ${PATH_DATA}ecoli.txt -p ${PATH_DATA}
```

> 💡 **Tip:**  
> Ensure that all accessions are approximately of the same size and verify that you haven't downloaded any accessions with suspicious lengths and number of chromosomes. The following command gives you three columns: the name of the file, the size, and the number of chromosomes.
> ```
> for file in "${PATH_DATA}"/*.fasta; do echo "$(basename "$file") $(du -h "$file" | cut -f1) $(grep -c "^>" "$file")"; done
> ```

## Preliminary mode

This mode might be useful for quickly reviewing the data.
As a result, you will get a visualization of draft reference-based alignments.
To set up the reference, please define the `REF_NAME` variable (basename without extension).
For example:
```
REF_NAME=GCA_000005845.2
```

And run Pannagram in preliminary mode:
```
conda activate pannagram
${PATH_PAN}pannagram.sh -path_in ${PATH_DATA} -path_out ${PATH_OUT} -ref ${REF_NAME} -cores 8 -pre 
```

<!-- pannagram.sh -path_in ${PATH_DATA} -path_out ${PATH_OUT} -cores 8 -nchr 1 -log 2 -ref ${REF_NAME} -accessions ${ACC_FILE} -->


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
> 3. If certain accessions appear suspicious based on visual inspection, and you decide not to analyze them, create a file called `ACC_FILE` that contains only the accessions you wish to analyze, with each accession listed on a separate line. For all further analyses, use the flag `-accessions ${ACC_FILE}` in every command. For example:
> ```
> ACC_FILE="${PATH_DATA}ecoli_good.txt"
> echo -e "GCA_000005845.2\nGCA_000008865.2\nGCA_042692195.1\nGCA_042259155.1\nGCA_042189615.1" > "${ACC_FILE}"
> ```

> ⚠️ **Warning:**  
> If you have modified the genome-files after the preliminary mode, please **delete** the `PATH_OUT` folder completely before building the alignment.

## Alignment modes

In alignment modes, it is expected that the number of chromosomes in the query accessions is the same. 
If they differ — for instance, if the genome files contain not only chromosomes but also a varying number of scaffolds — you need to use the `-nchr N` flag, where `N` represents the number first sequences in the genome FASTA files.

In the case of the *E. coli* test example, the first sequence is the chromosome, while the others are plasmids. To analyze only the chromosome, specify `-nchr 1`. To analyze the chromosome and the first plasmid, use `-nchr 2`, and so on (but please ensure that the homologous plasmids are located in the same order).

### Reference-based mode

This mode produces the alignmnet of all genomes
```
conda activate pannagram
${PATH_PAN}pannagram.sh -path_in ${PATH_DATA} -path_out ${PATH_OUT} -ref ${REF_NAME} -cores 8 -nchr 1
```

The result files are located at and called `ref_*.h5`. The prefix `ref` is crucial for the `analysis.sh` script.

```
PATH_ALN="${PATH_OUT}intermediate/consensus/"
cd ${PATH_ALN}
ls -lrt ref*
```

### REFERENCE-FREE mode

Number of chromosomes

```
conda activate pannagram
${PATH_PAN}pannagram.sh -path_in ${PATH_DATA} -path_out ${PATH_OUT} -cores 8 -nchr 1
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


## Paramenetrs for the reference genome

Attention:
- Names of genomes should not have symbol `.`


Reference genome can be from another species and being in another folder...

Number of chromosomes, path to the reference genome

Во всех референсах должно быть одно и число и во всех образцах должно быть одно число



