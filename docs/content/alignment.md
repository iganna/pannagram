
# Construct the Multiple Genome Alignment

Before running the example, define the key working paths in the command line:
- `PATH_DATA` - the path to the data folder
- `PATH_PAN` - the path to the `pannagram` folder
- `PATH_OUT` - the path to the output directory


## Test Dataset

The test dataset comprises bacterial genomes (Escherichia coli) from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=562).  
GCA_000005845.2, GCA_000008865.2, GCA_042692195.1, GCA_042259155.1, GCA_042189615.1,  
GCA_042192135.1, GCA_042192365.1, GCA_042161425.1, GCA_042138005.1, GCA_964264615.1,  
GCA_042016495.1, GCA_042017145.1, GCA_042017895.1, GCA_040964775.2, GCA_041950585.1,  
GCA_041954705.1, GCA_041927255.1, GCA_041897545.1, GCA_002853715.1, GCA_000013265.1,  
GCA_003697165.2, GCA_000210475.1, GCA_000010385.1, GCA_000019645.1, GCA_003018455.1,  
GCA_024300685.1, GCA_027925825.1, GCA_040571355.1, GCA_003018035.1, GCA_027925745.1,  
GCA_027925805.1, GCA_027925765.1, GCA_027925785.1, GCA_027925565.1, GCA_027925845.1

You can download these genomes manually or utilize a tool for automatic downloading.

### Automated Option (Recommended):

1. Create a list of genome IDs as follows:

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
> Ensure that all accessions are approximately of the same size and verify that you haven't downloaded any accessions with suspicious lengths and number of chromosomes. For example, this command gives you three columns: the name of the file, the size, and the number of chromosomes.
> ```
> for file in "${PATH_DATA}"/*.fasta; do echo "$(basename "$file") $(du -h "$file" | cut -f1) $(grep -c "^>" "$file")"; done
> ```

## Preliminary mode

This mode might be useful for quickly reviewing the data.
As a result, you will get a visualization of draft reference-based alignments.
To set up the reference, please define the `REF_NAME` variable.
For example:
```
REF_NAME=GCA_000005845_2
```

And run Pannagram in preliminary mode:
```
conda activate pannagram
${PATH_PAN}pannagram.sh -path_in ${PATH_DATA} -path_out ${PATH_OUT} -ref ${REF_NAME} -cores 8 -pre 
```

The result is pairwise dot plots in PDF format:
```
cd "${PATH_OUT}plots/plots_${REF_NAME}/"
```

These plots are useful for deciding in which mode Pannagram should be used:
- `one2one`, when all chromosomes in all genomes are sorted in the same order and should be aligned one-to-one
- `all2all`, when all chromosomes should be compared with each other.

> 💡 **Tips:**  
> 1. If the chromosomes are not sorted but there is a one-to-one correspondence, it is recommended to reorder the original genome-files so that the chromosomes are arranged in the correct order for more convenient further analysis.
> 2. If some chromosomes are in the reverse complement compared orientation to the reference, it is advisable to correct them for clearer analysis later on, but this is not strictly required.
> 3. If, based on visual inspection, certain accessions appear suspicious and you decide not to analyze them, create a file called `ACC_FILE` that contains only the accessions you wish to analyze, with each accession listed on a separate line. For all further analyses, use the flag `-accessions ${ACC_FILE}` in every command. For examples:
> ```
> ACC_FILE="${PATH_DATA}ecoli_good.txt"
> echo -e "GCA_000005845.2\nGCA_000008865.2\nGCA_042692195.1\nGCA_042259155.1\nGCA_042189615.1" > "${ACC_FILE}"
> ```


## Reference-based mode

To run the alignments in the reference-based mode, the number of chromocomes for the analysis should be the same.
```
conda activate pannagram
${PATH_PAN}pannagram.sh -path_in ${PATH_DATA} -path_out ${PATH_OUT} -ref ${REF_NAME} -cores 8
```

The result files are located at and called `ref_*.h5`. The prefix `ref` is crucial for the `analysis.sh` script.

```
PATH_ALN="${PATH_OUT}intermediate/consensus/"
cd ${PATH_ALN}
ls -lrt ref*
```

## REFERENCE-FREE mode

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
> You can specify the number of references to use with the -nref flag or directly provide specific references using the -refs flag.


## Paramenetrs for the reference genome

Attention:
- Names of genomes should not have symbol `.`


Reference genome can be from another species and being in another folder...

Number of chromosomes, path to the reference genome

Во всех референсах должно быть одно и число и во всех образцах должно быть одно число



