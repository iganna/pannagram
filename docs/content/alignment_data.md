# Test Dataset

The test dataset comprises bacterial genomes (Escherichia coli) from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=562). You can download some genomes manually and put fasta-files into the `PATH_DATA` folder or utilize a tool for automatic downloading.

## Automated Option (Recommended):

1. Create the list of genome IDs as follows:

```bash
echo -e "GCA_000005845.2\nGCA_000008865.2\nGCA_042189615.1" > "${PATH_DATA}ecoli.txt"
echo -e "GCA_042016495.1\nGCA_042017145.1\nGCA_042017895.1" >> "${PATH_DATA}ecoli.txt"
echo -e "GCA_041954705.1\nGCA_000013265.1\nGCA_000210475.1" >> "${PATH_DATA}ecoli.txt"
echo -e "GCA_000010385.1\nGCA_000019645.1\nGCA_027925825.1" >> "${PATH_DATA}ecoli.txt"
echo -e "GCA_027925745.1\nGCA_027925805.1\nGCA_027925765.1" >> "${PATH_DATA}ecoli.txt"
echo -e "GCA_027925785.1\nGCA_027925565.1\nGCA_027925845.1" >> "${PATH_DATA}ecoli.txt"
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

To proceed with further analysis using this dataset, the reference genome variable should be set as follows:
```
REF_NAME=GCA_000005845.2
```
Additionally, specify the file with accessions for analysis:
```
FILE_ACC="${PATH_DATA}ecoli_good.txt"
echo -e "GCA_000005845.2\nGCA_000008865.2\nGCA_042692195.1\nGCA_042259155.1\nGCA_042189615.1" > "${FILE_ACC}"
```

# Run Pannagram

## Preliminary mode
```
conda activate pannagram
${PATH_PAN}pannagram.sh -path_in ${PATH_DATA} -path_out ${PATH_OUT} -ref ${REF_NAME} -cores 8
```
The refults are in the folder:
```
cd "${PATH_OUT}plots/plots_${REF_NAME}/"
```

## Reference-free mode

In the case of the *E. coli* test example, the first sequence is the chromosome, while the others are plasmids. To analyze only the chromosome, specify `-nchr 1`. To analyze the chromosome and the first plasmid, use `-nchr 2`, and so on (but please ensure that the homologous plasmids are located in the same order).
```
conda activate pannagram
${PATH_PAN}pannagram.sh -path_in ${PATH_DATA} -path_out ${PATH_OUT} -accessions ${FILE_ACC} -cores 8 -nchr 1
```
The results are in the folder:
```
PATH_ALN="${PATH_OUT}"
ls -lrt ${PATH_ALN}msa*
```






