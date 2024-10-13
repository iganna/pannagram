# Test Dataset

The test dataset comprises bacterial genomes (Escherichia coli) from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=562). You can download some genomes manually and put fasta-files into the `PATH_DATA` folder or utilize a tool for automatic downloading.

## Automated Option (Recommended):

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

To proceed with further analysis using this dataset, the reference genome variable should be set as follows:
```
REF_NAME=GCA_000005845.2
```
Additionally, specify the file with accessions for analysis:
```
ACC_FILE="${PATH_DATA}ecoli_good.txt"
echo -e "GCA_000005845.2\nGCA_000008865.2\nGCA_042692195.1\nGCA_042259155.1\nGCA_042189615.1" > "${ACC_FILE}"
```
