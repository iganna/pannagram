# Test Dataset

The test dataset comprises bacterial genomes (Escherichia coli) from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=562). You can download some genomes manually and put fasta-files into the `PATH_DATA` folder or utilize a tool for automatic downloading.

## Automated Option (Recommended):

1. Create the list of genome IDs as follows:

```bash
echo -e "GCA_009914755.4\nGCF_009914755.1\nGCA_964212195.1\nGCA_964212945.1\nGCA_964212425.1\nGCA_964212325.1" > "${PATH_DATA}human.txt"
echo -e "GCA_964212965.1\nGCA_964213075.1\nGCA_964212445.1\nGCA_964212265.1\nGCA_964213085.1\nGCA_964213115.1" >> "${PATH_DATA}human.txt"
echo -e "GCA_964212375.1\nGCA_964212495.1\nGCA_964212475.1\nGCA_964212205.1\nGCA_964213035.1\nGCA_964212285.1" >> "${PATH_DATA}human.txt"
echo -e "GCA_964212335.1\nGCA_964212305.1\nGCA_964213105.1\nGCA_964212215.1\nGCA_964212515.1\nGCA_964212245.1" >> "${PATH_DATA}human.txt"
echo -e "GCA_964212535.1\nGCA_964213055.1\nGCA_964212255.1\nGCA_964212935.1\nGCA_964212505.1\nGCA_964212465.1" >> "${PATH_DATA}human.txt"
echo -e "GCA_964212235.1\nGCA_964212525.1\nGCA_964212175.1\nGCA_964212365.1\nGCA_964212345.1\nGCA_964212975.1" >> "${PATH_DATA}human.txt"
echo -e "GCA_964212295.1\nGCA_964212485.1\nGCA_964212355.1\nGCA_964212405.1\nGCA_964212275.1\nGCA_964212315.1" >> "${PATH_DATA}human.txt"
echo -e "GCA_964213125.1\nGCA_964213015.1\nGCA_964212985.1\nGCA_964213045.1\nGCA_964212185.1\nGCA_964212545.1" >> "${PATH_DATA}human.txt"
echo -e "GCA_964213005.1\nGCA_964212225.1\nGCA_964212455.1\nGCA_964212415.1\nGCA_964212955.1\nGCA_964212995.1" >> "${PATH_DATA}human.txt"
echo -e "GCA_964213095.1\nGCA_964212385.1\nGCA_964213025.1\nGCA_964213065.1\nGCA_964212925.1\nGCA_000001405.29" >> "${PATH_DATA}human.txt"

```


```bash
echo -e "GCA_009858895.3\nGCF_009858895.2\nGCA_011545035.1\nGCA_011537615.1\nGCA_011537235.1\nGCA_011537295.1" >> "${PATH_DATA}covid.txt"
echo -e "GCA_011537695.1\nGCA_011537705.1\nGCA_011545165.1\nGCA_011545245.1\nGCA_011545275.1\nGCA_011741995.1" >> "${PATH_DATA}covid.txt"
echo -e "GCA_011742015.1\nGCA_011537005.1\nGCA_011537135.1\nGCA_011545235.1\nGCA_011536935.1\nGCA_009937895.1" >> "${PATH_DATA}covid.txt"
echo -e "GCA_009948465.1\nGCA_011537565.1\nGCA_011545535.1\nGCA_011537015.1\nGCA_011545555.1\nGCA_011537065.1" >> "${PATH_DATA}covid.txt"
echo -e "GCA_011537455.1\nGCA_009937935.1\nGCA_011536975.1\nGCA_009937905.1\nGCA_009937915.1\nGCA_009937925.1" >> "${PATH_DATA}covid.txt"
echo -e "GCA_009937945.1\nGCA_011537075.1\nGCA_011537085.1\nGCA_011537145.1\nGCA_011537265.1\nGCA_011537355.1" >> "${PATH_DATA}covid.txt"
echo -e "GCA_011537475.1\nGCA_011537505.1\nGCA_011537515.1\nGCA_011537525.1\nGCA_011537555.1\nGCA_011537825.2" >> "${PATH_DATA}covid.txt"
echo -e "GCA_011537835.2\nGCA_011537865.2\nGCA_011537895.2\nGCA_011537945.2\nGCA_011537975.2\nGCA_011537985.2" >> "${PATH_DATA}covid.txt"
echo -e "GCA_011538025.2\nGCA_011538055.2\nGCA_011538105.2\nGCA_011544975.2\nGCA_011545025.2\nGCA_011545065.2" >> "${PATH_DATA}covid.txt"
echo -e "GCA_011545075.2\nGCA_011545085.2\nGCA_011545095.2\nGCA_011545125.2\nGCA_011545285.2\nGCA_011545295.1" >> "${PATH_DATA}covid.txt"
echo -e "GCA_011545325.2\nGCA_011545335.2\nGCA_011545345.1\nGCA_011545355.1\nGCA_011545425.1\nGCA_011742025.1" >> "${PATH_DATA}covid.txt"
echo -e "GCA_011742035.1\nGCA_009938055.1\nGCA_009938065.1\nGCA_011537155.1\nGCA_011537715.1\nGCA_011537625.1" >> "${PATH_DATA}covid.txt"
echo -e "GCA_011742005.2\nGCA_011537365.1\nGCA_011537395.1\nGCA_011538015.2\nGCA_011545135.1\nGCA_011545505.1" >> "${PATH_DATA}covid.txt"
echo -e "GCA_011537465.1\nGCA_011537325.1\nGCA_009948555.1\nGCA_009948525.1\nGCA_011537785.1\nGCA_009948495.1" >> "${PATH_DATA}covid.txt"
echo -e "GCA_011537815.1\nGCA_011545455.1\nGCA_009937885.1\nGCA_011545495.1\nGCA_011537225.1\nGCA_011545545.1" >> "${PATH_DATA}covid.txt"
echo -e "GCA_009948425.1\nGCA_011537425.1\nGCA_011545485.1" >> "${PATH_DATA}covid.txt"


```


```bash
echo -e "GCA_011100685.1\nGCF_011100685.1\nGCA_012044875.1\nGCA_012045015.1\nGCA_044048985.1\nGCA_043643935.1" >> "${PATH_DATA}dog.txt"
echo -e "GCA_013276365.2\nGCA_031010295.1\nGCA_000002285.4\nGCF_000002285.5\nGCA_031010765.1\nGCA_031165255.1" >> "${PATH_DATA}dog.txt"
echo -e "GCA_008641055.3\nGCA_014441545.1\nGCF_014441545.1\nGCA_031010635.1\nGCA_004886185.2\nGCA_005444595.1" >> "${PATH_DATA}dog.txt"
echo -e "GCF_005444595.1" >> "${PATH_DATA}dog.txt"
```


```bash
echo -e "GCA_000208745.2\nGCF_000208745.1\nGCA_035896635.1\nGCA_036851125.1\nGCA_958328385.1\nGCA_036851095.1" >> "${PATH_DATA}caco.txt"
echo -e "GCA_958329045.1\nGCA_036851155.1\nGCA_958329735.1\nGCA_041222385.1\nGCA_041222375.1\nGCA_000403535.1" >> "${PATH_DATA}caco.txt"
echo -e "GCF_000403535.1" >> "${PATH_DATA}caco.txt"
```


```bash
echo -e "GCA_905475465.2\nGCA_905231885.1" >> "${PATH_DATA}peris_napi.txt"
echo -e "GCA_028984075.1\nGCA_029001895.1" >> "${PATH_DATA}peris_mannii.txt"
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






