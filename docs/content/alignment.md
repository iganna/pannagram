
# Construct the Multiple Genome Alignment

## Test Dataset

The test dataset consists of bacterial genomes (*Rhizobium leguminosarum bv. viciae*) from NCBI. 
First, create a directory for the genomes:
```
mkdir data
echo -e "GCA_013391665.1\nGCA_013391685.1\nGCA_013391705.1\nGCA_012275595.1\nGCA_012275605.1\nGCA_012275615.1\nGCA_012275645.1\nGCA_012275665.1\nGCA_012275835.1\nGCA_012275875.1" > data/rhizobia.txt
```

For a quick download, one may use the script `genbank_download_list` from the repo [poputils](https://github.com/iganna/poputils).
To clone the repository, use the following command:

```bash
git clone https://github.com/iganna/poputils.git
cd poputils/genomes
```

Then, download the genomes using NCBI genome IDs. 
```
./genbank_download_list.sh -f ../../data/rhizobia.txt -p data
```

## Run Pannagram in the preliminary mode

This mode might be needed for a quick look at the data. As a result one whould have a visualisation of a draft reference-based alignments.

