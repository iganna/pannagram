# Input Data

The input data for **Pannagram** is a folder containing files with whole-genome assemblies at the chromosomal level.  
All genome files must be placed in a single folder.  
The following file extensions are supported: `(.fasta, .fna, .fa, .fas)`.  
This folder will be used as the input parameter for Pannagram.

### Genome Names  
Genome names are taken from the file names of the genome assemblies, without the file extensions.  
To simplify post-analysis steps, avoid using overly complex or lengthy file names.

### Using Reference Genomes  

If you plan to use Pannagram in **preliminary** or **reference-based** mode, you need to specify the name of the reference genome.  
Although we recommend storing the reference genome file in the same folder as the other genome assemblies, the reference genome can be located separately.

### Selecting Specific Genomes  

If you do **not** want to use all genomes from the folder, provide a text file listing only the desired ones.
Each line must contain a single accession, with no spaces or separators.

An example of such a file could be named `accessions.txt`:
```
GCA_000005845.2
GCA_000008865.2
GCA_042189615.1
```

## Quick Download and Create a Dataset from NCBI

If you don't already have a folder with genomes, you can quickly create a dataset by downloading genomes directly from NCBI.

### 0. Prepare Accessions List

- Open the [NCBI Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/) page.  
- Search for your species of interest and copy the GenBank accession IDs of the available assemblies. For alignment, you need at least two assemblies.
- Arrange these IDs in a file named `accessions.txt` — one accession per line.

> **Tips when selecting genomes:**
> - Check the year and sequencing technology — newer assemblies with long-read sequencing are usually better quality.  
> - Prefer Chromosomal-level assemblies over scaffold or contig-level.  
> - Ensure the genome sizes are roughly similar to avoid inconsistencies.

### 1. Define Paths

Open a Terminal and define the following variables with **absolute paths**:

- `PATH_DATA` — where you want to store genome files  
- `PATH_TOOLS` — where you usually cloned GitHub repositories  
- `FILE_ACCESSIONS` — the full path to your `accessions.txt` file  

Example command:
```bash
PATH_DATA=/absolute/path/to/data
PATH_TOOLS=/absolute/path/to/tools
FILE_ACCESSIONS=/absolute/path/to/accessions.txt
```

### 2. Clone the download helper [poputils](https://github.com/iganna/poputils)

```bash
cd ${PATH_TOOLS}
git clone https://github.com/iganna/poputils.git
```

### 3. Download Genomes from NCBI

Use the provided script to download genomes listed in your accessions file to the folder `PATH_DATA`:

```bash
cd ${PATH_TOOLS}/poputils/genomes
./genbank_download_list.sh -f ${FILE_ACCESSIONS} -p ${PATH_DATA}
```

### → Result
After running the commands, `PATH_DATA` will contain genome files for all accessions of interest, ready for downstream analysis.


## Test Dataset

If you just want to play with **Pannagram**, you can use a small test dataset with *E. coli* genomes.  
Specify the `FILE_ACCESSIONS` file as follows:

``` bash
echo -e "GCA_000005845.2\nGCA_000008865.2\nGCA_042189615.1" > "${FILE_ACCESSIONS}"
echo -e "GCA_042016495.1\nGCA_042017145.1\nGCA_042017895.1" >> "${FILE_ACCESSIONS}"
echo -e "GCA_041954705.1\nGCA_000013265.1\nGCA_000210475.1" >> "${FILE_ACCESSIONS}"
echo -e "GCA_000010385.1\nGCA_000019645.1\nGCA_027925825.1" >> "${FILE_ACCESSIONS}"
echo -e "GCA_027925745.1\nGCA_027925805.1\nGCA_027925765.1" >> "${FILE_ACCESSIONS}"
echo -e "GCA_027925785.1\nGCA_027925565.1\nGCA_027925845.1" >> "${FILE_ACCESSIONS}"
```

Then, run steps 1–3 of the pipeline described above to download the genomes from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=562).

## Working with Scaffold-Level Assemblies

If your genomes are available only at the scaffold level, you can create chimeric chromosomes by merging scaffolds in the corresponding order.  
We strongly recommend removing short or uninformative scaffolds beforehand.  
For example, when the goal is to work with Structural Variants (SVs) to analyse Mobile Elements, which are typically several kbp long, scaffolds shorter than 50–100 kbp will not provide sufficient resolution for accurate analysis.

### How to Merge Scaffolds Using Pannagram

1. Keep only the longest scaffolds from your assemblies.  
2. Run Pannagram in **preliminary mode**, using an accession with a chromosomal-level assembly as the reference.  
3. Execute the **Chromotools** module of Pannagram to generate chimeric chromosomes for all accessions.  

By following these steps, you can prepare chimeric genomes at the chromosomal level for comprehensive analysis with Pannagram.
