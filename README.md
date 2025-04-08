<!-- <p align="left">
<img src="images/pannagram_logo.png" width="30%" height="auto">
</p> -->

# Pannagram

<p align="left">
<img src="images/pannagram_scheme.png" width="100%" height="auto">
</p>

Pannagram is a package for constructing pan-genome alignments, analyzing structural variants, and translating annotations between genomes.
Additionally, Pannagram contains useful functions for visualization. The manual is available at the [pannagram-page](https://iganna.github.io/pannagram/).

## Setting Up the Working Environment

Follow these instructions to set up your Pannagram environment.

### Prerequisites

Make sure you have one of the following package managers installed:
- [Conda](https://docs.conda.io/projects/conda/en/latest/index.html)
- [Mamba](https://github.com/mamba-org/mamba)
- [Micromamba](https://github.com/mamba-org/mamba#micromamba)

Use your selected package manager by replacing ```<manager>``` with conda, mamba, or micromamba.

### Linux and macOS (Intel)

```bash
<manager> env create -f pannagram.yml
<manager> activate pannagram
```

### macOS (M-series chips)

```bash
<manager> env create --platform osx-64 -f pannagram_m4.yml
<manager> activate pannagram
```

### Alternative: Setting Up the Environment Without Explicit Versions

Use this option if you prefer an environment where package versions are not explicitly specified, and packages are installed with the latest compatible versions available:

**Linux and macOS (Intel)**

```bash
<manager> env create -f pannagram_min.yml
<manager> activate pannagram
```

**macOS (M-series chips)**

```bash
<manager> env create --platform osx-64 -f pannagram_min.yml
<manager> activate pannagram
```

### Running RStudio with the Environment

Make sure that [RStudio-Desktop](https://posit.co/download/rstudio-desktop/) is installed.
Then run the following in the command line:

```bash
<manager> activate pannagram
open -a RStudio
```

One may also create an alias:
```
alias panR="micromamba activate pannagram && open -a RStudio"
```

### Included Dependencies

The environment provides the following dependencies, each accessible directly via the command line:

- R interpreter (required version)
- [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/manual/manual.html)


### Windows users
Can try running code from this repo under [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (as Bash and `/` path separator are used extensively in the code). Nevertheless it was never tested in such environment.


## 1. Pangenome linear alignment

### 1.1 Building the alignment
Pangenome alignment can be built in two modes:
 - **reference-free**:
```sh
./pannagram.sh -path_in '<genome files directory path>' \
    -path_out '<output files path>' \
    -cores 8
```

 - **reference-based**:
```sh
./pannagram.sh -ref '<reference genome name>' \
    -path_in '<genome files directory path>' \
    -path_out '<output files path>' \
    -cores 8
```

 - **quick look**:
If there is no information on genomes and corresponding chromosomes available, one can run preparation steps:
```sh
./pannagram.sh -ref '<reference genome name>' \
    -path_in '<genome files directory path>' \
    -path_out '<output files path>' \
    -cores 8 -pre
```
 
An extended description of the parameters for all three scripts are avaliable by executing scripts with the flag `-help`.

### 1.2 Extract information from the pangenome alignment
Synteny blocks, SNPs, and sequence consensus (for the [IGV browser](https://igv.org)) can be extracted from the alignment:
```sh
./analys.sh -path_msa '<output path with consensus>' \
      -path_chr '<path with chromosomes>' \
      -blocks  \  # Find Synteny block inforamtion for visualisation
      -seq  \     # Create consensus sequence of the pangenome
      -snp        # SNP calling
```

### 1.3 Calling structural variants
When the pangenome linear alignment is built, SVs can be called using the following script:
```sh
./analys.sh -path_msa '<output path with consensus>' \
      -sv_call  \         # Create output .gff and .fasta files with SVs
      -sv_sim te.fasta \  # Compare with a set of sequences (e.g., TEs)
      -sv_graph           # Construct the graph of SVs
```

## 2. Visualisation
Pannagram contains a number of useful methods for visualization in R.

### 2.1 Visualisation of the pangenome alignment
All genomes together:
<p align="left">
<img src="images/pangenome_alignment.png" width="50%" height="auto">
</p>

A dotplot for a pair of genomes:
<p align="left">
<img src="images/syntenyplot.png" width="50%" height="auto">
</p>

### 2.2 Graph of Nestedness on Structural variants

Every node is an SV:
<p align="left">
<img src="images/graph_of_svs.png" width="50%" height="auto">
</p>

Every node is a unique sequence, size - the amount of this sequence in SVs:
<p align="left">
<img src="images/graph_of_svs_te.png" width="60%" height="auto">
</p>


### 2.3 Nucleotide plot for a fragment of the alignment

 - In the ACTG-mode:

<p align="left">
<img src="images/msaplot.png" width="50%" height="auto">
</p>

```r
# --- Quick start code ---
source('utils/utils.R')             # Functions to work with sequences
source('visualisation/msaplot.R')   # Visualisation
aln.seq = readFastaMy('aln.fasta')  # Vector of strings
aln.mx = aln2mx(aln.seq)            # Transfom into the matrix
msaplot(aln.mx)                     # ggplot object
```

- In the Polymorphism mode:

<p align="left">
<img src="images/msaplot_diff.png" width="50%" height="auto">
</p>


```r
# --- Quick start code ---
msadiff(aln.mx)                     # ggplot object
```
### 2.4 Dotplots of Sequences

Simultaneously on forward (dark color) and reverse complement (pink color) strands:
<p align="left">
<img src="images/dotplot.png" width="40%" height="auto">
</p>


```r
# --- Quick start code ---
source('utils/utils.R')             # Functions to work with sequences
source('visualisation/dotplot.R')   # Visualisation
s = sample(c("A","C","G","T"), 100, replace = T)
dotplot(s, s, 15, 9)                # ggplot object
```

### 2.5 ORF-finder and visualisation

<p align="left">
<img src="images/orfplot.png" width="40%" height="auto">
</p>

```r
# --- Quick start code ---
source('utils/utils.R')             # Functions to work with sequences
source('visualisation/orfplot.R')   # Visualisation
str = nt2seq(s)
orfs = orfFinder(str)
orfplot(orfs$pos)                   # ggplot object
```


## 3. Additional useful tools
### 3.1 Search for similar sequences

#### ...on the genome
The first approach involves searching against entire genomes or individual chromosomes. 
The quickstart toy-example is:
```sh
./simsearch.sh -in_seq genes.fasta -on_genome genome.fasta -out out.txt
```
The result is a GFF file with hits matching the similarity threshold.

#### ...on another set
The second approach, in contrast, is designed to search for similarities against another set of sequences. 
The quickstart toy-example is:
```sh
./simsearch.sh -in_seq genes.fasta -on_seq genome.fasta -out out.txt
```
The result is an RDS (R Data Structure) table. 
This table shows the coverage of one sequence over another and 
includes a flag column that indicates whether the sequences meet the similarity threshold. 
Additionally, the second script takes into account the coverage strand, 
determining not just if a sequence is covered, but also if it's covered in a specific orientation.

## Acknowledgements

**Development:**
- Anna Igolkina - Lead Developer and Project Initiator
- Alexander Bezlepsky - Assistant

**Testing:**
- Anna Igolkina: Lead Tester
- Anna Glushkevich: Testing the alignment on _A. lyrata_ genomes
- Elizaveta Grigoreva: Testing the alignment on _A. thaliana_ and _A. lyrata_ genomes
- Jilong Ma: Testing the SV-graph on spider genomes
- Alexander Bezlepsky: Testing the Pannagram's functionality on Rhizobial genomes
- Gregoire Bohl-Viallefond: Testing the annotation converter on _A. thaliana_ alignment

**Resources:**
- Logo was generated with the help of DALL-E
- Parallel Processing Tool: O. Tange (2018): GNU Parallel 2018, ISBN 9781387509881, DOI [https://doi.org/10.5281/zenodo.1146014](https://doi.org/10.5281/zenodo.1146014).
