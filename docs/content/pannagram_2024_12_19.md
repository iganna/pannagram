# Workshop 2024.12.19

## 0. Install conda environment

### Download pannagram
```
git clone git@github.com:iganna/pannagram.git
git checkout dev
git pull
```

Create the environment:
```
cd pannagram
conda env create -f pannagram.yaml
./user.sh
```

## 1. Download genomes 
### Setup paths

Specity the working directory:
```
PATH_BASE=""
```

Specify other directories:
```
PATH_TOOLS="${PATH_BASE}"
PATH_DATA="${PATH_BASE}genomes/"
PATH_OUT="${PATH_BASE}alignment/"
mkdir -p ${PATH_DATA}

FILE_IDS="${PATH_DATA}ids.txt"
```

### Create the files with IDs
```
nano ${FILE_IDS}
```

### Download poputils
```
cd ${PATH_TOOLS}
git clone https://github.com/iganna/poputils.git
```

### Download genomes
```
cd ${PATH_TOOLS}poputils/genomes
./genbank_download_list.sh -f ${FILE_IDS} -p ${PATH_DATA}

ls -l ${PATH_DATA}*fasta
```

## 2. Perform the alignment in all modes
### Preliminary mode

Specify the reference genome:
```
REF_NAME=""
```

Run the Pannagram:
```
pannagram -path_in ${PATH_DATA} -path_out ${PATH_OUT} -ref ${REF_NAME} -pre -cores 8
```

Check the vvisualisation:
```
cd ${PATH_OUT}plots/plot_${REF_NAME}
```


### Preliminary mode

Remove the preliminary run:
```
cd ${PATH_BASE}
mv ${PATH_OUT} $(basename "$PATH_OUT")_pre
# rm -rf ${PATH_OUT}
```

Run the Pannagram:
```
pannagram -path_in ${PATH_DATA} -path_out ${PATH_OUT} -cores 8 -nchr 1
```

Check the vvisualisation:
```
cd ${PATH_OUT}plots/plot_${REF_NAME}
```

## 3. Get SNPs and SVs

Setup folders:
```
PATH_CHR="${PATH_OUT}intermediate/chromosomes/"
PATH_MSA="${PATH_OUT}intermediate/consensus/"
ALN_TYPE="msa_"
```

Run the analysis:
```
analys -path_msa ${PATH_MSA} -path_chr ${PATH_CHR} -seq -snp -sv -sv_graph
```

Results are here:
```
cd ${PATH_MSA}
```

<!-- Correspondence plot is here:
```
TODO
``` -->

## 4. Analyse-Visualise
Open R/RStudio.

### What is in the alignment file?

Setup:
```
library(hdf5)
library(ggplot2)
library(pannagram)

path.msa <- ''     # SPECIFY!
file.msa <- paste0(path, 'msa_1_1.h5')
```

Observe the alignment file:
```
h5ls(file.msa)
```

Get the coordinates of one accession:
```
acc <- ''   #SPECIFY!
v = h5read(file.msa, paste0("accs/", acc))
```

### Cut out a part of the alignment

Specify positions:

```
p.beg <- 10000
p.end <- p.beg + 5000
```

Get alignment and visualise:
```
aln.seq = cutAln(path.msa, 1, p, p+ 5000)
msaplot(aln.seq)
msadiff(aln.seq)
```


### gff2gff

```
gff1 = 
gff2 = gff2gff(path.msa, acc1 = '0', acc2 = 'pangen', gff1 = gff.tair, n.chr = 5, exact.match = F, aln.type = aln.type, s.chr = 'Chr')
```

## 5. Simsearch


### Against the set of sequences

```
simsearch -in_seq ${PATH_MSA}sv/seq_sv_big.fasta -on_seq <file_with_sequences> -out ${PATH_MSA}sv/on_seq/
```

### Against the folder with genomes

```
simsearch -in_seq "${PATH_MSA}sv/seq_sv_big.fasta" -on_path ${PATH_DATA}  -out "${PATH_MSA}sv/on_path/"
cd "${PATH_MSA}sv/on_path/"
ls -lrt
```

## 6. Find new Mobile elements

👻 Separate Workshop, I guess

<!-- 
### Components of similar SVs

Specify:
```
path.sv <- ''
```

```
res.cover.file = 'seq_sv_big_on_sv_cover.rds'
res.nest = readRDS(paste(path.sv, res.cover.file, sep = ''))
g.content = getGraphFromBlast(res.nest = res.nest, sim.cutoff = sim.cutoff, collapse = F)
g.comp <- getGraphComponents(g.content$edges)
```
 -->


































