# Workshop

## 0. Get pannagram

### Clone github repo
```
git clone https://github.com/iganna/pannagram.git
cd pannagram
```

### Replicate the environment
1. Platform-independent dependencies are given in `pannagram_min.yml`. To create the working environment, run one of the following commands depending on which tool you have installed:
    ```shell
    conda env create -f pannagram_min.yml
    conda activate pannagram
    ```

2. Then run script, which installs **pannagram** as an R package and generates docstrings:
    ```shell
    ./user.sh
    ```

## 1. Download genomes 
### Setup paths

Specity the absolute path to working directory:
```shell
PATH_BASE="<enter absolute path>/"
```

Specify other directories:
```shell
PATH_TOOLS="${PATH_BASE}"
PATH_DATA="${PATH_BASE}genomes/"
PATH_PROJECT="${PATH_BASE}alignment/"
mkdir -p ${PATH_DATA}

FILE_IDS="${PATH_DATA}ids.txt"
```

### Create the files with IDs
```shell
nano ${FILE_IDS}
```

### Download poputils
```shell
cd ${PATH_TOOLS}
git clone https://github.com/iganna/poputils.git
```

### Download genomes
```shell
cd ${PATH_TOOLS}poputils/genomes
./genbank_download_list.sh -f ${FILE_IDS} -p ${PATH_DATA}

ls -l ${PATH_DATA}*fasta
```

## 2. Perform the alignment in all modes
### Preliminary mode

Specify the reference genome (one id from `FILE_IDS`):
```shell
REF_NAME=""
```

Run the Pannagram:
```shell
pannagram -path_in ${PATH_DATA} -path_out ${PATH_PROJECT} -ref ${REF_NAME} -pre -cores 8
```

Check the visualisation:
```
cd ${PATH_PROJECT}plots/synteny_pairwise/${REF_NAME}
```


### Reference-free mode

Remove the preliminary run:
```shell
cd ${PATH_BASE}
mv ${PATH_PROJECT} $(basename "$PATH_PROJECT")_pre
# rm -r ${PATH_PROJECT}
```

Run the Pannagram:
```
pannagram -path_in ${PATH_DATA} -path_out ${PATH_PROJECT} -cores 8 -nchr 1
```

Check the visualisation:
```shell
cd ${PATH_PROJECT}plots/synteny_pairwise/${REF_NAME}
```

## 3. Get all features: synteny blocks, SNPs, SVs and SV-graph

Run the features script:
```
features -path_in ${PATH_PROJECT} -blocks -seq -snp -sv -sv_graph
```

Locations of some important result files:
```sh
workshop/alignment
├── features/
│   ├── msa/
│   │   └── msa_1_1.h5
│   ├── seq/
│   │   ├── seq_1_1.h5
│   │   └── seq_cons_1.fasta
│   ├── snp/
│   │   ├── snps_1_1_GCA_000005845.2.vcf
│   │   ├── snps_1_1_pangen.txt
│   │   └── snps_1_1_pangen.vcf
│   └── sv/
│       ├── gff/
│       ├── edges_solved.txt
│       ├── g_content_sim90.rds
│       ├── seq_sv_big.fasta
│       ├── seq_sv_big_on_sv_blast.txt
│       ├── seq_sv_big_on_sv_cover.rds
│       ├── seq_sv_small.fasta
│       ├── sv_in_graph_orfs.fasta
│       ├── sv_pangen_beg.rds
│       ├── sv_pangen_end.rds
│       ├── sv_pangen_pos.rds
│       └── sv_partition_solved.txt
└── plots/
    ├── sv/
    │   ├── graph_01_init.png
    │   ├── graph_02_refined_4_400.png
    │   ├── graph_03_louvain.png
    │   ├── graph_05_solved.png
    │   ├── graph_06_colored.png
    │   ├── graph_07_label.png
    │   ├── sv_chr_minlen15_pangen.pdf
    │   ├── sv_freq_hist_length_minlen15_abs.pdf
    │   ├── sv_freq_hist_length_minlen15_norm.pdf
    │   ├── sv_freq_hist.pdf
    │   └── sv_pie_chart.pdf
    ├── synteny_pairwise/
    └── synteny_pangenome/
        └── fig_synteny_chr1.pdf
```

<!-- Correspondence plot is here:
```
TODO
``` -->

## 4. Analyze-Visualize

Before entering R session manually copy to clipboard the current value from `PATH_PROJECT` variable from your terminal:

```shell
echo $PATH_PROJECT
```

### Open R/RStudio and prepare the workspace:

```R
library(pannagram)
library(rhdf5)
library(ggplot2)

path.project <- "<paste the path here>"
path.analisys <- file.path(path.project, 'analisys/')
if (!file.exists(path.analisys)) dir.create(path.analisys)
```

### Cut out a part of the alignment

Specify accession, chromosome number and an arbitrary window:
```R
acc <- "<put your accession id here>"
i.chr <- 1
window.beg <- 11000
window.end <- window.beg + 1000
```

Cut the specified window from the Multiple Genome Alignment:
```R
aln.seq <- cutAln(path.proj=path.project,
                  i.chr=i.chr,
                  p.beg=window.beg,
                  p.end=window.end,
                  acc=acc)
```

Visualize alignment in the window:

```R
p.nucl <- msaplot(aln.seq)
p.diff <- msadiff(aln.seq)

savePDF(p.nucl,
        path=path.analisys,
        name="msa_nucl",
        width=7,
        height=5)
savePDF(p.diff,
        path=path.analisys,
        name="msa_diff",
        width=7,
        height=5)
```

Save your window into FASTA format:
```R
sequences <- mx2aln(aln.seq)
writeFasta(sequences, file.path(path.analisys, "msa_window.fasta"))
```

### gff2gff

```R
file.gff1 <- "<your GFF file>"
acc1 <- "<set accession 1>"
acc2 <- "<set accession 2>"

gff1 <- read.table(file.gff1, stringsAsFactors=F)
gff1 <- gff1[gff1$V1 == gff1$V1[1],]
gff1$V1 <- paste0(gff1$V1, "Chr1")

gff2 <- gff2gff(path.proj, 
                acc1 = acc1,
                acc2 = acc2,
                gff1 = gff1,
                n.chr = 1,
                exact.match = F,
                s.chr = 'Chr')

write.table(gff2,
            file.path(path.analisys, "gff2.gff"),
            rownames=FALSE,
            colnames=FALSE,
            quotes=FALSE,
            sep='\t')
```

## 5. Simsearch

### Against the folder with genomes

```shell
simsearch \
    -in_seq "${PATH_PROJECT}/features/sv/seq_sv_big.fasta" \
    -on_path ${PATH_DATA} \
    -out "${PATH_BASE}/simsearch"
cd "${PATH_BASE}/simsearch"
ls -lrt
less -S simsearch.total_cnt_85_0.85.txt
```

### Against the set of sequences

```shell
simsearch \
    -in_seq "${PATH_PROJECT}/features/sv/seq_sv_big.fasta" \
    -on_seq '<file_with_sequences>' \
    -out "${PATH_BASE}/simsearch"
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


































