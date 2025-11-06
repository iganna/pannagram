## Output Structure and Results

After running the alignment and feature detection steps, all results will be stored in the directory `${PATH_PROJECT}`.  
This folder contains two main subdirectories:

```
PATH_PROJECT/  
├── features/     ← main analysis outputs  
└── plots/        ← visualizations and figures  
```

###  Main Output Data `features/`

This directory contains all processed data generated during the analysis.

```

features/
├── msa/    ← Pangenome alignmetns
├── seq/    ← Consensus pangenome sequences                        
├── snp/    ← SNP-calls
└── sv/     ← Structural Variants (SVs) and Families of Mobile Elements
  └── gff/  ← Annotation of SVs in every genome and pangenome coordinate
```



#### **1. Pangenome Alignments — `msa/`**
- **Files:** `msa*.h5`  
  Contain pangenome alignments per chromosomes.

#### **2. Consensus Sequences — `seq/`**
- **Files:** `seq_cons_*.fasta`  
  Consensus pangenome sequences where positions correspond to pangenome coordinates.

#### **3. SNP Variants — `snp/`**
- **Files:** `snps*pangen.vcf`  
  VCF files containing SNP calls for each chromosome in pangenome coordinates.

#### **4. Structural Variants — `sv/`**
Contains detected Structural Variants, their sequences, and Families of Mobile Elements.

- **SV Positions:**  
  - Start coordinates: `sv_pangen_beg.rds`  
  - End coordinates: `sv_pangen_end.rds`  

- **SV Sequences:**  
  - Large SVs (≥50 bp):   `seq_sv_large.fasta`  
  - Small SVs (15–50 bp): `seq_sv_small.fasta`  

- **Families of Mobile Elements:**  
  - `sv/edges_solved.txt` — adjacency (edge) matrix describing relationships between SVs after filtering.

- **SV Annotations:**  
  - `sv/gff/` — GFF files annotating SVs in individual genomes and in the pangenome coordinates.




## Results

The results of the alignment and feature-detection will ba avaliabe at `PATH_PROJECT` and this folder will have the following structure:


.
├── features
│ ├── msa
│ ├── seq
│ ├── snp
│ └── sv
│   ├── gff
│   ├── seq_sv_large.fasta
│   ├── seq_sv_small.fasta
│   ├── sv_pangen_beg.rds
│   ├── sv_pangen_end.rds
│   ├── sv_pangen_pos.rds
│   ├── seq_sv_large_85_85.txt
│   ├── edges_solved.txt
│   └── sv_partition_solved.txt
└── plots
  ├── snp
  ├── sv
  ├── synteny_pairwise
  └── synteny_pangenome



The results will be saved in the following folders:  
- `${PATH_PROJECT}/features` — main output data  
- `${PATH_PROJECT}/plots` — visualizations of selected steps  

### Main Features

All paths below are relative to `${PATH_PROJECT}/features/`:  

- **Consensus sequences:** `seq/seq_cons_*.fasta` — consensus pangenome alignment (positions match pangenome coordinates).  
- **SNPs:** `snp/snps*pangen.vcf` — VCF-files in pangenome coordinates for each chromosome.  
- **SV positions:**  
  - Start: `sv/sv_pangen_beg.rds`  
  - End: `sv/sv_pangen_end.rds`  
- **SV sequences:**  
  - Large SVs (>=50 bp): `sv/seq_sv_large.fasta`  
  - Small SVs (15–50 bp): `sv/seq_sv_small.fasta`
- **Graph of Nestedness on SVs** `sv/edges_solved.txt` - edge matrix after filtarion steps

### Main Figures

All paths below are relative to `${PATH_PROJECT}/plots/`:  








