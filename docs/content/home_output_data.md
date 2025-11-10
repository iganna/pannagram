# Output Structure and Results

After running the pangenome alignment and feature calling steps, all results are saved in the output directory, which contains two main subfolders:

```
./  
├── features/     ← main analysis outputs  
└── plots/        ← visualizations and figures  
```

## Main Output Data in `features/`

This directory contains all processed data generated during the analysis.

```
features/
├── msa/    ← Pangenome alignmetns
├── seq/    ← Consensus pangenome sequences                        
├── snp/    ← SNP-calls
└── sv/     ← Structural Variants (SVs) and Families of Mobile Elements
  └── gff/  ← Annotation of SVs in every genome and pangenome coordinate
```

#### 1. Alignments — `msa/`
Files `msa*.h5`: Pangenome alignments per chromosomes.

#### 2. Consensus Sequences — `seq/`
Files: `seq_cons_*.fasta`: Consensus pangenome sequences where positions correspond to pangenome coordinates.

#### 3. SNPs — `snp/`
Files `snps*pangen.vcf`: SNP variant calls per chromosome.

#### 4. Structural Variants — `sv/`
Contains detected Structural Variants, their sequences, and Families of Mobile Elements.

- **Positions:**  
  `sv_pangen_beg.rds` Beginnings  
  `sv_pangen_end.rds` Ends

- **Sequences:**  
  `seq_sv_large.fasta` Large SVs (≥50 bp)  
  `seq_sv_small.fasta` Short SVs (15–50 bp)

- **Families:**  
  `edges_solved.txt` — adjacency (edge) matrix describing the Graph of SV Families. Each connected component in the graph represents a family.  
  `sv_partition_solved.txt` — table assigning SVs to their corresponding families.

- **Annotations:**  
  `gff/` — path with GFF files annotating SVs in individual genomes and in the pangenome coordinates.


## Output Figures in `plots/`

This directory contains all plots data generated during the analysis.

```
plots/
├── snp/                   ← SNP-based π-diversity plots per chromosome
├── sv/                    ← Visualizations of SV-related statistics
├── synteny_pairwise/      ← Dot-plots of pairwise synteny between accessions 
│                            (generated after reference-based alignment steps)
└── synteny_pangenome/     ← Pangenome synteny plots between accessions per chromosome
```








