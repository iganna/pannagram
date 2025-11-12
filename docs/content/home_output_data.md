# Output Structure and Results

After running the pangenome alignment and feature calling steps, all results are saved in the output directory, which contains two main subfolders:

```
./  
├── features/     ← main output data
└── plots/        ← visualizations and figures  
```

> 💡 **Note:**  
> Output files may also include two numbers in the format `_N_M_`, where:  
> - `N` represents the chromosome number of the **query genome**  
> - `M` represents the chromosome number of the **reference genome**.


## Main Output Data in `features/`

This directory contains all processed data generated during the analysis.

```
features/
├── alignments/    ← Pangenome alignmetns
├── consensus/     ← Consensus pangenome sequences                        
├── snp/           ← SNP-calls
└── sv/            ← Structural Variants (SVs) and Families of Mobile Elements
  └── gff/         ← Annotation of SVs in every genome and pangenome coordinate
```

#### 1. Alignments — `alignments/`
  `ref_*.h5`: Reerence-based alignments per chromosomes.
  `pan_*.h5`: Pangenome alignments per chromosomes.

#### 2. Consensus Sequences — `consensus/`
  `seq_cons_*.fasta`: Consensus pangenome sequences where positions correspond to pangenome coordinates.

#### 3. SNPs — `snp/`
  `snps_*_pangen.vcf`: SNP variant calls per chromosome.  
  `snps_*_pangen_output.sites.pi`: Per-site nucleotide diversity (π) calculated using [VCFtools](https://vcftools.sourceforge.net).  
  `snps_*_pangen_dist.txt`: π-diversity between accessions calculated using [Plink](https://www.cog-genomics.org/plink/).

#### 4. Structural Variants — `sv/`
Contains detected Structural Variants, their sequences, and Families of Mobile Elements.

- **Positions:**  
  `sv_pangen_beg.rds`: Beginnings  
  `sv_pangen_end.rds`: Ends

- **Sequences:**  
  `seq_sv_large.fasta`: Large SVs (≥50 bp)  
  `seq_sv_short.fasta`: Short SVs (15–50 bp)

- **Families:**  
  `edges_families.txt`: adjacency (edge) matrix describing the Graph of SV Families. Each connected component in the graph represents a family.  
  `sv_families.txt`: table assigning SVs to their corresponding families.

- **Annotations:**  
  `gff/`: path with GFF files annotating SVs in individual genomes and in the pangenome coordinates.


## Output Figures in `plots/`

This directory contains all plots data generated during the analysis.

```
plots/
├── snp/                   ← SNP-based π-diversity plots per chromosome
├── sv/                    ← Visualizations of SV-related statistics
├── synteny_ref/           ← Dot-plots of pairwise synteny between reference accessions and the rest
│                            (generated after reference-based alignment steps)
└── synteny_pan/           ← Pangenome synteny plots between accessions per chromosome
```

#### 1. SNP-based π-diversity — `snp/`
  `snps_*_pangen_output.sites.pi_smooth.pdf`: π-diversity along chromosomes using a 200 kb window.

#### 2. Statistics on SVs and SV Families — `sv/`

- **Statistics:**  
  `sv_pie_chart.pdf`: Pie chart showing the distribution of simple and complex SVs.  
  `sv_chr_minlen15_pangen.pdf`: Distribution of simple and complex SVs across chromosomes (SV length ≥ 15 nt).  
  `sv_freq_hist.pdf`: Histogram of SV counts versus presence frequency (number of genomes with the presence allele).  
  `sv_freq_hist_length_minlen15_abs.pdf`: Bar-plot of SV counts (length ≥ 15 nt) versus presence frequency, colored by SV length.  
  `sv_freq_hist_length_minlen15_norm.pdf`: Normalized bar-plot (to 1) of SV counts (length ≥ 15 nt) versus presence frequency, colored by SV length.

- **Families:**  
  `graph_*_families_colored.png`:  Graph of SV Families colored by length of SVs. Every node - one place in the Pangenome.
  `graph_*_labeled.png`: Graph of SV Families showing the labels of Families.


#### 3. Reference-based synteny — `synteny_ref/`

  `NameRef-NameQuery.png`: Dot plots showing the synteny between the reference genomes (Y-axis) and the query genome (X-axis).

#### 4. Pangenome representation — `synteny_pan/`
  `pan_synteny_*.pdf`: Referene-free Pangenome Synteny between all genomes per chromosome.

