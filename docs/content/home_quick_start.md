# Quick start

Before running the example, define the key working paths in the command line:
- `PATH_GENOMES` - the path to the data folder with genomes
- `PATH_PROJECT` - the path to the output directory

## Reference-Free Pangenome Alignment

Use this command to quickly generate a reference-free pangenome alignment:
```shell
pannagram  -path_genomes ${PATH_GENOMES} \
           -path_project ${PATH_PROJECT} \
           -cores 8
```

## Feature-calling

To call all available features, run:
```shell
features  -path_project ${PATH_PROJECT} \
          -blocks \
          -seq \
          -snp \
          -sv_calling \
          -sv_graph \
          -cores 8
```

## Results

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
  - Large SVs (>50 bp): `sv/seq_sv_large.fasta`  
  - Small SVs (15–50 bp): `sv/seq_sv_small.fasta`
- **Graph of Nestedness on SVs** `sv/edges_solved.txt` - edge matrix after filtarion steps

### Main Figures

All paths below are relative to `${PATH_PROJECT}/plots/`:  
- 







