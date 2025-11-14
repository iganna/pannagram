# SNPs and Diversity

To generate SNP-related features, run the following command:

```bash
features -path_project '${PATH_PROJECT}' -snp -snp_pi
```

Below is a detailed description of flags.

## SNP-Calling

After generating the consensus sequences for each genome, SNPs can be called by running:
```bash
features -path_project '${PATH_PROJECT}' -snp
```

Output `VCF` files are at `${PATH_PROJECT}/features/snp/`.  
Each file `snps_*_pangen.vcf` corresponds to a chromosome.

## π-diversity

This option uses VCFtools and PLINK to estimate nucleotide diversity (π) per site and between genomes.
To run it, use the following command:

```bash
features -path_project '${PATH_PROJECT}' -snp_pi
```
#### Output Data

Files in `${PATH_PROJECT}/features/snp/`:  
`snps_*_pangen_output.sites.pi` – Per-site π-diversity calculated using [VCFtools](https://vcftools.sourceforge.net).   
`snps_*_pangen_dist.txt` – Pairwise π-diversity between accessions calculated using [Plink](https://www.cog-genomics.org/plink/).

#### Visual Output
Files in `${PATH_PROJECT}/plots/snp/`:  
`snps_*_pangen_output.sites.pi_smooth.pdf` – π-diversity along chromosomes (200 kb window)


Example:
<div style="width: 70%;">
<p align="left">
  <img src="images/snp_pi.png" style="width:100%; object-fit:cover;"/>
</p>
</div>

