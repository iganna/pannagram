# Manupulations with the alignment

After running `pannagram` pipeline you are able to get more features of your data! Pay attention, some of flags are independent from each other, others need to be passed together:
* **Extract information** from the pangenome alignment:
    ```sh
    features -path_in '${PATH_PROJECT}' \
        -blocks  \  # Find Synteny block inforamtion for visualisation
        -seq  \     # Create consensus sequence of the pangenome
        -snp \      # SNP calling
        -cores 8
    ```

* **Structural variants** calling. When the pangenome linear alignment is built, SVs can be called using the following command:
    ```sh
    features -path_in '${PATH_PROJECT}' \
        -sv_call  \         # Create output .gff and .fasta files with SVs
        -sv_sim te.fasta \  # Compare with a set of sequences (e.g., TEs)
        -sv_graph  \        # Construct the graph of SVs
        -cores 8
    ```

New files will be generated in respective subdirectories:
```shell 
PATH_PROJECT/
├── features/
│   ├── msa/               # .h5 files with MGA
│   ├── seq/               # Consensus sequence of the pangenome
│   ├── snp/               # VCF files with SNPs
│   └── sv/                # GFF and FASTA files with SVs
└── plots/
    ├── snp/               # Pi diversity plots
    ├── sv/                # Plots with SV statistics and SV graphs
    ├── synteny_pairwise/  # dotplot with REF-based synteny
    └── synteny_pangenome/ # synteny plots of all accessions
```