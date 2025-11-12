# SV-calling

When the pangenome linear alignment is built, SVs can be called using the following command:
```sh
features -path_in '${PATH_PROJECT}' \
         -sv_call  \         # Create output .gff and .fasta files with SVs
         -sv_graph  \        # Construct the graph of SVs
         -sv_sim te.fasta \  # Compare with a set of sequences (e.g., TEs)
         -cores 8
```


* **Structural variants** calling. When the pangenome linear alignment is built, SVs can be called using the following command:
    ```sh
    features -path_project '${PATH_PROJECT}' \
        -sv_call  \         # Create output .gff and .fasta files with SVs
        -sv_sim te.fasta \  # Compare with a set of sequences (e.g., TEs)
        -sv_graph  \        # Construct the graph of SVs
        -cores 8
    ```


The data output can be found in `${PATH_PROJECT}/features/sv`
Newly generated figures will appear in `${PATH_PROJECT}/plots/sv`


## Main SVs data. the Ouput from -sv_call
Results can be found in `${PATH_PROJECT}/features/sv`
This folder contains the following main files:

** Files with coordinates of all SVs **
sv_pangen_pos.rds
sv_pangen_beg.rds
sv_pangen_end.rds

seq_sv_small.fasta
seq_sv_big.fasta

** Figures with statisticks on SVs are generated **



## Ouput data from -sv_graph

Figures cofferpond to all steps of graph genetating

Newly generated images will appear in `${PATH_PROJECT}/plots/sv`:

<img
    src="images/sv_graph.png"
    style="width: 90%; object-fit: cover;"
/>

## Generathe the graph of nestedness from scratch

If you have aset of sequences