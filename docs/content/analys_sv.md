# SV-calling

When the pangenome linear alignment is built, SVs can be called using the following command:
```sh
features -path_in '${PATH_PROJECT}' \
    -sv_call  \         # Create output .gff and .fasta files with SVs
    -sv_sim te.fasta \  # Compare with a set of sequences (e.g., TEs)
    -sv_graph  \        # Construct the graph of SVs
    -cores 8
```
Results can be found in `${PATH_PROJECT}/features/sv`

Newly generated images will appear in `${PATH_PROJECT}/plots/sv`:

<img
    src="images/sv_graph.png"
    style="width: 90%; object-fit: cover;"
/>
