# Alignment Summary

To generate an alignment summary, run the following command:

```bash
features -path_project '${PATH_PROJECT}' -synteny -consensus
```

Below is a detailed description of flags.

## Synteny Blocks

When the `apannagram` alignment is completed, you can generate synteny blocks and visualize the alignment by the following command:
```bash
features -path_project '${PATH_PROJECT}' -synteny
```

Results you will find in `${PATH_PROJECT}/plots/synteny_pan/`:

<div style="width: 70%;">
<p align="left">
  <img src="images/pangenome_alignment.png" style="width:70%; object-fit:cover;"/>
</p>
</div>

Each row represents a genome.  
The alignment is colored relative to the bottom genome: black indicates the same orientation as in the first genome, while pink marks inversions.  
Gray shading serves as a grid to help visually trace correspondence across genomes.

## Consensus sequence of the Pangenome

To generate the consensus sequence for the pangenome, run the following command:
```bash
features -path_project '${PATH_PROJECT}' -consensus
```

Results you will find in `${PATH_PROJECT}/features/consensus/`.
This folder contains `FASTA` files named `seq_cons_*.fasta`, each representing a consensus sequence for a corresponding chromosome.  
These sequences include the complete aligned content of the pangenome, excluding only highly unaligned regions (e.g., centromeres or assembly gaps).  
The consensus sequences can be used as input for various downstream analyses, including visualization in the `IGV` browser or use in genome annotation pipelines ([EDTA](https://github.com/oushujun/EDTA), [Helixer](https://github.com/usadellab/Helixer), etc.).


<!-- ## Growth of Pangenome coordinate -->

