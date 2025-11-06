# Build a Pangenome Alignment

The main module in the Pannagram package for building a pangenome linear alignment is also called `pannagram` :)
The alignment can be built in three modes:

- **Preliminary** – builds only a skeleton alignment to visually check the correspondence between aligned genomes.  
- **Reference-based** – aligns all genomes to a chosen reference genome.  
- **Reference-free** – creates a pangenome alignment independent of any reference genome.


The overall scheme of the alignment pipeline is presented below:

<div style="width: 70%;">
<p align="center">
  <img src="images/scheme_pipeline.png" style="width:70%; object-fit:cover;"/>
</p>
</div>