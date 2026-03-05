# Build a Pangenome Alignment

The main module in the Pannagram package for building a pangenome linear alignment is also called `pannagram` 🤓   
Alignments are performed on a per-chromosome basis.  
Being linear, this alignment also accounts for inversions and translocations within a single chromosome.  

The alignment can be built in three modes:

- **`I` Preliminary** – builds only a skeleton alignment to visually check the correspondence between aligned genomes.  
- **`II` Reference-based** – aligns all genomes to a chosen reference genome.  
- **`III` Reference-free** – creates a pangenome alignment independent of any reference genome.

The three modes are not independent — each subsequent mode builds upon the results of the previous one: the Preliminary mode precedes the Reference-based mode, and several Reference-based alignments are utilised to construct the Reference-free mode.  
The overall scheme of the alignment pipeline is shown below:

<div style="width: 70%;">
<p align="center">
  <img src="images/scheme_pipeline.png" style="width:70%; object-fit:cover;"/>
</p>
</div>

The alignment pipeline is coded in a modular way so that some components can be replaced with faster modules in the future, or new blocks can be added to improve alignment.  
For example, in the near future, an Ancestral Recombination Graph (ARG) is planned to be introduced to achieve better alignment of complex regions and to phase complex Structural Variants.

🤝 **If you have suggestions** or methods to be incorporated into the Pannagram alignment pipeline, please feel free to contact Anna Igolkina via [email](mailto:igolkinaanna11@gmail.com).
