# Pangenome visualisation

In the Pannagram R package, there are three main functions for visualising genome comparisons:

1. **`panplot`** — provides a Pangenome-level correspondence between multiple accessions simultaneously in a stacked layout.
2. **`pangrowth`** — shows a pairwise correspondence between a single accession and the pangenome, visualised as a dotplot.  
This allows you to assess how the pangenome coordinate expands relative to the coordinate of the real genome.
3. **`syntenyplot`** — displays pairwise correspondence between a single accession and a reference genome, also shown as a dotplot.


## 1. `panplot` - Pangenome alignment 

Example:
```r
panplot(path.project, i.chr = 1)
```

Parameters:
- `path.project` — path to the Pannagram project directory.
- `i.chr` — chromosome/contig index to plot.
- `accessions` — (optional) subset of accessions and the order to display .
- `wnd.size` — (optional) binning window size along the pangenome, default: 1000000.

Example output:

<div style="width: 40%;">
<p align="left">
  <img src="images/panplot.png" style="width:100%; object-fit:cover;"/>
</p>
</div>


Description:
- Each horizontal row corresponds to a genome (accession).
- The x-axis shows the accession’s own genome coordinates.
- Connectors between accessions indicate sequence correspondence.
- Shades of grey act as a grid marking shared pangenome coordinates.
- White gaps indicate missing correspondence.
- Pink regions highlight inversions, referenced to the bottom accession.


On the plot one can see, for example, an inversion in accession 6069 and a stable translocation in the population on the right half of the chromosome.


## 2. `pangrowth` - Pangenome growth

Example:
```r
pangrowth(path.project, acc='6909', i.chr=1)
```

Parameters:
- `path.project` — path to the Pannagram project directory.
- `acc` — the accessions to compate with the Pangenome.
- `i.chr` — chromosome/contig index to plot.

Example output:
<div style="width: 30%;">
<p align="left">
  <img src="images/pangrowth.png" style="width:100%; object-fit:cover;"/>
</p>
</div>

In the middle part of the plot the pangenome coordinate correspondence shows an increased growth rate,  
indicating many structural variants in this centromeric region.


## 3. `syntenyplot`: Pairwise correspondence

Example:
```r
syntenyplot(path.project = path.project, acc = 'GCA_000008865.2', ref='GCA_000005845.2')
```

Parameters:
- `path.project` — path to the Pannagram project directory.
- `acc` — the accessions to compate with the Pangenome.
- `ref` — The name of the reference genome.

Example output:
<div style="width: 30%;">
<p align="left">
  <img src="images/syntenyplot.png" style="width:100%; object-fit:cover;"/>
</p>
</div>

This plot displays all chromosomes aligned between the accession and reference genomes.
