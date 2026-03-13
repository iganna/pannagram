# Structural Variants (SVs)

To generate SV-related features, run the following command:

```bash
features -path_project '${PATH_PROJECT}' -sv -sv_family -sv_orf
```

Below is a detailed description of flags.

## SV-calling


After generating the consensus sequences for each genome, SVs can be called by running:
```bash
features -path_project '${PATH_PROJECT}' -sv
```

### SV name format

Example: `SVgr_1_id_0022|18`  
Meaning: `SV group` • `chromosome 1` • `ID 0022` | `length = 18 bp`

### FASTA files with SV sequences
are located at `${PATH_PROJECT}/features/sv/`

  - `seq_sv_large.fasta`: Large SVs (≥ 50 bp)
  - `seq_sv_short.fasta`: Short SVs (15–50 bp)  

### SV-Annotation per accession
are located at `${PATH_PROJECT}/features/sv/gff/`

Standard 9-column **GFF3** format.  
The **3rd column** (`type`) uses one of:

- `complex` — complex structural variant (SV)  
- `deletion` — simple SV, likely a deletion (absent in ≤10% of accessions)  
- `insertion` — simple SV, likely an insertion (present in ≤10% of accessions)  
- `indel` — simple SV with intermediate presence frequency among accessions

### General SV information

Stored in `sv_pangen_pos.rds`.

This table has the following columns:

- `gr` — SV group identifier
- `beg` — start position in pangenome coordinates
- `end` — end position in pangenome coordinates
- `len` — SV length in the pangenome
- `single` — indicator of SV type  
  - `1` — simple SV  
  - `0` — complex SV
- accession columns — length of the presence allele in each accession  
  - `0` — absence allele  
  - positive number — allele length  
  - `-1` — missing information
- `chr` — chromosome number  


### Positions of SVs in genomes
Two RDS tables provide per-genome coordinates for each SV:
- `sv_pangen_beg.rds`: Start positions  
- `sv_pangen_end.rds`: End positions

Each table has the following columns:
- `gr` — SV group ID (e.g., `SVgr_1_id_0100`)
- Genome accessions (`GCA_*`) — coordinate positions
- `chr` — chromosome number

Example:
```
            gr GCA_000005845.2 GCA_000008865.2 GCA_000010385.1 GCA_000013265.1 chr
SVgr_1_id_0101          226309          229642          228065          231053   1
SVgr_1_id_0102          226932          230264          228688          231675   1
SVgr_1_id_0103          227628          230959          229383          232370   1
SVgr_1_id_0104          227630          230961          229385          232372   1
SVgr_1_id_0105          228687          232018          230442          233431   1
```

### Output Figures

All figures are located at `${PATH_PROJECT}/plots/sv/`.
These figures are shown together below:

**a.** `sv_pie_chart.pdf`: Pie chart showing the distribution of simple and complex SVs.  
**b.** `sv_chr_minlen15_pangen.pdf`: Distribution of simple and complex SVs across chromosomes (SV length ≥ 15 nt).  
**c.** `sv_freq_hist.pdf`: Histogram of SV counts vs. frequency of presence (number of genomes with the presence allele).  
**d.** `sv_freq_hist_length_minlen15_abs.pdf`: Bar-plot of SV counts (length ≥ 15 nt) vs. frequency of presence, colored by SV length.  
**e.** `sv_freq_hist_length_minlen15_norm.pdf`: Normalized to 1 bar-plot of SV counts (length ≥ 15 nt) vs. frequency of presence, colored by SV length.


<div style="width: 50%;">
<p align="left">
  <img src="images/sv_stat.png" style="width:100%; object-fit:cover;"/>
</p>
</div>

## Families of Mobile Elements in SVs

After extracting SVs, Families of Mobile Elements can be identified by running:
```bash
features -path_project '${PATH_PROJECT}' -sv_family
```

### Families
are located at `${PATH_PROJECT}/features/sv/`.

- `edges_families.txt`: adjacency (edge) matrix describing the Graph of SV Families.  
Each connected component in the graph represents a family.  
Each edge indicates that the SV in the first column has a nested similarity to the SV in the second column.
- `sv_families.txt`: table assigning SVs to their corresponding families.


### Graphs of Families
are located at `${PATH_PROJECT}/plots/sv/`

**a-b.** `graph_*_families_colored.png`:  Graph of SV Families colored by length of SVs. Every node is an SV sequence - one place in the Pangenome.  
**c-d.** `graph_*_labeled.png`: Graph of SV Families showing the labels of Families. Labeled families make it easier to identify and locate families with interesting graph structures.

In the figure below, **a** and **c** correspond to *E. coli*, while **b** and **d** correspond to a butterfly.
<div style="width: 70%;">
<p align="left">
  <img src="images/sv_graph.png" style="width:100%; object-fit:cover;"/>
</p>
</div>



## ORFs in SVs

After extracting SVs, one can gen ORFs in long SVs with the following:
```bash
features -path_project '${PATH_PROJECT}' -sv_orf
```

The ourput file `seq_sv_large_orfs.fasta` is located at `${PATH_PROJECT}/features/sv/`.

### ORF name format

Example: `SVgr_1_id_0064|1336|ORF|940|20|-|aaLEN|307`  
Meaning: `<SV_name>|ORF|<nt_start>|<nt_end>|<strand>|aaLEN|<aa_length>`


## Compare SVs with a set of sequences of interest

To compare SVs against a specific set of sequences (for example, known transposons) run the following command:
```bash
features -path_project '${PATH_PROJECT}' -sv_sim <sequences_of_interest.fasta>
```

The results will be generated in `${PATH_PROJECT}/features/sv/`.  
The output file will be named `nestedness_sv_large_on_<sequences_of_interest>.txt`.  
It is a table with the following columns:

- **name.query** — Name of SV (query).
- **name.target** — Name of sequence from the set of interest (target).
- **strand** — Alignment strand (`+` or `-`).  
- **length.query** — Length (bp) of the SV sequence.  
- **length.target** — Length (bp) of the target sequence.  
- **coverage.query** — Fraction of the SV sequence covered by the target.  
- **coverage.target** — Fraction of the target sequence covered by the SV.  


The example and meaning are the following:

```
name.query          name.target strand  length.query  length.target  coverage.query  coverage.target
SVgr_1_id_0011|1349 M11300.1    +       1349          1341           0.9933283914    0.9977628635
SVgr_1_id_0026|712  AF439538.1  +       712           2047           1.0000000000    0.3483146067
SVgr_1_id_0026|712  DQ523609.1  +       712           3906           1.0000000000    0.1825396825
SVgr_1_id_0026|712  FJ609808.1  +       712           3906           1.0000000000    0.1825396825
SVgr_1_id_0026|712  L25844.1    +       712           601            0.8441011236    1.0000000000
SVgr_1_id_0064|1336 AY225450.1  +       1336          3833           0.9977544910    0.3477693712
SVgr_1_id_0064|1336 AY627636.1  +       1336          2785           0.9962574850    0.4779174147
SVgr_1_id_0064|1336 FJ664545.1  -       1336          3888           0.9992514970    0.3436213991
SVgr_1_id_0064|1336 KC806221.1  -       1336          3088           0.9962574850    0.4310233160
SVgr_1_id_0064|1336 KJ825883.1  +       1336          2548           0.9977544910    0.5231554160

```


