# Common annotation in pangenome coordinates

Gene annotations of the aligned genomes are arranged by the pangenome coordinate:
genes that occupy the same place form an **annotation group**.
The groups give a correspondence table between the annotations,
and every accession gets a full annotation, where the genes it lacks are projected from the other accessions.

```
features_n -path_project ${PATH_PROJECT} -annogroup ${PATH_ANNOTATION}
```

`${PATH_ANNOTATION}` is a folder with GFF files named after the accessions: `<accession>.gff`.
Not every accession needs an annotation: an accession without a file gets all its genes by projection.
With a single GFF file the step works as a liftover of this annotation to all the genomes.

## Options

| Option | Meaning |
|---|---|
| `-annogroup PATH` | Folder with `<accession>.gff` files |
| `-annot_chr FORMAT` | Sequence names in GFF are `<anything>FORMAT<number>`, the number is the chromosome. Default: `_Chr` (`acc_Chr1`). For names like `contig_1` use `-annot_chr _` |
| `-annot_feature TYPES` | Types of the top-level features to group, comma-separated. Default: all |
| `-annot_snap` | Keep a gene whose end is not aligned, moving the end inwards. Default: the gene is reported as lost |

Top-level features are those without `Parent`: `gene` in most annotations, `CDS`/`tRNA`/`rRNA` in Prokka output.
Their children (`mRNA`, `exon`, `CDS`, ...) follow them.

## How groups are formed

1. Every gene is transferred to the pangenome coordinates. A gene is lost if one of its ends
is not aligned or if the ends fall into different synteny blocks.
2. Two genes on the same strand are the same gene if they overlap by at least a half of the longer one.
A group is a set of genes connected by this relation.
Neighbouring genes that overlap by a few nucleotides stay in different groups.
3. A gene annotated as one model in some accessions and as two in the others
belongs to one group, but covers another one. It is shown in the covered group in brackets,
and nothing is projected to its accession there.

## Output

All files are in `${PATH_PROJECT}/features/annotation/`.

| File | Content |
|---|---|
| `annogroups.txt` | Groups: ID, position in the pangenome, number of accessions and genes |
| `correspondence.txt` | Group × accession: IDs of the genes. `-` – no gene, `(ID)` – a gene of another group covers this one |
| `lost_genes.txt` | Genes that were not placed, with the reason |
| `projection.txt` | Projected gene models with their status |
| `gff/<accession>.gff` | Annotation of the accession: own genes and projected ones, with common IDs |
| `pangen.gff` | All annotations in the pangenome coordinates, to view along the consensus sequence |

### Columns of `annogroups.txt`

- **group** – ID of the group: `AG<chromosome>G<number along the chromosome>`
- **chr, strand, beg, end, len** – position in the pangenome
- **n.acc, n.genes** – accessions that have a gene of the group, and the number of these genes
- **n.acc.multi** – accessions with more than one gene in the group
- **n.shared** – genes of other groups that cover this group
- **n.split** – number of groups in the cluster of overlapping genes the group belongs to
- **n.projected** – accessions that got the gene by projection
- **types** – types of the top-level features

### Reasons in `lost_genes.txt`

- **unaligned** – an end of the gene is not in the alignment
- **diff_blocks** – the ends are in different synteny blocks
- **chr_format** – the sequence name does not match the chromosome format
- **chr_not_aligned** – no alignment for this chromosome or accession

### Status of a projected model

A model is projected to an accession that has no gene of the group.
For coding models the sequence of the accession is checked against the model it was projected from:

- **intact** – the coding sequence keeps the frame, the start and the stop codons
- **frameshift** – the length of the coding sequence differs by a number not divisible by three
- **stop** – a premature stop codon
- **start_lost**, **stop_lost** – the start or the stop codon is gone
- **truncated** – a part of the coding sequence is absent or its boundary is not aligned
- **rearranged** – the ends of the gene are in different synteny blocks (not written to GFF)
- **absent** – the accession has no sequence at this place (not written to GFF)
- **unchecked** – the chromosome sequence was not available

Models without a coding sequence are **intact** if all their parts are in place and **truncated** otherwise.

### IDs in GFF files

A gene gets the ID `<group>.<accession>`, its parts – `<gene ID>.<type><number>`.
The original ID is kept in `orig_id`, the groups of the gene in `annogroup`.
Projected features have `pannagram` in the second column and the attributes
`projected_from`, `status` and `cov` (share of the model length present in the accession).
Genes that were not placed keep their original records.
