# Simsearch

Simsearch is a Pannagram module designed to compare query sequences with different sequence sources, applying similarity and coverage cutoffs.

It can search for query sequences in:
- **`I.`** another set of sequences
- **`II.`** a genome
- **`III.`** a collection of genomes
The module uses BLAST to identify initial matches and then merges neighboring hits to produce the final search results.

### Common Parameters for all modes:

- `-query_seq <query_sequences.fasta>`  
  **Required.** Input FASTA file containing the query sequences.
- `-out <output_path>`  
  **Required.**
- `-similarity <similarity_value>`  
  Optional. Similarity cutoff (default: `85`).
- `-coverage <coverage_value>`  
  Optional. Coverage cutoff (default: same as `<similarity_value>`).

## `I.` Search for query sequences in another (target) set of sequences

To run this mode, use the following command:
```bash
simsearch \
    -query_seqs <query_sequences.fasta> \
    -target_seqs <target_sequences.fasta> \
    -out "<output_path>"
```

The output file is named `target_sequences_X_Y.txt`, where `X` represents the similarity cutoff and `Y` represents the coverage cutoff.  
This file contains the table with the following colmns:

- **name.query** — Name of a query sequence.
- **name.target** — Name of a target sequence.
- **strand** — Alignment strand (`+` or `-`).  
- **length.query** — Length (bp) of the query sequence.  
- **length.target** — Length (bp) of the target sequence.  
- **coverage.query** — Fraction of the query sequence covered by the target.  
- **coverage.target** — Fraction of the target sequence covered by the query.  


## `II.` Search for query sequences in a genome

To run this mode, use the following command:
```sh
simsearch \
    -query_seqs <query_sequences.fasta> \
    -target_genome <target_genome.fasta> \
    -out "<output_path>"
```

**Output Files:**
- `target_genome_X_Y.gff` — GFF file containing hits that meet the specified similarity and coverage cutoffs.  
- `target_genome_X_Y.table` — Same results as in the GFF file, but in a tabular format with descriptive column names.  
- `target_genome_X_Y.cnt` — Summary file with counts of query sequences found per chromosome and in total across the target genome.

`X` represents the similarity cutoff and `Y` represents the coverage cutoff.

## `III.` Search for query sequences in a collection of genomes


To run this mode, use the following command:
```sh
simsearch \
    -query_seqs <query_sequences.fasta> \
    -target_path <path_with_genomes> \
    -out "<output_path>"
```

Output files:
- The same files as described above at **`II`** for each genome in the collection.
- A summary file `total_counts_X_Y.txt`, containing total counts of query sequences for every genome.

