# Structural variants (SVs) output

We do not provide VCF files for SVs.  
VCF requires the `SVTYPE` field (deletion, insertion, inversion, duplication) defined relative to a reference genome.  

Because our dataset is based on a pangenome representation, such reference-dependent classification may be misleading.  
In addition, complex structural variants cannot be meaningfully described using standard VCF SV types.

Instead, SVs are provided as:

1. FASTA files containing SV sequences
2. GFF files for each accession describing presence alleles  
3. Three tables describing SV coordinates and metadata  

The SV dataset includes:

| File                 | Description                      |
|----------------------|----------------------------------|
| `seq_sv_large.fasta` | sequences of SVs ≥50 bp          |
| `seq_sv_short.fasta` | sequences of SVs 15–50 bp        |
| `sv_pangen_pos.rds`  | general information about SVs    |
| `sv_pangen_beg.rds`  | start coordinates in each genome |
| `sv_pangen_end.rds`  | end coordinates in each genome   |
| `*.gff`              | SV annotations per accession     |

## SV classification

SVs are divided into two categories:

- **Simple SVs** — contain two alleles: presence and absence
- **Complex SVs** — contain multiple alleles with different sequence lengths

## SV name format

Example: `SVgr_1_id_0022|18`  
Meaning: `SV group` • `chromosome 1` • `ID 0022` | `length = 18 bp`


## FASTA files with SV sequences

  - `seq_sv_large.fasta`: Large SVs (≥ 50 bp)
  - `seq_sv_short.fasta`: Short SVs (15–50 bp)  

## SV annotation per accession

Standard 9-column **GFF3** format.  
The **3rd column** (`type`) uses one of:

- `multi` — complex structural variant (SV)  
- `deletion` — simple SV, likely a deletion (absent in ≤10% of accessions)  
- `insertion` — simple SV, likely an insertion (present in ≤10% of accessions)  
- `indel` — simple SV with intermediate presence frequency among accessions


## General SV information and positions in the pangenome

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


## Positions of SVs in genomes

Two RDS tables provide per-genome coordinates for each SV:
- `sv_pangen_beg.rds`: Start positions  
- `sv_pangen_end.rds`: End positions

Each table has the following columns:
- `gr` — SV group ID (e.g., `SVgr_1_id_0100`)
- Genome accessions name — coordinate positions in the accession
- `chr` — chromosome number


