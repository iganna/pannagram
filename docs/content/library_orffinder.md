# ORF-Finder 

ORF-Finder scans the input DNA sequence in all six reading frames and extracts all stretches of amino acids that do not contain stop codons, regardless of whether they start with a canonical start codon (M).

## Input Sequences

The ORF-Finder accepts plain nucleotide strings:
```
seq <- "ACGTACGTACGT"
```

### Example nucleotide sequences
**[View](../examples/seqs_orf_finder.txt)** · **[Download](../examples/seqs_orf_finder.fasta)**  
This example contains a *Arabidopsis thaliana* LTR/Copia transposable element.  
To load these sequences in R, run:

```r
library(pannagram)

seq <- readFasta("seqs_orf_finder.fasta")
```

## Running the ORF-Finder

Use the function:
```R
orf.info = orfFinder(seq)
```
The result orf.info is a list with two main components:
- `orf.info$orf`: extracted ORF protein sequences
- `orf.info$pos`: a table describing the genomic positions of each ORF

The number of elements in `orf.info$orf` corresponds to the number of rows in `orf.info$pos`.  
The ORFs are sorted by length: `orf.info$orf[1]` is the longest ORF.

### Columns in `orf.info$pos`:
- **beg**: Start nt position of the ORF
- **end**: End nt position of the ORF 
- **shift**: Reading frame (0, 1, or 2) 
- **len**: ORF length in nucleotides 
- **aalen**: ORF length in amino acids 
- **strand**: Strand (“+” or “–”)

<!-- **Example output:**

```yaml
beg   end    shift   len   aalen   strand
39    2037   2       2631  877     +
38    276    2       1758  586     +
13    2692   0       444   148     +
6     1039   0       339   113     +
94    2032   2       327   109     -
83    4501   2       285   95      -
``` -->

## Visualization of ORFs

To visualize ORF positions along the sequence, run:
```R
orfplot(orf.info$pos)
```

<div style="width: 30%;">
<p align="left">
  <img src="images/orf_orfplot.png" style="width:100%; object-fit:cover;"/>
</p>
</div>

**Interpretation:**
- Red – forward strand
- Blue – reverse complement

There are two long ORFs, either of which can be extracted and, for example, analyzed using [protein BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp):
```
orfs.interest = c(orf.info$orf[1],
                  orf.info$orf[2])

writeFasta(orfs.interest, <your_file_to_store_orfs.fasta>)
```

### Recommended additional step
If the sequence you are analyzing is supposed to be a transposable element, it is helpful to inspect self-similarity via a dotplot:

```R
dotplot(seq, seq)
```
The sequence in our example is LTR/Copia retrotransposon, which often produces characteristic angled corner patterns in self-dotplots due to long terminal repeats (LTRs):
<div style="width: 20%;">
<p align="left">
  <img src="images/orf_dotplot.png" style="width:100%; object-fit:cover;"/>
</p>
</div>




