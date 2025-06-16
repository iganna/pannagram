# Simsearch

* **...in set of sequences** This approach, in contrast, is designed to search for similarities against another set of sequences. 
    ```sh
    simsearch \
        -in_seq genes.fasta \
        -on_seq genome.fasta \
        -sim 90 \
        -out "<out path>"
    ```

* **...in the genome** This approach involves searching against entire genomes or individual chromosomes:
    ```sh
    simsearch \
        -in_seq genes.fasta
        -on_genome genome.fasta \
        -out "<out path>"
    ```
    The result is a GFF file with hits matching the similarity threshold.

* **...on all genomes in directory** Here we learch in all genomes in given directory:
    ```sh
    simsearch \
        -in_seq genes.fasta \
        -on_path "<path to genomes>" \
        -out "<out path>"
    ```