# Alignment Summary

## Synteny Blocks

After running `pannagram` pipeline you are able to get more features of your data! For example you can get synteny blocks by running:
```sh
features -path_in '${PATH_PROJECT}' -blocks
```
Results you will find in `${PATH_PROJECT}/plots/synteny_pangenome`:

<div style="width: 70%;">
<p align="left">
  <img src="images/pangenome_alignment.png" style="width:70%; object-fit:cover;"/>
</p>
</div>



* **Extract Alignment Summary** from the pangenome alignment:
    ```sh
    features -path_project '${PATH_PROJECT}' \
        -blocks  \  # Find Synteny block inforamtion for visualisation
        -seq  \     # Create consensus sequence of the pangenome
        -snp \      # SNP calling
        -cores 8
    ```

    

## Consensus sequence


## Growth of Pangenome coordinate

