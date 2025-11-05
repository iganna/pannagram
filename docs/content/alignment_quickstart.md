# Quick start

## Preliminary mode
```shell
pannagram -pre \
    -path_in ${PATH_DATA} \
    -path_out ${PATH_PROJECT} \
    -ref ${REF_NAME} \
    -cores 8
```
The refults are in the folder: `"${PATH_PROJECT}/plots/synteny_pairwise/${REF_NAME}/"`

## Reference-free mode

In the case of the *E. coli* test example, the first sequence is the chromosome, while the others are plasmids. To analyze only the chromosome, specify `-nchr 1`. To analyze the chromosome and the first plasmid, use `-nchr 2`, and so on (but please ensure that the homologous plasmids are located in the same order).
```shell
pannagram \
    -path_in ${PATH_DATA}\
    -path_out ${PATH_PROJECT} \
    -cores 8 \
```

The results are in the folder:
```
PATH_ALN="${PATH_PROJECT}/features/msa/"
ls -lrt ${PATH_ALN}/ref*
```

> 💡 **Tip:**  
> If you want to analyze only a subset of the genomes, specify the file with accessions for analysis:
> ```
> FILE_ACC="${PATH_DATA}ecoli_good.txt"
> echo -e "GCA_000005845.2\nGCA_000008865.2\nGCA_042692195.1\nGCA_042259155.1\nGCA_042189615.1" > "${FILE_ACC}"
> ```
> and then specisy the -accessions parameter:
> ```
> pannagram -path_in ${PATH_DATA} -path_out ${PATH_PROJECT} -accessions ${FILE_ACC} -ref ${REF_NAME} -cores 8 
> pannagram -path_in ${PATH_DATA} -path_out ${PATH_PROJECT} -accessions ${FILE_ACC} -cores 8 -nchr 1 
> ```




