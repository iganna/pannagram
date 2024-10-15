# Manupulations with the alignment

The entire workflow with the results is intended to be done in R. 
To access the important functions, please load the Pannagram library along with other required libraries from the `pannagram` conda environment.

```
library(pannagram)
library(rhdf5)
library(ggplot2)
```

# Visualisation of the alignment

Please define the working directory, which corresponded to ${PATH_OUT}, and some others to store results

```
path.work <- "my_output_forder"
path.figures <- paste0(path.work, 'figures/')
```

