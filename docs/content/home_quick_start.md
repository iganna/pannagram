# Quick Start

Before running the example, specify the key working directories (preferably using absolute paths) in the command line:

- `PATH_GENOMES` - Path to the folder containing input genomes
- `PATH_PROJECT` - Path to the output project directory


## Reference-Free Pangenome Alignment

Run the following command to perform a reference-free pangenome alignment:

```bash
pannagram  -path_genomes ${PATH_GENOMES} \
           -path_project ${PATH_PROJECT} \
           -cores 8
```

## Feature Calling

After the alignment step is complete, run the feature-calling module to identify all available genomic features:

```bash
features  -path_project ${PATH_PROJECT} \
          -blocks \
          -seq \
          -snp \
          -snp_pi \
          -sv \
          -sv_families \
          -cores 8
```

All results will be saved under ${PATH_PROJECT} after both steps are complete:

```
PATH_PROJECT/  
├── features/     ← main analysis outputs  
└── plots/        ← visualizations and figures  
```

A detailed description of all output files and their formats is provided in the “Output Data” section.

