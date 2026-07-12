# Pannagram

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)
![R](https://img.shields.io/badge/R-%E2%89%A5%204.0-blue)
![install: conda](https://img.shields.io/badge/install-conda-green)
![Platform](https://img.shields.io/badge/platform-linux--64%20%7C%20osx--64-lightgrey)

📖 [Documentation](https://iganna.github.io/pannagram/) · 📄 [Paper](https://doi.org/10.1101/2025.02.07.637071) · 🐛 [Issues](https://github.com/iganna/pannagram/issues)

<img src="docs/images/pannagram_scheme.png" width="90%" alt="Pannagram framework overview">

## Overview

Pannagram is a toolkit for building reference-free linear pangenome alignments and analyzing genomic polymorphisms.  
It consists of a command-line interface (CLI) for alignment construction, feature extraction, and sequence search, and an R library for downstream analysis and visualization.

Key capabilities:
- Reference-free pangenome alignment
- SNP and structural variant detection
- Mobile element family discovery
- Search for sequences in genomes
- Annotation liftover between genomes
- Visualization and sequence analysis

The CLI provides four commands:
- `pannagram` – reference-free (or reference-based) whole-genome alignment
- `features` – extract SNPs, structural variants, and mobile element families from an alignment
- `simsearch` – similarity search of sequences against sequence sets or whole genomes
- `chromotools` – reorder and rearrange chromosomes to match a common reference structure

Documentation can be found at [Pannagram-page](https://iganna.github.io/pannagram/).

## Requirements

- A Conda-compatible package manager: [conda](https://docs.conda.io/), [mamba](https://github.com/mamba-org/mamba), or [micromamba](https://github.com/mamba-org/mamba#micromamba)
- Linux (x86-64) or macOS (Intel; Apple Silicon via `--platform osx-64`)
- R ≥ 4.0 (installed automatically into the Conda environment)

Platform-specific details are in the [installation guide](https://iganna.github.io/pannagram/).

## Quick Installation

Clone the repository and create the conda environment:

```bash
git clone https://github.com/iganna/pannagram.git
cd pannagram
conda env create -f pannagram.yml
conda activate pannagram
./user.sh
./verify_installation.sh  # Verify the successful installation
```

## Quick Start

The typical workflow consists of two steps:
1. Build the pangenome alignment
2. Call genomic features from the alignment

Before running the example, set the following variables (preferably absolute paths) in the command line:

- `PATH_GENOMES` – directory containing input genome FASTA files
- `PATH_PROJECT` – directory where the project output will be stored

### Reference-Free Pangenome Alignment

Run the following command to perform a reference-free pangenome alignment:

```bash
pannagram  -path_genomes ${PATH_GENOMES} \
           -path_project ${PATH_PROJECT} \
           -cores 8
```

### Feature Calling

After the alignment step is complete, run the feature-calling module to identify all available genomic features:

```bash
features  -path_project ${PATH_PROJECT} \
          -synteny \
          -consensus \
          -snp \
          -snp_pi \
          -sv \
          -sv_families \
          -cores 8
```

All results will be saved under ${PATH_PROJECT} after both steps are complete:

```text
PATH_PROJECT/
├── features/     ← main analysis outputs
└── plots/        ← visualizations and figures
```

A detailed description of all output files and their formats is available in the documentation under **Getting Started → Output Data**.


## Pannagram R Library

In your R session, load the library:
```R
library(pannagram)
```

Pannagram R library provides functions for:
- working with FASTA files
- annotation liftover
- extracting specific pangenome regions as multiple sequence alignments
- ORF finding and visualization
- dot plots
- multiple sequence alignment visualization
- pangenome plots

For detailed documentation, visit the [Pannagram-page](https://iganna.github.io/pannagram/).

## Contributing

Questions, bug reports, and feature requests are welcome via [GitHub Issues](https://github.com/iganna/pannagram/issues).
For building Pannagram from source, see the [developer guide](DEV_BUILD.md).

## Citation

If you use Pannagram, please cite:

- **Pannagram: unbiased pangenome alignment and Mobilome calling**  
  *Anna A. Igolkina et al.*, *bioRxiv*, 2025. [**Link**](https://doi.org/10.1101/2025.02.07.637071)

To explore Pannagram applications, we recommend:

- **A comparison of 27 *Arabidopsis thaliana* genomes and the path toward an unbiased characterization of genetic polymorphism**  
  *Anna A. Igolkina et al.*, *Nature Genetics*, 2025. [**Link**](https://doi.org/10.1038/s41588-025-02293-0)

## License

Pannagram is released under the [MIT License](LICENSE.md).

## Acknowledgements

**Development:** Anna A. Igolkina (lead) and Alexander Bezlepsky.

**Testing (first release):** Anna Glushkevich, Elizaveta Grigoreva, Jilong Ma, and Grégoire Bohl-Viallefond.

**Tools:** GNU Parallel — O. Tange (2018), *GNU Parallel 2018*, [10.5281/zenodo.1146014](https://doi.org/10.5281/zenodo.1146014).
