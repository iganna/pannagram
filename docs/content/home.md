
# Pannagram

Pannagram is a package for constructing pan-genome alignments, analyzing structural variants, and translating annotations between genomes.
Additionally, Pannagram contains useful functions for visualization. The manual is available [here](https://iganna.github.io/pannagram/).

## Installation

Recreating working environment.

### Linux users
Make sure you have [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) or [Mamba](https://github.com/mamba-org/mamba) installed. To create and activate the package environment run:
```sh
conda env create -f pannagram.yaml
conda activate pannagram
# OR
mamba env create -f pannagram.yaml
mamba activate pannagram
```
The environment downloads required R interpreter version and all needed libraries, including [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/), [MAFFT](https://mafft.cbrc.jp/alignment/software/manual/manual.html) and others.

### MacOS users
should also run:
```sh
brew install coreutils
```
to make sure all the needed shell commands are installed.

### Windows users
Can try running code from this repo under [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (as Bash and `/` path separator are used extensively in the code). However, it was never tested in such environment..

## Acknowledgements

**Development:**
- Anna Igolkina - Lead Developer and Project Initiator
- Alexander Bezlepsky - Assistant

**Testing:**
- Anna Igolkina: Lead Tester
- Anna Glushkevich: Testing the alignment on _A. lyrata_ genomes
- Elizaveta Grigoreva: Testing the alignment on _A. thaliana_ and _A. lyrata_ genomes
- Jilong Ma: Testing the SV-graph on spider genomes
- Alexander Bezlepsky: Testing the Pannagram's functionality on Rhizobial genomes
- Gregoire Bohl-Viallefond: Testing the annotation converter on _A. thaliana_ alignment

**Resources:**
- Logo was generated with the help of DALL-E
- Parallel Processing Tool: O. Tange (2018): GNU Parallel 2018, ISBN 9781387509881, DOI [https://doi.org/10.5281/zenodo.1146014](https://doi.org/10.5281/zenodo.1146014).