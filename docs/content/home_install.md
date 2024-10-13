# Installation

Recreating working environment.

## Linux users
Make sure you have [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) or [Mamba](https://github.com/mamba-org/mamba) installed. To create and activate the package environment run:
```sh
conda env create -f pannagram.yaml
conda activate pannagram
# OR
mamba env create -f pannagram.yaml
mamba activate pannagram
```
The environment downloads required R interpreter version and all needed libraries, including [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/), [MAFFT](https://mafft.cbrc.jp/alignment/software/manual/manual.html) and others.

## MacOS users
should also run:
```sh
brew install coreutils
```
to make sure all the needed shell commands are installed.

## Windows users
Can try running code from this repo under [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (as Bash and `/` path separator are used extensively in the code). However, it was never tested in such environment..
