# Installation

## Setting Up the Working Environment

After cloning the repo follow these instructions to set up your **pannagram** environment.

### Prerequisites

It is expected, that you have one of the following package managers installed on your machine:
- [Conda](https://docs.conda.io/projects/conda/en/latest/index.html)
- [Mamba](https://github.com/mamba-org/mamba)
- [Micromamba](https://github.com/mamba-org/mamba#micromamba)

Below package manager of your choice is refered to as `<manager>`.

### Linux

```bash
<manager> env create -f pannagram.yml
<manager> activate pannagram
```

### macOS (M-series chips)

```bash
<manager> env create --platform osx-64 -f pannagram_m4.yml
<manager> activate pannagram
```

### Alternative: Setting Up the Environment Without Explicit Versions

If you want to resolve package dependencies by yourself use `pannagram_min.yml` where only direct dependencies are specified with no explicit versions. Packages will be thus installed with the latest compatible versions available:

**Linux and macOS (Intel)**

```bash
<manager> env create -f pannagram_min.yml
<manager> activate pannagram
```

**macOS (M-series chips)**

```bash
<manager> env create --platform osx-64 -f pannagram_min.yml
<manager> activate pannagram
```

### Running RStudio with the Environment

Make sure that [RStudio-Desktop](https://posit.co/download/rstudio-desktop/) is installed.
Then run the following in the command line:

```bash
<manager> activate pannagram
open -a RStudio
```

One may also create an alias:
```bash
alias panR="<manager> activate pannagram && open -a RStudio"
```

### Installing Pannagram in Already Running RStudio Environment

We encourage you to run RStudio from activated package environment, but if RStudio server is already running and you want to install **pannagram** there as an R package run in console:
```R
setwd("<path to pannagram repo>")
source("install_in_rstudio.R")
```