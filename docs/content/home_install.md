# Installation

This guide explains how to set up the **Pannagram** environment on various operating systems using popular package managers.
Once the setup is complete and the pannagram environment is activated, you can **run the following commands from anywhere** without having to remember the installation path:

```bash
pannagram
features
simsearch
chromotools

R && library(pannagram)  # Start R and load the library
```

## Quick Setup of the Working Environment

Follow these steps to quickly set up the Pannagram working environment using Conda.

```bash
git clone https://github.com/<user>/pannagram.git
cd pannagram
conda env create -f pannagram.yml
conda activate pannagram
./user.sh
./verify_installation.sh  # Verify the successful installation
```

## Detailed Setup of the Working Environment

This guide explains how to set up the **Pannagram** working environment on different platforms.


### 1. Clone the Repository

```bash
git clone https://github.com/<user>/pannagram.git
cd pannagram
```


### 2. Prerequisites

Before proceeding, make sure you have one of the following package managers installed:

- [conda](https://docs.conda.io/projects/conda/en/latest/index.html)
- [mamba](https://github.com/mamba-org/mamba)
- [micromamba](https://github.com/mamba-org/mamba#micromamba)

In the commands below, replace `<manager>` with the package manager you’re using.


### 3. Environment Setup

#### **Linux**
```bash
<manager> env create -f pannagram.yml
<manager> activate pannagram
```

#### **macOS (Apple Silicon / M-series)**
```bash
<manager> env create --platform osx-64 -f pannagram_m4.yml
<manager> activate pannagram
```


### 4. Alternative: Setup Without Explicit Versions

If you prefer to resolve package dependencies manually, use `pannagram_min.yml`, which includes only direct dependencies (no pinned versions).

#### **Linux and macOS (Intel)**
```bash
<manager> env create -f pannagram_min.yml
<manager> activate pannagram
```

#### **macOS (Apple Silicon / M-series)**
```bash
<manager> env create --platform osx-64 -f pannagram_min.yml
<manager> activate pannagram
```


### 5. Running RStudio Within the Environment

Make sure [RStudio Desktop](https://posit.co/download/rstudio-desktop/) is installed.

Then start RStudio from within the activated environment:
```bash
<manager> activate pannagram
open -a RStudio
```

Optionally, you can create a shortcut alias for convenience. For example:
```bash
alias panR="<manager> activate pannagram && open -a RStudio"
```

### 6. Installing Pannagram in an Already Running RStudio Session

While we recommend launching RStudio from the activated environment,  
you can also install **Pannagram** into an already running RStudio instance:

```R
setwd("<path to pannagram repo>")
source("install_in_rstudio.R")
```


### 7. Verifying Installation

After activating the environment, verify that Pannagram is available by running:

```bash
./verify_installation.sh
```
