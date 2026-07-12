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
git clone https://github.com/iganna/pannagram.git
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
git clone https://github.com/iganna/pannagram.git
cd pannagram
```


### 2. Prerequisites

Before proceeding, make sure you have one of the following package managers installed:

- [conda](https://docs.conda.io/projects/conda/en/latest/index.html)
- [mamba](https://github.com/mamba-org/mamba)
- [micromamba](https://github.com/mamba-org/mamba#micromamba)

In the commands below, replace `<manager>` with the package manager you’re using.


### 3. Environment Setup

#### **Linux and macOS (Intel)**
```bash
<manager> env create -f pannagram.yml
<manager> activate pannagram
```

#### **macOS (Apple Silicon / M-series)**
```bash
<manager> env create --platform osx-64 -f pannagram.yml
<manager> activate pannagram
```


### 4. Alternative: Fully Pinned Environment (Linux)

`pannagram.yml` lets the package manager resolve dependencies for your platform.
If you need an exact, reproducible build on **Linux (x86_64)**, use `pannagram_pinned.yml`,
which lists every package with a fixed version and build.

Install it with flexible channel priority — a fully-pinned spec otherwise trips the
solver's strict channel-priority check with spurious conflicts:

```bash
<manager> env create --channel-priority flexible -f pannagram_pinned.yml
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

### 6. Install Pannagram in Active RStudio Session

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
