# Building Conda package is currently under rework

## Current workflow
1. Create conda env
    ```sh
    conda env create -f pannagram.yaml
    ```
2. Activate conda env
    ```sh
    conda activate pannagram
    ```
3. 
    * Run `./user.sh` to build R package normally and create symlinks to `bash` scripts

        or
    * Run `./developer.sh` to build R package in a quick way (no documentation) and create symlinks to `bash` scripts


<!-- # Developer mode
1. Shell:
    ```sh
    conda env create -f pannagram_dev.yaml
    conda activate pannagram_dev
    ./developer.sh
    ```
    This will reinstall `pannagram` into your R environment with all its documentation.
2. Now you can call your `.sh` scripts:
    ```sh
    ./inst/pannagram.sh -path_in ... -path_out ... -cores 8 # ... etc.
    ```
    > **NOTE!** If you make any changes in `inst/**/*.R` files in order to see the results you must re-run `./developer.sh` to update internal `pannagram` package files accordingly!

# Builder mode:
1. Shell:
    ```sh
    conda create -n pannagram_build conda-build -c conda-forge
    conda activate pannagram_build
    conda-build -c bioconda -c conda-forge . --output-folder ./build --no-test
    zip -r build.zip build/
    ```
2. Now `build.zip` archive is ready to be shared.

# User mode
1. Shell:
    ```sh
    unzip build.zip
    $ conda env create -n pannagram -c bioconda -c conda-forge build/noarch/pannagram-*.tar.bz2
    conda activate pannagram
    ```
2. Now `.sh` scripts are linked to PATH and you are ready to use `pannagram` package! -->