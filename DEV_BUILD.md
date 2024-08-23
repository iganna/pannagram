# Developer mode
1. Shell:
    ```sh
    conda env create -f pannagram_dev.yaml
    conda activate pannagram_dev
    Rscript -e "devtools::install()"
    ```
2. Now you can call your `.sh` scripts as usual.

# Builder mode:
1. Shell:
    ```sh
    conda env create -n pannagram_build -c conda-forge conda-build
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
2. Now you are ready to use pannagram package!