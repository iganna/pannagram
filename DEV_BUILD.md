# Building Conda package is currently under rework

## Current workflow
1. Create conda env
    ```sh
    conda env create -f pannagram.yml
    ```
2. Activate conda env
    ```sh
    conda activate pannagram
    ```
3. 
    * Run `./user.sh` to build R package normally and create symlinks to `bash` scripts

        or
    * Run `./developer.sh` to build R package in a quick way (no documentation) and create symlinks to `bash` scripts