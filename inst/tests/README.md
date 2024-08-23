# Testing Pannagram

Currently only some functions from `utils/utils.R` are covered with tests. The core testing package `testthat` was not added to `environment.yml`, so in order to run existing tests you should:

1. Activate the pannagram environment as explained in [root README](../README.md#recreating-working-environment)
2. Manually add `testthat` package to it:
    ```sh
    conda install r-testthat
    # OR
    mamba install r-testthat
    ```
3. Run tests (from the root of the repo):
    ```sh
    Rscript tests/testthat/test.R   
    ```