## usage
This conda recipe compiles and installs the `r-saigeqtl` conda package. It is designed to be run by the [bioconda.yaml GitHub action workflow](../.github/workflows/bioconda.yaml).

The workflow will generate conda packages for any commit to main or a PR that targets main. You can download and test any of the packages:
1. First, navigate to the Actions tab and click on the most recent successful execution:
    https://github.com/weizhou0/SAIGEQTL/actions/workflows/bioconda.yaml?query=branch%3Amain+is%3Asuccess
2. Now scroll down to the section titled "Artifacts" and download the zip archive file that is most appropriate for your platform. At the moment, the recipe supports the following.
    - Linux (x86) aka _Ubuntu_
    - Linux (aarch64) aka _Ubuntu ARM_
    - macOS 15 (x86)
    - the latest macOS (arm64)
    <img width="200" alt="Artifacts" src="https://github.com/user-attachments/assets/a2c46c35-08ff-46b0-bb88-c70d376833f7" />
3. To install the package, run the following via the command-line
    ```bash
    unzip -d r-saigeqtl r-saigeqtl-*.zip # unzip the artifact from GitHub
    conda install -c conda-forge -c bioconda -c "file://$PWD/r-saigeqtl" r-saigeqtl
    rm -r r-saigeqtl # delete the artifact, which is no longer needed
    ```
4. Done! Now, you can run `step1_fitNULLGLMM_qtl.R --help`, for example.

## maintenance
Occasional maintenance of the bioconda recipe must be performed. This may include:

1. Updating the list of dependencies in the `requirements` section of the [meta.yaml](r-saigeqtl/meta.yaml) file when new packages are added or old ones removed
2. Updating the patch files in the [patches/](r-saigeqtl/patches) directory whenever a substantial change is made to the `DESCRIPTION` or `src/Makevars` files
3. Updating the list of `.R` scripts in the `test` section of the [meta.yaml](r-saigeqtl/meta.yaml) file and in the [build.sh](r-saigeqtl/build.sh) file when new scripts are added or old ones removed. In this case, the shebang of the `.R` script might also need to be patched.

You can search online for help with any of these.
For example, you can refer to [this tutorial](https://bioconda.github.io/tutorials/2024-debugging-bioinformatic-software-to-bioconda.html) for debugging broken recipes on the Bioconda website.
You can refer to [this blog](https://www.anaconda.com/blog/patching-source-code-to-conda-build-recipes) for help creating and maintaining patch files.
ChatGPT can also be very helpful if you share the error message from the logs and the changes you've made to your `DESCRIPTION` file, and ask it how to update the recipe files.

Typically, maintenance of a bioconda recipe would also include updates to the metadata of the bioconda package at the top of the [meta.yaml](r-saigeqtl/meta.yaml) file (such as the latest version number, commit, and SHA256). However, the [bioconda.yaml GitHub action workflow](../.github/workflows/bioconda.yaml) should handle this automatically.

If you are a user of SAIGEQTL and you'd like to contribute to this recipe, please fork this repository and create a pull request. You can [enable GitHub actions on your fork](https://github.com/orgs/community/discussions/176785#discussioncomment-14667305) to test out your changes. The most common maintenance task is #2 (updating the patch files).

## publishing this recipe on Bioconda
Once SAIGEQTL is ready, this recipe should be [copied over to the Bioconda repository](https://bioconda.github.io/contributor/index.html) so that users can install it via a single command.
```
conda install -c conda-forge -c bioconda r-saigeqtl
```
If you are a user of SAIGEQTL, you can [click here](https://anaconda.org/search?q=saigeqtl) to check whether someone has uploaded the recipe yet.

When copying the recipe, you should make sure to manually update the metadata at the top of the file.

All Bioconda packages are automatically converted into Docker images and published to the [Biocontainers Registry](https://biocontainers-edu.readthedocs.io/en/latest/conda_integration.html#defining-a-conda-package), so users would be able to use the Docker images in that repository at that point.

## one caveat
All dependencies of a conda package must also be published to conda. The `fastSave` package is the only dependency that isn't already published to conda, so this recipe doesn't install `fastSave`, unfortunately. In the future, someone could create a recipe for `fastSave` and publish it to conda, as well.
