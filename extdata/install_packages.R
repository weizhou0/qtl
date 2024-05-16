#!/usr/bin/env Rscript
# install required R packages, from Finnge/SAIGE-IT

req_packages <- c(
    "BH",
    "data.table",
    "dplyr",
    "furrr",
    "Matrix",
    "methods",
    "optparse",
    "R.utils",
    "Rcpp",
    "RcppArmadillo",
    "RcppEigen",
    "RcppNumerical",
    "RcppParallel",
    "remotes",
    "RhpcBLASctl",
    "RSQLite",
    "SKAT",
    "SPAtest"
)
for (pack in req_packages) {
    if (!require(pack, character.only = TRUE)) {
        install.packages(pack, repos = "https://cloud.r-project.org", dependencies = TRUE)
        print(packageVersion(pack))
    }
}

github_packages <- c(
    "cysouw/qlcMatrix",
    "leeshawn/MetaSKAT"
)
for (pack in github_packages) {
    if (!require(pack, character.only = TRUE)) {
        remotes::install_github(pack)
    }
}
