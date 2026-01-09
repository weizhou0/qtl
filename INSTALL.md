# SAIGEQTL Installation Guide

## Quick Installation (Recommended)

The easiest way to install SAIGEQTL is using the `remotes` package:

```r
# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install SAIGEQTL from GitHub
remotes::install_github("weizhou0/qtl", dependencies = TRUE)
```

## System Requirements

- **R version**: >= 3.5.0
- **System dependencies**: 
  - GNU make
  - C++14 compatible compiler
  - OpenMP support (for parallel processing)
  - BLAS/LAPACK libraries

### Installing System Dependencies

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install build-essential libblas-dev liblapack-dev libomp-dev
```

#### CentOS/RHEL/Fedora
```bash
sudo yum install gcc-c++ openblas-devel lapack-devel libgomp
# or for newer versions:
sudo dnf install gcc-c++ openblas-devel lapack-devel libgomp
```

#### macOS
```bash
# Using Homebrew
brew install llvm libomp openblas lapack

# You may need to set compiler paths:
export PATH="/usr/local/opt/llvm/bin:$PATH"
export LDFLAGS="-L/usr/local/opt/llvm/lib"
export CPPFLAGS="-I/usr/local/opt/llvm/include"
```

## Alternative Installation Methods

### Method 1: Standard R Installation (Above)
Uses standard R package management - recommended for most users.

### Method 2: Pixi Environment Installation
Uses pixi for dependency management - useful for reproducible environments:

```bash
# Install GitHub dependencies first
CONDA_OVERRIDE_GLIBC=2.28 pixi run install-github-deps

# Install package to custom library location  
CONDA_OVERRIDE_GLIBC=2.28 pixi run --manifest-path=/path/to/qtl/pixi.toml \
    R CMD INSTALL /path/to/qtl \
    --library=/path/to/R/library \
    --no-lock --preclean
```

**Note**: Both methods provide identical functionality. Choose based on your preferred dependency management approach.

### Method 3: Manual Step-by-Step Installation

If you prefer to install dependencies manually:

```r
# 1. Install GitHub dependencies
remotes::install_github("weizhou0/SPAtest@v3.1.2")
remotes::install_github("cysouw/qlcMatrix")
remotes::install_github("leeshawn/MetaSKAT")
remotes::install_github("barkasn/fastSave")

# 2. Install SAIGEQTL
remotes::install_github("weizhou0/qtl")
```

## Verification

After installation, verify that SAIGEQTL works correctly:

```r
library(SAIGEQTL)

# Check that the package loads without errors
packageVersion("SAIGEQTL")

# Test basic functionality (optional)
# This will run a quick test if sample data is available
# SAIGEQTL::test_installation()
```

## Addressing Reviewer Concerns

### ✅ No cget Dependency
The package has been updated to **eliminate the problematic cget dependency**. The standard installation method (Method 1) uses only standard R build tools and system libraries. The pixi method is provided as an alternative for users who prefer conda environments, but is **not required**.

### ✅ No conda/Python Requirements  
Standard R installation works without any Python dependencies or conda environments. The package follows standard R package conventions and installs like any CRAN or Bioconductor package.

### ✅ Repository Naming
The GitHub repository name (`qtl`) differs from the package name (`SAIGEQTL`) for historical reasons. This is common practice in R packages. The package name to use in R is `SAIGEQTL`.

## Troubleshooting

### Common Issues

1. **"cget not found" errors**: Use Method 1 instead of conda/pixi methods
   - The standard installation no longer requires cget

2. **OpenMP not found**: Ensure your system has OpenMP support
   - Linux: Install `libgomp` or `libomp-dev`
   - macOS: Install `libomp` via Homebrew

3. **BLAS/LAPACK errors**: Install system BLAS/LAPACK libraries
   - See system-specific instructions above

4. **Compilation errors**: Ensure you have a C++14 compatible compiler
   - Update your compiler or install build tools

5. **Docker installation failures**: Use Method 1 within standard R Docker containers
   - No special environment setup needed

### Getting Help

- **GitHub Issues**: https://github.com/weizhou0/qtl/issues
- **Documentation**: https://weizhou0.github.io/SAIGE-QTL-doc/

## For Package Maintainers

### Bioconductor Submission

This package is being prepared for Bioconductor submission. Key compliance points:

- ✅ Standard R package structure
- ✅ No Python dependencies
- ✅ No conda environment requirements
- ✅ Proper DESCRIPTION file
- ✅ System requirements clearly documented
- ✅ Uses standard R installation mechanisms

### Development Installation

For developers working on the package:

```bash
# Clone the repository
git clone https://github.com/weizhou0/qtl.git
cd qtl

# Install in development mode
R CMD INSTALL . --preclean

# Or using devtools
Rscript -e "devtools::install('.', dependencies = TRUE)"

# Or with vignettes
Rscript -e "devtools::install('.', dependencies = TRUE, build_vignettes = TRUE)"
```

## Package Overview

SAIGEQTL is an R package for scalable and accurate expression quantitative trait locus (eQTL) mapping for single-cell studies. It implements a Generalized Poisson mixed model to handle:

- Repeat and complex data structure (multiple cells per individual, relatedness between individuals)
- Discrete read counts
- Large-scale data (20k genes, 200 cell types, millions of cells, millions of variants)
- Rare variant testing

### Two-Step Analysis Workflow

The package follows a two-step approach:

1. **Step 1**: Fit NULL GLMM model (`fitNULLGLMM_multiV`)
   - Estimate variance components
   - Calculate model parameters
   - Output: `.rda` model file and `.varianceRatio.txt`

2. **Step 2**: Association testing (`SAIGE_SPATest`)
   - Single variant tests
   - Gene/region-based tests
   - Uses SPA (Saddlepoint Approximation) for better p-value calibration

### Key Functions

#### Core Analysis Functions
- `fitNULLGLMM_multiV()` - Main null model fitting function
- `fitNULLGLMM_multiV_preprocess()` - Shared preprocessing for multiple phenotypes  
- `fitNULLGLMM_multiV_phenotype()` - Individual phenotype processing
- `fitNULLGLMM_multiV_batch()` - Batch processing for multiple phenotypes
- `SAIGE_SPATest()` - Association testing with SPA

#### Executable Scripts (in `extdata/`)
- `step1_fitNULLGLMM_qtl.R` - Command-line Step 1 script
- `step2_tests_qtl.R` - Command-line Step 2 script
- `cmd_qtl.sh` - Example analysis workflow

### Input Data Formats
- **PLINK**: bed/bim/fam files
- **VCF**: Compressed VCF with tabix index (.vcf.gz + .vcf.gz.tbi/.vcf.gz.csi)
- **BGEN**: BGEN v1.2 with .bgi index
- **Phenotypes**: Tab/space-delimited text files with headers

### Optimization Features
The package includes recent optimizations for multi-phenotype analysis:
- Shared preprocessing reduces computation time by 40-60% for multiple phenotypes
- Memory-efficient batch processing
- Optional parallel processing with OpenMP
- Pre-computed matrix operations reused across phenotypes

## Testing the Package

### Basic Functionality Test
```r
# Load the package
library(SAIGEQTL)

# Check that main functions exist
exists("fitNULLGLMM_multiV")  # Should return TRUE
exists("SAIGE_SPATest")       # Should return TRUE

# Check package version
packageVersion("SAIGEQTL")
```

### Regression Testing (For Developers)
```bash
# Run comprehensive regression tests (all formats)
./run_regression_test.sh comprehensive

# Run tests for specific input formats
./run_regression_test.sh test    # PLINK format
./run_regression_test.sh vcf     # VCF format
./run_regression_test.sh bgen    # BGEN format

# Compare existing outputs with baseline
./run_regression_test.sh compare

# Clean test outputs
./run_regression_test.sh clean
```

## Multi-Phenotype Analysis

Use optimized functions for multiple phenotypes:

```r
# Batch processing (recommended for multiple phenotypes)
results <- fitNULLGLMM_multiV_batch(
  phenoCols = gene_list,
  phenoFile = "data/gene_expression.txt",
  plinkFile = "data/genotypes", 
  nCores = 4  # Parallel processing
)

# Manual control
preprocess_data <- fitNULLGLMM_multiV_preprocess(...)
for (gene in gene_list) {
  result <- fitNULLGLMM_multiV_phenotype(preprocess_data, phenoCol = gene, ...)
}
```