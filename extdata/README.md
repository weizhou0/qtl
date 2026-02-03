# SAIGEQTL Example Data and Scripts

This directory contains essential wrapper scripts and test data for SAIGEQTL.

## Essential Wrapper Scripts

These are the main command-line tools for SAIGEQTL analysis:

- **`step1_fitNULLGLMM_qtl.R`** - Step 1: Fit NULL GLMM model
- **`step2_tests_qtl.R`** - Step 2: Perform association tests  
- **`step3_gene_pvalue_qtl.R`** - Step 3: Calculate gene-level p-values
- **`makeGroupFile.R`** - Utility to create group files for gene-based tests

## Test Data

- **`input/`** - Example input data files for testing
- **`expected_output/`** - Expected output files for validation
- **`INSTALLATION_VALIDATION.md`** - Installation validation instructions

## Matrix Files

- **`*.mtx`** - Sparse matrix files (kinship, block-diagonal matrices)
- **`*_sampleID.txt`** - Sample ID files corresponding to matrices

## Usage

After installing SAIGEQTL, you can use these wrapper scripts directly:

```bash
# Step 1: Fit null model
Rscript extdata/step1_fitNULLGLMM_qtl.R --help

# Step 2: Run association tests  
Rscript extdata/step2_tests_qtl.R --help

# Step 3: Calculate gene p-values
Rscript extdata/step3_gene_pvalue_qtl.R --help

# Create group file
Rscript extdata/makeGroupFile.R --help
```

Or from Docker:
```bash
docker run --rm -v $(pwd):/data weizhou0/saigeqtl:latest step1_fitNULLGLMM_qtl.R --help
```