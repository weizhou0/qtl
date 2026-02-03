# SAIGEQTL Installation Validation

This directory contains expected output files and validation tools to verify that your SAIGEQTL installation is working correctly.

## Quick Validation

After installing SAIGEQTL, run this command from the package root directory:

```bash
./testing/run_regression_test.sh
```

This will:
1. Run Step 1 (NULL GLMM fitting) 
2. Run Step 2 (association testing) for PLINK format
3. Compare results against expected outputs
4. Report success/failure status

## Comprehensive Validation

To test all supported input formats (PLINK, VCF, BGEN):

```bash
./testing/run_regression_test.sh comprehensive
```

This runs the full validation suite and compares results across all formats.

## Validation Criteria

A successful validation should show:
- ✅ All steps complete without errors
- ✅ Output files are generated with expected content
- ✅ Cross-format correlations > 0.999 (for comprehensive test)
- ✅ Variance ratio estimates match reference values within tolerance

## Expected Output Files

The `expected_output/` directory contains reference results:

- `nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_plink_cis` - PLINK format association results  
- `nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_vcf_cis` - VCF format association results
- `nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_bgen_cis` - BGEN format association results
- `nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared.varianceRatio.txt` - Variance ratio estimates

## Troubleshooting

### Test Fails with Library Path Error

The validation script will automatically try to detect your SAIGEQTL installation. If detection fails, you'll see helpful instructions for finding your installation.

**Find your R library path:**
```bash
# Check where R packages are installed
R -e '.libPaths()'

# Look in common locations
ls -d ~/R/library/SAIGEQTL 2>/dev/null
ls -d /usr/local/lib/R/site-library/SAIGEQTL 2>/dev/null

# Search for SAIGEQTL installation
find /usr -name SAIGEQTL -type d 2>/dev/null | head -5
find ~ -name SAIGEQTL -type d 2>/dev/null | head -5
```

**Set library path and run validation:**
```bash
# Method 1: Export environment variable
export SAIGEQTL_LIBRARY_PATH=/path/to/your/R/library
./testing/run_regression_test.sh

# Method 2: Specify directly
SAIGEQTL_LIBRARY_PATH=/path/to/R/library ./testing/run_regression_test.sh
```

### Reset Configuration
```bash
./testing/run_regression_test.sh reconfig
```

### Clean Test Outputs  
```bash
./testing/run_regression_test.sh clean
```

### View Test Options
```bash
./testing/run_regression_test.sh help
```

## Version Information

These reference outputs were generated with SAIGEQTL version 0.3.2+ with optimization features.

Key features validated:
- Multi-phenotype preprocessing optimization
- Cross-format consistency (PLINK/VCF/BGEN)
- Variance ratio estimation
- SPA (Saddlepoint Approximation) correction
- Mixed model convergence

## Test Data

The validation uses simulated data included in `extdata/input/`:
- 100 individuals, 100 cells per individual  
- Single gene expression phenotype (gene_1)
- ~300 genetic variants on chromosome 2
- Covariates: X1, X2 (sample-level), pf1, pf2 (cell-level)

This represents a realistic but small-scale eQTL analysis scenario for testing purposes.