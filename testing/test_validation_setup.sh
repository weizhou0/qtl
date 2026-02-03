#!/bin/bash

# Test script to verify validation setup is ready

echo "=== SAIGEQTL Installation Validation Setup Test ==="
echo

# Check required files exist
echo "1. Checking required files..."

# Test data
if [[ -d "extdata/input" ]]; then
    echo "✓ Test input data found"
else
    echo "✗ Test input data missing"
    exit 1
fi

# Expected output
if [[ -d "extdata/expected_output" ]]; then
    echo "✓ Expected output directory found"
    echo "  Files: $(ls extdata/expected_output/ | wc -l) reference files"
else
    echo "✗ Expected output directory missing"
    exit 1
fi

# Scripts
if [[ -f "run_regression_test.sh" ]]; then
    echo "✓ Main test runner found"
else
    echo "✗ Main test runner missing"
    exit 1
fi

if [[ -f "test_package_regression.R" ]]; then
    echo "✓ R test script found"
else
    echo "✗ R test script missing"
    exit 1
fi

# Documentation
if [[ -f "extdata/INSTALLATION_VALIDATION.md" ]]; then
    echo "✓ Validation documentation found"
else
    echo "✗ Validation documentation missing"
    exit 1
fi

echo

echo "2. Checking expected output content..."
EXPECTED_FILES=(
    "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_plink_cis"
    "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_vcf_cis"  
    "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_bgen_cis"
    "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared.varianceRatio.txt"
)

all_found=true
for file in "${EXPECTED_FILES[@]}"; do
    if [[ -f "extdata/expected_output/$file" ]]; then
        size=$(stat -f%z "extdata/expected_output/$file" 2>/dev/null || stat -c%s "extdata/expected_output/$file" 2>/dev/null || echo "unknown")
        echo "✓ $file ($size bytes)"
    else
        echo "✗ $file (missing)"
        all_found=false
    fi
done

if ! $all_found; then
    echo
    echo "✗ Some expected output files are missing!"
    exit 1
fi

echo

echo "3. Testing script configuration..."
# Test the expected output directory detection logic
if [[ -d "extdata/expected_output" ]]; then
    echo "✓ Script will use extdata/expected_output for validation"
else
    echo "⚠ Script will fall back to expected_output directory"
fi

echo

echo "=== Validation Setup Complete ==="
echo
echo "Users can now validate their installation by running:"
echo "  ./run_regression_test.sh validate"
echo
echo "Or for comprehensive testing:"  
echo "  ./run_regression_test.sh comprehensive"
echo
echo "For help:"
echo "  ./run_regression_test.sh help"
echo
echo "Documentation available at:"
echo "  extdata/INSTALLATION_VALIDATION.md"