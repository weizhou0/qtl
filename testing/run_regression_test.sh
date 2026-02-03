#!/bin/bash

# SAIGEQTL Package Regression Test Runner
# This script provides a convenient interface for running regression tests

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PACKAGE_ROOT"

CONFIG_FILE="$PACKAGE_ROOT/.saigeqtl_config"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored messages
print_message() {
    local color=$1
    local message=$2
    echo -e "${color}${message}${NC}"
}

# Function to detect SAIGEQTL installation
detect_saigeqtl_library() {
    local library_path=""
    
    # First try using R to find SAIGEQTL in the current .libPaths()
    if command -v pixi &> /dev/null && [[ -f "$SAIGEQTL_PIXI_MANIFEST" ]]; then
        library_path=$(pixi run --manifest-path "$SAIGEQTL_PIXI_MANIFEST" Rscript -e "
            paths <- .libPaths()
            for(p in paths) {
                if(file.exists(file.path(p, 'SAIGEQTL'))) {
                    cat(p)
                    break
                }
            }
        " 2>/dev/null | head -n1)
    fi
    
    # If pixi not available or didn't find it, try with system R
    if [[ -z "$library_path" ]] && command -v Rscript &> /dev/null; then
        library_path=$(Rscript -e "
            paths <- .libPaths()
            for(p in paths) {
                if(file.exists(file.path(p, 'SAIGEQTL'))) {
                    cat(p)
                    break
                }
            }
        " 2>/dev/null | head -n1)
    fi
    
    # If still not found, try common installation locations (but no hardcoded paths)
    if [[ -z "$library_path" ]]; then
        local common_paths=(
            "$HOME/R/library"
            "/usr/local/lib/R/site-library"
            "/usr/lib/R/site-library" 
            "/usr/lib/R/library"
        )
        
        for path in "${common_paths[@]}"; do
            if [[ -d "$path/SAIGEQTL" ]]; then
                library_path="$path"
                break
            fi
        done
    fi
    
    echo "$library_path"
}

# Function to validate library path
validate_library_path() {
    local library_path="$1"
    
    if [[ -z "$library_path" ]]; then
        print_message "$RED" "Library path cannot be empty"
        return 1
    fi
    
    if [[ ! -d "$library_path" ]]; then
        print_message "$RED" "Library path does not exist: $library_path"
        return 1
    fi
    
    if [[ ! -d "$library_path/SAIGEQTL" ]]; then
        print_message "$RED" "SAIGEQTL package not found in: $library_path"
        print_message "$YELLOW" "Expected to find: $library_path/SAIGEQTL"
        return 1
    fi
    
    # Check if DESCRIPTION file exists (confirms it's a real R package)
    if [[ ! -f "$library_path/SAIGEQTL/DESCRIPTION" ]]; then
        print_message "$RED" "SAIGEQTL directory found but appears incomplete: $library_path/SAIGEQTL"
        print_message "$YELLOW" "Missing DESCRIPTION file - this may not be a properly installed R package"
        return 1
    fi
    
    print_message "$GREEN" "✓ Valid SAIGEQTL installation found at: $library_path"
    return 0
}

# Function to save config
save_config() {
    local library_path="$1"
    cat > "$CONFIG_FILE" << EOF
# SAIGEQTL Configuration File
# Auto-generated on $(date)
SAIGEQTL_LIBRARY_PATH="$library_path"
EOF
    print_message "$GREEN" "Configuration saved to $CONFIG_FILE"
}

# Function to load config
load_config() {
    if [[ -f "$CONFIG_FILE" ]]; then
        source "$CONFIG_FILE"
        return 0
    fi
    return 1
}

# Function to setup library path
setup_library_path() {
    # First, try to load from config file
    if load_config && [[ -n "$SAIGEQTL_LIBRARY_PATH" ]] && [[ -d "$SAIGEQTL_LIBRARY_PATH/SAIGEQTL" ]]; then
        print_message "$GREEN" "Using library path from config: $SAIGEQTL_LIBRARY_PATH"
        export SAIGEQTL_LIBRARY_PATH
        return 0
    fi
    
    # If not in config or config is invalid, try to detect
    print_message "$YELLOW" "Detecting SAIGEQTL installation..."
    local detected_path=$(detect_saigeqtl_library)
    
    if [[ -n "$detected_path" ]] && validate_library_path "$detected_path"; then
        export SAIGEQTL_LIBRARY_PATH="$detected_path"
        save_config "$detected_path"
        return 0
    else
        print_message "$RED" "Could not automatically detect SAIGEQTL installation."
        echo ""
        print_message "$YELLOW" "To find your SAIGEQTL installation, try one of these methods:"
        echo ""
        echo "1. Check where R packages are installed:"
        echo "   R -e '.libPaths()'"
        echo ""
        echo "2. Look for SAIGEQTL directory in common locations:"
        echo "   ls -d ~/R/library/SAIGEQTL 2>/dev/null"
        echo "   ls -d /usr/local/lib/R/site-library/SAIGEQTL 2>/dev/null"
        echo ""
        echo "3. Search for SAIGEQTL installation:"
        echo "   find /usr -name SAIGEQTL -type d 2>/dev/null | head -5"
        echo "   find ~ -name SAIGEQTL -type d 2>/dev/null | head -5"
        echo ""
        print_message "$YELLOW" "Once you find the library path, set it and run again:"
        echo "   export SAIGEQTL_LIBRARY_PATH=/path/to/your/R/library"
        echo "   $0"
        echo ""
        print_message "$YELLOW" "Or specify it directly:"
        echo "   SAIGEQTL_LIBRARY_PATH=/path/to/R/library $0"
        echo ""
        return 1
    fi
}

# Default configuration - can be overridden by environment variables
export SAIGEQTL_PACKAGE_ROOT="${SAIGEQTL_PACKAGE_ROOT:-$PACKAGE_ROOT}"
export SAIGEQTL_PIXI_MANIFEST="${SAIGEQTL_PIXI_MANIFEST:-$PACKAGE_ROOT/pixi.toml}"
export SAIGEQTL_TEST_DIR="${SAIGEQTL_TEST_DIR:-test_output}"
export SAIGEQTL_EXTDATA_DIR="${SAIGEQTL_EXTDATA_DIR:-extdata}"

# Set expected output directory - prefer extdata/expected_output if it exists, fall back to expected_output
if [[ -d "$SAIGEQTL_EXTDATA_DIR/expected_output" ]]; then
    export SAIGEQTL_EXPECTED_DIR="${SAIGEQTL_EXPECTED_DIR:-$SAIGEQTL_EXTDATA_DIR/expected_output}"
else
    export SAIGEQTL_EXPECTED_DIR="${SAIGEQTL_EXPECTED_DIR:-expected_output}"
fi

# Setup library path automatically
if [[ -z "$SAIGEQTL_LIBRARY_PATH" ]]; then
    setup_library_path || exit 1
else
    # Validate user-provided library path
    if ! validate_library_path "$SAIGEQTL_LIBRARY_PATH"; then
        print_message "$RED" "Invalid library path provided: $SAIGEQTL_LIBRARY_PATH"
        exit 1
    fi
fi

# Function to show usage
show_usage() {
    echo "SAIGEQTL Regression Test Runner"
    echo "Usage: $0 [OPTION]"
    echo ""
    echo "Options:"
    echo "  test          Run PLINK format regression test (default)"
    echo "  vcf           Run VCF format test only"
    echo "  bgen          Run BGEN format test only"
    echo "  comprehensive Run all format tests (PLINK, VCF, BGEN) with comparison"
    echo "  all           Same as comprehensive"
    echo "  compare       Compare existing Step 2 outputs across formats"
    echo "  clean         Clean test output directory"
    echo "  baseline      Save current test results as baseline"
    echo "  reconfig      Reconfigure SAIGEQTL library path"
    echo "  help          Show this help message"
    echo "  validate      Validate package installation (same as 'test')"
    echo ""
    echo "Environment Variables (optional):"
    echo "  SAIGEQTL_PACKAGE_ROOT   - Path to SAIGEQTL package root (default: script directory)"
    echo "  SAIGEQTL_LIBRARY_PATH   - R library path for package installation (auto-detected if not set)"
    echo "  SAIGEQTL_PIXI_MANIFEST  - Path to pixi.toml file (default: ./pixi.toml)"
    echo ""
    echo "Examples:"
    echo "  $0                                    # Run PLINK test"
    echo "  $0 test                               # Run PLINK test"
    echo "  $0 vcf                                # Run VCF test only (uses GT field)"
    echo "  $0 bgen                               # Run BGEN test only (uses dosages)"
    echo "  $0 comprehensive                      # Run all format tests + compare results"
    echo "  $0 compare                            # Compare existing outputs only"
    echo "  $0 clean                              # Clean test outputs"
    echo "  $0 baseline                           # Save baseline results"
    echo "  $0 reconfig                           # Reconfigure library path"
    echo "  $0 validate                           # Validate installation"
    echo ""
    echo "  # Custom library path:"
    echo "  SAIGEQTL_LIBRARY_PATH=/custom/path $0 comprehensive"
}

# Function to check prerequisites
check_prerequisites() {
    print_message "$YELLOW" "Checking prerequisites..."
    
    # Check if pixi is available
    if ! command -v pixi &> /dev/null; then
        print_message "$RED" "Error: pixi not found. Please install pixi."
        exit 1
    fi
    
    # Check if pixi.toml exists
    if [[ ! -f "$SAIGEQTL_PIXI_MANIFEST" ]]; then
        print_message "$RED" "Error: pixi.toml not found at $SAIGEQTL_PIXI_MANIFEST"
        exit 1
    fi
    
    # Check if required R scripts exist
    if [[ ! -f "$SAIGEQTL_EXTDATA_DIR/step1_fitNULLGLMM_qtl.R" ]]; then
        print_message "$RED" "Error: step1_fitNULLGLMM_qtl.R not found in $SAIGEQTL_EXTDATA_DIR/"
        exit 1
    fi
    
    if [[ ! -f "$SAIGEQTL_EXTDATA_DIR/step2_tests_qtl.R" ]]; then
        print_message "$RED" "Error: step2_tests_qtl.R not found in $SAIGEQTL_EXTDATA_DIR/"
        exit 1
    fi
    
    # Check if test data exists
    if [[ ! -d "$SAIGEQTL_EXTDATA_DIR/input" ]]; then
        print_message "$RED" "Error: Test data directory $SAIGEQTL_EXTDATA_DIR/input not found"
        exit 1
    fi
    
    # Check if test R script exists
    if [[ ! -f "testing/test_package_regression.R" ]]; then
        print_message "$RED" "Error: testing/test_package_regression.R not found"
        exit 1
    fi
    
    # Verify SAIGEQTL package is available
    if [[ -n "$SAIGEQTL_LIBRARY_PATH" ]] && [[ -d "$SAIGEQTL_LIBRARY_PATH/SAIGEQTL" ]]; then
        print_message "$GREEN" "SAIGEQTL package found at: $SAIGEQTL_LIBRARY_PATH"
    else
        print_message "$RED" "Error: SAIGEQTL package not found at: $SAIGEQTL_LIBRARY_PATH"
        exit 1
    fi
    
    print_message "$GREEN" "Prerequisites check passed"
}

# Function to run the regression test
run_test() {
    print_message "$GREEN" "Starting PLINK format regression test..."
    echo "=========================================="
    
    check_prerequisites
    
    # Show configuration
    print_message "$YELLOW" "Configuration:"
    echo "  Package Root: $SAIGEQTL_PACKAGE_ROOT"
    echo "  Pixi Manifest: $SAIGEQTL_PIXI_MANIFEST"
    echo "  Library Path: ${SAIGEQTL_LIBRARY_PATH:-<not set>}"
    echo "  Test Dir: $SAIGEQTL_TEST_DIR"
    echo ""
    
    # Run shared Step 1
    if ! run_shared_step1; then
        exit 1
    fi
    
    # Step 2: PLINK association tests  
    print_message "$YELLOW" "Step 2: Running PLINK format association tests"
    
    # Create region file for cis eQTL analysis
    region_file="$SAIGEQTL_TEST_DIR/gene_1_cis_region.txt"
    echo -e "2\t1\t9810000" > "$region_file"
    
    step2_args="--bedFile=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.bed --bimFile=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.bim --famFile=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.fam --SAIGEOutputFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_plink_cis --chrom=2 --minMAF=0 --minMAC=20 --LOCO=FALSE --GMMATmodelFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared.rda --SPAcutoff=2 --varianceRatioFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared.varianceRatio.txt --rangestoIncludeFile=$region_file --markers_per_chunk=10000"
    
    # Build step2 command with library path
    step2_cmd="pixi run --manifest-path $SAIGEQTL_PIXI_MANIFEST Rscript $SAIGEQTL_EXTDATA_DIR/step2_tests_qtl.R"
    if [[ -n "$SAIGEQTL_LIBRARY_PATH" ]]; then
        step2_cmd="$step2_cmd --library=$SAIGEQTL_LIBRARY_PATH"
    fi
    step2_cmd="$step2_cmd $step2_args"
    
    if eval $step2_cmd; then
        print_message "$GREEN" "✓ PLINK format test completed successfully!"
    else
        print_message "$RED" "✗ PLINK format test failed!"
        exit 1
    fi
}

# Function to clean test outputs
clean_outputs() {
    print_message "$YELLOW" "Cleaning test outputs..."
    pixi run --manifest-path "$SAIGEQTL_PIXI_MANIFEST" Rscript testing/test_package_regression.R clean
    print_message "$GREEN" "✓ Test outputs cleaned"
}

# Function to save baseline
save_baseline() {
    print_message "$YELLOW" "Saving current results as baseline..."
    pixi run --manifest-path "$SAIGEQTL_PIXI_MANIFEST" Rscript testing/test_package_regression.R baseline
    print_message "$GREEN" "✓ Baseline saved"
}

# Function to reconfigure library path
reconfigure_library() {
    print_message "$YELLOW" "Reconfiguring SAIGEQTL library path..."
    
    if [[ -f "$CONFIG_FILE" ]]; then
        rm "$CONFIG_FILE"
        print_message "$YELLOW" "Removed existing configuration"
    fi
    
    # Reset the environment variable
    unset SAIGEQTL_LIBRARY_PATH
    
    # Force re-detection
    setup_library_path || {
        print_message "$RED" "Failed to auto-detect SAIGEQTL installation"
        print_message "$YELLOW" "Please install SAIGEQTL or set SAIGEQTL_LIBRARY_PATH manually"
        exit 1
    }
    
    print_message "$GREEN" "✓ Library path reconfigured"
}

# Shared Step 1 function - runs NULL GLMM fitting once for all formats
run_shared_step1() {
    local step1_prefix="$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared"
    
    # Check if Step 1 outputs already exist
    if [[ -f "$step1_prefix.rda" ]] && [[ -f "$step1_prefix.varianceRatio.txt" ]]; then
        print_message "$GREEN" "✓ Step 1 outputs already exist, skipping..."
        echo "  - Model file: $step1_prefix.rda"
        echo "  - Variance ratio file: $step1_prefix.varianceRatio.txt"
        return 0
    fi
    
    print_message "$YELLOW" "Step 1: Fitting shared NULL GLMM for all formats"
    step1_args="--useSparseGRMtoFitNULL=FALSE --useGRMtoFitNULL=FALSE --phenoFile=./$SAIGEQTL_EXTDATA_DIR/input/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt --phenoCol=gene_1 --covarColList=X1,X2,pf1,pf2 --sampleCovarColList=X1,X2 --sampleIDColinphenoFile=IND_ID --traitType=count --outputPrefix=$step1_prefix --skipVarianceRatioEstimation=FALSE --isRemoveZerosinPheno=FALSE --isCovariateOffset=FALSE --isCovariateTransform=TRUE --skipModelFitting=FALSE --tol=0.00001 --plinkFile=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1_01.step1 --IsOverwriteVarianceRatioFile=TRUE"
    
    # Build step1 command with library path
    step1_cmd="pixi run --manifest-path $SAIGEQTL_PIXI_MANIFEST Rscript $SAIGEQTL_EXTDATA_DIR/step1_fitNULLGLMM_qtl.R"
    if [[ -n "$SAIGEQTL_LIBRARY_PATH" ]]; then
        step1_cmd="$step1_cmd --library=$SAIGEQTL_LIBRARY_PATH"
    fi
    step1_cmd="$step1_cmd $step1_args"
    
    if eval $step1_cmd; then
        print_message "$GREEN" "✓ Shared Step 1 (NULL GLMM) completed successfully!"
        return 0
    else
        print_message "$RED" "✗ Shared Step 1 (NULL GLMM) failed!"
        return 1
    fi
}

# Function to run VCF format test
run_vcf_test() {
    print_message "$GREEN" "Starting VCF format regression test..."
    echo "=========================================="
    
    check_prerequisites
    
    # Show configuration
    print_message "$YELLOW" "Configuration:"
    echo "  Package Root: $SAIGEQTL_PACKAGE_ROOT"
    echo "  Pixi Manifest: $SAIGEQTL_PIXI_MANIFEST"
    echo "  Library Path: ${SAIGEQTL_LIBRARY_PATH:-<not set>}"
    echo "  Test Dir: $SAIGEQTL_TEST_DIR"
    echo ""
    
    # Run shared Step 1
    if ! run_shared_step1; then
        exit 1
    fi
    
    # Step 2: VCF association tests
    print_message "$YELLOW" "Step 2: Running VCF format association tests"
    
    # Create region file for cis eQTL analysis
    region_file="$SAIGEQTL_TEST_DIR/gene_1_cis_region_vcf.txt"
    echo -e "2\t300001\t610001" > "$region_file"
    
    step2_args="--vcfFile=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.vcf.gz --vcfFileIndex=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.vcf.gz.csi --vcfField=GT --chrom=2 --SAIGEOutputFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_vcf_cis --minMAF=0 --minMAC=20 --LOCO=FALSE --GMMATmodelFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared.rda --SPAcutoff=2 --varianceRatioFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared.varianceRatio.txt --rangestoIncludeFile=$region_file --markers_per_chunk=10000"
    
    # Build step2 command with library path
    step2_cmd="pixi run --manifest-path $SAIGEQTL_PIXI_MANIFEST Rscript $SAIGEQTL_EXTDATA_DIR/step2_tests_qtl.R"
    if [[ -n "$SAIGEQTL_LIBRARY_PATH" ]]; then
        step2_cmd="$step2_cmd --library=$SAIGEQTL_LIBRARY_PATH"
    fi
    step2_cmd="$step2_cmd $step2_args"
    
    if eval $step2_cmd; then
        print_message "$GREEN" "✓ VCF format test completed successfully!"
    else
        print_message "$RED" "✗ VCF format test failed!"
        exit 1
    fi
}

# Function to run BGEN format test
run_bgen_test() {
    print_message "$GREEN" "Starting BGEN format regression test..."
    echo "=========================================="
    
    check_prerequisites
    
    # Show configuration
    print_message "$YELLOW" "Configuration:"
    echo "  Package Root: $SAIGEQTL_PACKAGE_ROOT"
    echo "  Pixi Manifest: $SAIGEQTL_PIXI_MANIFEST"
    echo "  Library Path: ${SAIGEQTL_LIBRARY_PATH:-<not set>}"
    echo "  Test Dir: $SAIGEQTL_TEST_DIR"
    echo ""
    
    # Run shared Step 1
    if ! run_shared_step1; then
        exit 1
    fi
    
    # Step 2: BGEN association tests
    print_message "$YELLOW" "Step 2: Running BGEN format association tests"
    step2_args="--bgenFile=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.bgen --bgenFileIndex=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.bgen.bgi --sampleFile=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.sample --AlleleOrder=ref-first --chrom=2 --SAIGEOutputFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_bgen_cis --minMAF=0 --minMAC=20 --LOCO=FALSE --GMMATmodelFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared.rda --SPAcutoff=2 --varianceRatioFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared.varianceRatio.txt --markers_per_chunk=10000"
    
    # Build step2 command with library path
    step2_cmd="pixi run --manifest-path $SAIGEQTL_PIXI_MANIFEST Rscript $SAIGEQTL_EXTDATA_DIR/step2_tests_qtl.R"
    if [[ -n "$SAIGEQTL_LIBRARY_PATH" ]]; then
        step2_cmd="$step2_cmd --library=$SAIGEQTL_LIBRARY_PATH"
    fi
    step2_cmd="$step2_cmd $step2_args"
    
    if eval $step2_cmd; then
        print_message "$GREEN" "✓ BGEN format test completed successfully!"
    else
        print_message "$RED" "✗ BGEN format test failed!"
        exit 1
    fi
}

# Function to compare Step 2 outputs across formats
compare_format_outputs() {
    print_message "$GREEN" "Comparing Step 2 outputs across formats..."
    echo "============================================"
    
    local plink_output="$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_plink_cis"
    local vcf_output="$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_vcf_cis"
    local bgen_output="$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_bgen_cis"
    
    # Check if all output files exist
    local missing_files=()
    if [[ ! -f "$plink_output" ]]; then missing_files+=("PLINK"); fi
    if [[ ! -f "$vcf_output" ]]; then missing_files+=("VCF"); fi
    if [[ ! -f "$bgen_output" ]]; then missing_files+=("BGEN"); fi
    
    if [[ ${#missing_files[@]} -gt 0 ]]; then
        print_message "$YELLOW" "Warning: Some output files missing for comparison:"
        for format in "${missing_files[@]}"; do
            print_message "$YELLOW" "  - $format format output not found"
        done
        return 1
    fi
    
    print_message "$YELLOW" "Output files found:"
    echo "  - PLINK: $plink_output"
    echo "  - VCF:   $vcf_output"
    echo "  - BGEN:  $bgen_output"
    echo ""
    
    # Basic file size comparison
    print_message "$YELLOW" "File size comparison:"
    echo "  - PLINK: $(wc -l < "$plink_output") lines"
    echo "  - VCF:   $(wc -l < "$vcf_output") lines"
    echo "  - BGEN:  $(wc -l < "$bgen_output") lines"
    echo ""
    
    # Create R script for detailed comparison
    local comparison_script="$SAIGEQTL_TEST_DIR/compare_outputs.R"
    cat > "$comparison_script" << 'EOF'
#!/usr/bin/env Rscript

# SAIGEQTL Output Comparison Script

# Handle library path from environment
library_path <- Sys.getenv("SAIGEQTL_LIBRARY_PATH", "")
if (library_path != "") {
  .libPaths(c(library_path, .libPaths()))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: Rscript compare_outputs.R <plink_file> <vcf_file> <bgen_file>")
}

plink_file <- args[1]
vcf_file <- args[2]
bgen_file <- args[3]

cat("Reading output files...\n")

# Read the output files
plink_data <- tryCatch({
    read.table(plink_file, header = TRUE, stringsAsFactors = FALSE)
}, error = function(e) {
    cat("Error reading PLINK file:", e$message, "\n")
    return(NULL)
})

vcf_data <- tryCatch({
    read.table(vcf_file, header = TRUE, stringsAsFactors = FALSE)
}, error = function(e) {
    cat("Error reading VCF file:", e$message, "\n")
    return(NULL)
})

bgen_data <- tryCatch({
    read.table(bgen_file, header = TRUE, stringsAsFactors = FALSE)
}, error = function(e) {
    cat("Error reading BGEN file:", e$message, "\n")
    return(NULL)
})

if (is.null(plink_data) || is.null(vcf_data) || is.null(bgen_data)) {
    stop("Failed to read one or more input files")
}

cat("Comparing outputs...\n\n")

# Basic statistics
cat("=== BASIC STATISTICS ===\n")
cat(sprintf("PLINK: %d variants\n", nrow(plink_data)))
cat(sprintf("VCF:   %d variants\n", nrow(vcf_data)))
cat(sprintf("BGEN:  %d variants\n", nrow(bgen_data)))
cat("\n")

# Column names comparison
cat("=== COLUMN COMPARISON ===\n")
cat("PLINK columns:", paste(colnames(plink_data), collapse = ", "), "\n")
cat("VCF columns:  ", paste(colnames(vcf_data), collapse = ", "), "\n")
cat("BGEN columns: ", paste(colnames(bgen_data), collapse = ", "), "\n")
cat("\n")

# Find common variants (assuming CHR:POS or similar identifier)
# Try different potential ID columns
id_cols <- c("MarkerID", "SNPID", "CHR_POS", "rsID", "ID")
id_col <- NULL

for (col in id_cols) {
    if (col %in% colnames(plink_data) && col %in% colnames(vcf_data) && col %in% colnames(bgen_data)) {
        id_col <- col
        break
    }
}

if (is.null(id_col)) {
    # Try to create CHR:POS identifier if CHR and POS columns exist
    if ("CHR" %in% colnames(plink_data) && "POS" %in% colnames(plink_data)) {
        plink_data$ID <- paste(plink_data$CHR, plink_data$POS, sep = ":")
        vcf_data$ID <- paste(vcf_data$CHR, vcf_data$POS, sep = ":")
        bgen_data$ID <- paste(bgen_data$CHR, bgen_data$POS, sep = ":")
        id_col <- "ID"
    } else {
        cat("Warning: No common identifier column found. Using row order for comparison.\n")
        cat("This may not be accurate if variants are in different orders.\n\n")
        
        # Show first few rows for manual inspection
        cat("=== FIRST 5 ROWS COMPARISON ===\n")
        cat("PLINK (first 5 rows):\n")
        print(head(plink_data, 5))
        cat("\nVCF (first 5 rows):\n")
        print(head(vcf_data, 5))
        cat("\nBGEN (first 5 rows):\n")
        print(head(bgen_data, 5))
        
        return()
    }
}

# Compare common variants
cat(sprintf("Using '%s' as identifier column\n", id_col))

common_plink_vcf <- intersect(plink_data[[id_col]], vcf_data[[id_col]])
common_plink_bgen <- intersect(plink_data[[id_col]], bgen_data[[id_col]])
common_vcf_bgen <- intersect(vcf_data[[id_col]], bgen_data[[id_col]])
common_all <- intersect(common_plink_vcf, common_plink_bgen)

cat(sprintf("Common variants PLINK-VCF:  %d\n", length(common_plink_vcf)))
cat(sprintf("Common variants PLINK-BGEN: %d\n", length(common_plink_bgen)))
cat(sprintf("Common variants VCF-BGEN:   %d\n", length(common_vcf_bgen)))
cat(sprintf("Common variants (all):      %d\n", length(common_all)))
cat("\n")

if (length(common_all) == 0) {
    cat("No common variants found across all formats!\n")
    return()
}

# Subset to common variants and compare key statistics
plink_common <- plink_data[plink_data[[id_col]] %in% common_all, ]
vcf_common <- vcf_data[vcf_data[[id_col]] %in% common_all, ]
bgen_common <- bgen_data[bgen_data[[id_col]] %in% common_all, ]

# Sort by ID for proper comparison
plink_common <- plink_common[order(plink_common[[id_col]]), ]
vcf_common <- vcf_common[order(vcf_common[[id_col]]), ]
bgen_common <- bgen_common[order(bgen_common[[id_col]]), ]

# Compare p-values if available
p_cols <- c("p.value", "Pvalue", "P", "pval")
p_col <- NULL
for (col in p_cols) {
    if (col %in% colnames(plink_common)) {
        p_col <- col
        break
    }
}

if (!is.null(p_col)) {
    cat("=== P-VALUE COMPARISON ===\n")
    
    plink_p <- plink_common[[p_col]]
    vcf_p <- vcf_common[[p_col]]
    bgen_p <- bgen_common[[p_col]]
    
    # Correlation between p-values
    cor_plink_vcf <- cor(plink_p, vcf_p, use = "complete.obs")
    cor_plink_bgen <- cor(plink_p, bgen_p, use = "complete.obs")
    cor_vcf_bgen <- cor(vcf_p, bgen_p, use = "complete.obs")
    
    cat(sprintf("P-value correlation PLINK-VCF:  %.6f\n", cor_plink_vcf))
    cat(sprintf("P-value correlation PLINK-BGEN: %.6f\n", cor_plink_bgen))
    cat(sprintf("P-value correlation VCF-BGEN:   %.6f\n", cor_vcf_bgen))
    
    # Max differences
    max_diff_plink_vcf <- max(abs(plink_p - vcf_p), na.rm = TRUE)
    max_diff_plink_bgen <- max(abs(plink_p - bgen_p), na.rm = TRUE)
    max_diff_vcf_bgen <- max(abs(vcf_p - bgen_p), na.rm = TRUE)
    
    cat(sprintf("Max p-value diff PLINK-VCF:  %.2e\n", max_diff_plink_vcf))
    cat(sprintf("Max p-value diff PLINK-BGEN: %.2e\n", max_diff_plink_bgen))
    cat(sprintf("Max p-value diff VCF-BGEN:   %.2e\n", max_diff_vcf_bgen))
    cat("\n")
}

# Compare effect sizes if available
beta_cols <- c("BETA", "beta", "Beta", "Effect")
beta_col <- NULL
for (col in beta_cols) {
    if (col %in% colnames(plink_common)) {
        beta_col <- col
        break
    }
}

if (!is.null(beta_col)) {
    cat("=== EFFECT SIZE COMPARISON ===\n")
    
    plink_beta <- plink_common[[beta_col]]
    vcf_beta <- vcf_common[[beta_col]]
    bgen_beta <- bgen_common[[beta_col]]
    
    # Correlation between effect sizes
    cor_plink_vcf <- cor(plink_beta, vcf_beta, use = "complete.obs")
    cor_plink_bgen <- cor(plink_beta, bgen_beta, use = "complete.obs")
    cor_vcf_bgen <- cor(vcf_beta, bgen_beta, use = "complete.obs")
    
    cat(sprintf("Effect size correlation PLINK-VCF:  %.6f\n", cor_plink_vcf))
    cat(sprintf("Effect size correlation PLINK-BGEN: %.6f\n", cor_plink_bgen))
    cat(sprintf("Effect size correlation VCF-BGEN:   %.6f\n", cor_vcf_bgen))
    cat("\n")
}

cat("=== COMPARISON SUMMARY ===\n")
if (!is.null(p_col) && all(c(cor_plink_vcf, cor_plink_bgen, cor_vcf_bgen) > 0.999, na.rm = TRUE)) {
    cat("✓ EXCELLENT: All formats produce highly consistent results (r > 0.999)\n")
} else if (!is.null(p_col) && all(c(cor_plink_vcf, cor_plink_bgen, cor_vcf_bgen) > 0.99, na.rm = TRUE)) {
    cat("✓ GOOD: All formats produce consistent results (r > 0.99)\n")
} else if (!is.null(p_col)) {
    cat("⚠ WARNING: Some differences detected between formats\n")
} else {
    cat("? Unable to assess consistency - no p-value column found\n")
}

cat("Analysis complete.\n")
EOF

    # Run the R comparison script
    print_message "$YELLOW" "Running detailed statistical comparison..."
    
    # Build comparison command with library path
    if [[ -n "$SAIGEQTL_LIBRARY_PATH" ]]; then
        comparison_cmd="SAIGEQTL_LIBRARY_PATH=$SAIGEQTL_LIBRARY_PATH pixi run --manifest-path $SAIGEQTL_PIXI_MANIFEST Rscript $comparison_script $plink_output $vcf_output $bgen_output"
    else
        comparison_cmd="pixi run --manifest-path $SAIGEQTL_PIXI_MANIFEST Rscript $comparison_script $plink_output $vcf_output $bgen_output"
    fi
    
    if eval $comparison_cmd; then
        print_message "$GREEN" "✓ Output comparison completed successfully!"
    else
        print_message "$RED" "✗ Output comparison failed!"
        return 1
    fi
    
    return 0
}

# Function to run comprehensive test (all formats)
run_comprehensive_test() {
    print_message "$GREEN" "Starting comprehensive SAIGEQTL format testing..."
    echo "================================================="
    
    check_prerequisites
    
    # Show configuration
    print_message "$YELLOW" "Configuration:"
    echo "  Package Root: $SAIGEQTL_PACKAGE_ROOT"
    echo "  Pixi Manifest: $SAIGEQTL_PIXI_MANIFEST"
    echo "  Library Path: ${SAIGEQTL_LIBRARY_PATH:-<not set>}"
    echo "  Test Dir: $SAIGEQTL_TEST_DIR"
    echo ""
    
    # Run shared Step 1 once for all formats
    print_message "$GREEN" "Phase 1: Running shared NULL GLMM fitting for all formats..."
    if ! run_shared_step1; then
        exit 1
    fi
    echo ""
    
    # Track results for each format
    local all_passed=true
    
    print_message "$GREEN" "Phase 2: Running association tests for each format..."
    echo ""
    
    # PLINK format test (Step 2 only)
    print_message "$YELLOW" "Testing PLINK format..."
    region_file_plink="$SAIGEQTL_TEST_DIR/gene_1_cis_region_plink.txt"
    echo "2	1	9810000" > "$region_file_plink"
    
    step2_plink_args="--bedFile=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.bed --bimFile=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.bim --famFile=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.fam --SAIGEOutputFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_plink_cis --chrom=2 --minMAF=0 --minMAC=20 --LOCO=FALSE --GMMATmodelFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared.rda --SPAcutoff=2 --varianceRatioFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared.varianceRatio.txt --rangestoIncludeFile=$region_file_plink --markers_per_chunk=10000"
    
    step2_cmd="pixi run --manifest-path $SAIGEQTL_PIXI_MANIFEST Rscript $SAIGEQTL_EXTDATA_DIR/step2_tests_qtl.R"
    if [[ -n "$SAIGEQTL_LIBRARY_PATH" ]]; then
        step2_cmd="$step2_cmd --library=$SAIGEQTL_LIBRARY_PATH"
    fi
    
    if eval "$step2_cmd $step2_plink_args"; then
        print_message "$GREEN" "✓ PLINK format test passed"
    else
        print_message "$RED" "✗ PLINK format test failed"
        all_passed=false
    fi
    echo ""
    
    # VCF format test (Step 2 only)
    print_message "$YELLOW" "Testing VCF format..."
    step2_vcf_args="--vcfFile=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.vcf.gz --vcfFileIndex=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.vcf.gz.csi --vcfField=GT --chrom=2 --SAIGEOutputFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_vcf_cis --minMAF=0 --minMAC=20 --LOCO=FALSE --GMMATmodelFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared.rda --SPAcutoff=2 --varianceRatioFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared.varianceRatio.txt --markers_per_chunk=10000"
    
    if eval "$step2_cmd $step2_vcf_args"; then
        print_message "$GREEN" "✓ VCF format test passed"
    else
        print_message "$RED" "✗ VCF format test failed"
        all_passed=false
    fi
    echo ""
    
    # BGEN format test (Step 2 only)
    print_message "$YELLOW" "Testing BGEN format..."
    step2_bgen_args="--bgenFile=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.bgen --bgenFileIndex=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.bgen.bgi --sampleFile=./$SAIGEQTL_EXTDATA_DIR/input/n.indep_100_n.cell_1.sample --AlleleOrder=ref-first --chrom=2 --SAIGEOutputFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_bgen_cis --minMAF=0 --minMAC=20 --LOCO=FALSE --GMMATmodelFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared.rda --SPAcutoff=2 --varianceRatioFile=./$SAIGEQTL_TEST_DIR/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_shared.varianceRatio.txt --markers_per_chunk=10000"
    
    if eval "$step2_cmd $step2_bgen_args"; then
        print_message "$GREEN" "✓ BGEN format test passed"
    else
        print_message "$RED" "✗ BGEN format test failed"
        all_passed=false
    fi
    echo ""
    
    # Summary
    print_message "$YELLOW" "Comprehensive Test Summary:"
    print_message "$YELLOW" "=========================="
    if $all_passed; then
        print_message "$GREEN" "✓ All format tests completed successfully!"
        print_message "$GREEN" "  - PLINK format: PASSED"
        print_message "$GREEN" "  - VCF format: PASSED"
        print_message "$GREEN" "  - BGEN format: PASSED"
        print_message "$YELLOW" "Optimization: Step 1 (NULL GLMM) was run only once and shared across all formats"
        echo ""
        
        # Run comparison analysis
        print_message "$GREEN" "Phase 3: Comparing outputs across formats..."
        if compare_format_outputs; then
            print_message "$GREEN" "✓ Cross-format comparison completed!"
        else
            print_message "$YELLOW" "⚠ Cross-format comparison had issues (see details above)"
        fi
    else
        print_message "$RED" "✗ Some format tests failed!"
        exit 1
    fi
}

# Function to run specific test
run_specific_test() {
    local test_type=$1
    print_message "$GREEN" "Starting SAIGEQTL $test_type test..."
    echo "=========================================="
    
    check_prerequisites
    
    # Show configuration
    print_message "$YELLOW" "Configuration:"
    echo "  Package Root: $SAIGEQTL_PACKAGE_ROOT"
    echo "  Pixi Manifest: $SAIGEQTL_PIXI_MANIFEST"
    echo "  Library Path: ${SAIGEQTL_LIBRARY_PATH:-<not set>}"
    echo "  Test Dir: $SAIGEQTL_TEST_DIR"
    echo ""
    
    # Run the R test script with specific argument
    if pixi run --manifest-path "$SAIGEQTL_PIXI_MANIFEST" Rscript testing/test_package_regression.R "$test_type"; then
        print_message "$GREEN" "✓ $test_type test completed successfully!"
    else
        print_message "$RED" "✗ $test_type test failed!"
        exit 1
    fi
}

# Main execution
case "${1:-test}" in
    "test"|"validate")
        run_test
        ;;
    "vcf")
        run_vcf_test
        ;;
    "bgen")
        run_bgen_test
        ;;
    "comprehensive"|"all")
        run_comprehensive_test
        ;;
    "compare")
        compare_format_outputs
        ;;
    "clean")
        clean_outputs
        ;;
    "baseline")
        save_baseline
        ;;
    "reconfig")
        reconfigure_library
        ;;
    "help")
        show_usage
        ;;
    *)
        print_message "$RED" "Unknown option: $1"
        show_usage
        exit 1
        ;;
esac
