#!/bin/bash

# SAIGEQTL Test Configuration Script
# Sets up environment variables for testing with custom paths

# Example configuration for your specific setup
export SAIGEQTL_PACKAGE_ROOT="/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/tool_dev/code_optimize/qtl_install_pixi_youngde/qtl"
#export SAIGEQTL_LIBRARY_PATH="/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/tool_dev/installs_0.3.2_beforeoptimization"
export SAIGEQTL_LIBRARY_PATH="/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/tool_dev/installs_0.3.2_afteroptimization_testing"
export SAIGEQTL_PIXI_MANIFEST="$SAIGEQTL_PACKAGE_ROOT/pixi.toml"

# Optional: Override default directories
# export SAIGEQTL_TEST_DIR="custom_test_output"
# export SAIGEQTL_EXPECTED_DIR="custom_expected_output"
# export SAIGEQTL_EXTDATA_DIR="extdata"

echo "SAIGEQTL Test Configuration Set:"
echo "  Package Root: $SAIGEQTL_PACKAGE_ROOT"
echo "  Library Path: $SAIGEQTL_LIBRARY_PATH"
echo "  Pixi Manifest: $SAIGEQTL_PIXI_MANIFEST"
echo ""
echo "Usage:"
echo "  source config_test.sh"
echo "  ./run_regression_test.sh"
echo ""
echo "Or run directly:"
echo "  source config_test.sh && ./run_regression_test.sh"
