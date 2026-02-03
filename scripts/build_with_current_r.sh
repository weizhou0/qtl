#!/bin/bash

# Simple builder using your working pixi environment

set -e

echo "=== SAIGEQTL Builder (Using Current Pixi R) ==="

# Use your existing working pixi environment
OUTPUT_DIR="current_r_binaries"
mkdir -p "$OUTPUT_DIR"

echo "Checking current R version in pixi..."
R_VERSION=$(CONDA_OVERRIDE_GLIBC=2.28 pixi run R --slave -e "cat(paste(R.Version()\$major, R.Version()\$minor, sep='.'))")
R_MAJOR_MINOR=$(CONDA_OVERRIDE_GLIBC=2.28 pixi run R --slave -e "cat(paste(R.Version()\$major, strsplit(R.Version()\$minor, '\\\\.')[[1]][1], sep='.'))")

echo "Found R version: $R_VERSION"
echo "Using R: $R_MAJOR_MINOR"

echo ""
echo "Installing dependencies..."
if CONDA_OVERRIDE_GLIBC=2.28 pixi run install-github-deps; then
    echo "✓ GitHub dependencies installed"
else
    echo "! GitHub dependencies failed (continuing anyway)"
fi

echo ""
echo "Building package with current R ($R_VERSION)..."
if CONDA_OVERRIDE_GLIBC=2.28 pixi run R -e "source('build_binary.R'); main('$OUTPUT_DIR')"; then
    echo "✓ Package built successfully"
    
    # Find and rename the package
    LATEST_PKG=$(ls -t "$OUTPUT_DIR"/SAIGEQTL_*.t*gz 2>/dev/null | head -1)
    if [ -n "$LATEST_PKG" ]; then
        NEW_NAME="$OUTPUT_DIR/SAIGEQTL_0.3.4_R-${R_MAJOR_MINOR}_linux-x86_64_pixi.tar.gz"
        mv "$LATEST_PKG" "$NEW_NAME"
        echo "Package renamed to: $(basename "$NEW_NAME")"
        
        echo ""
        echo "Testing package..."
        if CONDA_OVERRIDE_GLIBC=2.28 pixi run R -e "
            temp_lib <- tempfile('test_')
            dir.create(temp_lib)
            install.packages('$NEW_NAME', lib=temp_lib, repos=NULL, type='source')
            library(SAIGEQTL, lib.loc=temp_lib)
            cat('✓ Package test successful - R', '$R_VERSION', '\\n')
            cat('Version:', as.character(packageVersion('SAIGEQTL')), '\\n')
            unlink(temp_lib, recursive=TRUE)
        "; then
            echo "✓ Package test passed"
        else
            echo "✗ Package test failed"
        fi
    fi
else
    echo "✗ Package build failed"
    exit 1
fi

echo ""
echo "=== BUILD COMPLETE ==="
echo "Package created: $(basename "$NEW_NAME")"
echo "Location: $OUTPUT_DIR/"
echo "Size: $(du -h "$NEW_NAME" | cut -f1)"

echo ""
echo "=== TESTING INSTRUCTIONS ==="
echo "To test this package on different R versions:"
echo ""
echo "# With R 4.3:"
echo "module load R/4.3"
echo "R -e \"install.packages('$NEW_NAME', repos=NULL, type='source'); library(SAIGEQTL)\""
echo ""
echo "# With R 4.4:"
echo "module load R/4.4"  
echo "R -e \"install.packages('$NEW_NAME', repos=NULL, type='source'); library(SAIGEQTL)\""
echo ""
echo "# In pixi environment:"
echo "CONDA_OVERRIDE_GLIBC=2.28 pixi run R -e \"install.packages('$NEW_NAME', repos=NULL, type='source'); library(SAIGEQTL)\""

echo ""
echo "=== NEXT STEPS ==="
echo "1. Test the package with different R versions using commands above"
echo "2. If tests pass, you're ready to push to GitHub"
echo "3. Create a release tag to trigger GitHub Actions for multi-platform binaries"