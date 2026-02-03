#!/bin/bash

# Build SAIGEQTL for multiple R versions
# This ensures compatibility across different R environments

set -e

echo "=== Multi-Version SAIGEQTL Builder ==="

# Available R versions (modify based on your system)
R_VERSIONS=("R/4.1" "R/4.2" "R/4.3" "R/4.4")

# Output directory
OUTPUT_DIR="multi_version_binaries"
mkdir -p "$OUTPUT_DIR"

# Function to build for specific R version
build_for_r_version() {
    local r_module="$1"
    echo ""
    echo "=== Building for $r_module ==="
    
    # Load specific R version
    if module load "$r_module" 2>/dev/null; then
        echo "✓ Loaded module: $r_module"
    else
        echo "✗ Failed to load module: $r_module"
        return 1
    fi
    
    # Get R version info
    local r_version=$(R --slave -e "cat(paste(R.Version()\$major, R.Version()\$minor, sep='.'))")
    local r_major_minor=$(R --slave -e "cat(paste(R.Version()\$major, strsplit(R.Version()\$minor, '\\\\.')[[1]][1], sep='.'))")
    echo "R Version: $r_version"
    
    # Create version-specific output directory
    local version_dir="$OUTPUT_DIR/R-$r_major_minor"
    mkdir -p "$version_dir"
    
    # Build package
    echo "Building package..."
    R --slave -e "
    source('build_binary.R')
    main(output_dir = '$version_dir', test_installation = TRUE)
    "
    
    if [ $? -eq 0 ]; then
        echo "✓ Build successful for $r_module"
        
        # List created packages
        echo "Created packages:"
        ls -la "$version_dir"/*.t*gz "$version_dir"/*.tar.gz 2>/dev/null || echo "No packages found"
        
    else
        echo "✗ Build failed for $r_module"
        return 1
    fi
    
    # Unload module for next iteration
    module unload "$r_module" 2>/dev/null || true
}

# Function to test cross-compatibility
test_cross_compatibility() {
    echo ""
    echo "=== Testing Cross-Compatibility ==="
    
    for r_module in "${R_VERSIONS[@]}"; do
        if ! module load "$r_module" 2>/dev/null; then
            echo "Skipping $r_module (not available)"
            continue
        fi
        
        local r_major_minor=$(R --slave -e "cat(paste(R.Version()\$major, strsplit(R.Version()\$minor, '\\\\.')[[1]][1], sep='.'))")
        echo ""
        echo "Testing with $r_module (R-$r_major_minor):"
        
        # Test packages from all R versions
        for pkg_dir in "$OUTPUT_DIR"/R-*; do
            if [ ! -d "$pkg_dir" ]; then continue; fi
            
            local built_version=$(basename "$pkg_dir")
            echo "  Testing $built_version packages..."
            
            # Find source packages
            local source_pkg=$(ls "$pkg_dir"/*source*.tar.gz 2>/dev/null | head -n1)
            if [ -n "$source_pkg" ]; then
                echo "    Testing: $(basename "$source_pkg")"
                R --slave -e "
                tryCatch({
                  temp_lib <- tempfile('test_')
                  dir.create(temp_lib)
                  install.packages('$source_pkg', lib = temp_lib, repos = NULL, type = 'source', quiet = TRUE)
                  library(SAIGEQTL, lib.loc = temp_lib, quietly = TRUE)
                  cat('      ✓ Compatible\\n')
                  unlink(temp_lib, recursive = TRUE)
                }, error = function(e) {
                  cat('      ✗ Incompatible:', e\$message, '\\n')
                })
                "
            fi
        done
        
        module unload "$r_module" 2>/dev/null || true
    done
}

# Main execution
echo "Available R modules:"
module avail R 2>&1 | grep -E "R/[0-9]" || echo "No R modules found"

echo ""
read -p "Continue with building? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled"
    exit 0
fi

# Build for each R version
successful_builds=0
for r_module in "${R_VERSIONS[@]}"; do
    if build_for_r_version "$r_module"; then
        ((successful_builds++))
    fi
done

echo ""
echo "=== BUILD SUMMARY ==="
echo "Successful builds: $successful_builds/${#R_VERSIONS[@]}"
echo "Output directory: $OUTPUT_DIR"

if [ $successful_builds -gt 0 ]; then
    echo ""
    echo "Testing cross-compatibility..."
    test_cross_compatibility
    
    echo ""
    echo "=== USAGE INSTRUCTIONS ==="
    echo "To install for a specific R version:"
    echo "  R-4.1: install.packages('$OUTPUT_DIR/R-4.1/*source*.tar.gz', repos=NULL, type='source')"
    echo "  R-4.2: install.packages('$OUTPUT_DIR/R-4.2/*source*.tar.gz', repos=NULL, type='source')"
    echo "  R-4.3: install.packages('$OUTPUT_DIR/R-4.3/*source*.tar.gz', repos=NULL, type='source')"
fi