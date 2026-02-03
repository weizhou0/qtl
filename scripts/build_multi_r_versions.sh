#!/bin/bash

# Build SAIGEQTL binaries for multiple R versions using pixi

set -e

echo "=== Multi-R Version Binary Builder (Pixi) ==="

# R versions to try
R_VERSIONS=("4.1" "4.2" "4.3" "4.4")

# Create output directory
OUTPUT_DIR="multi_r_binaries"
mkdir -p "$OUTPUT_DIR"

# Function to test if R version is available
test_r_version() {
    local r_version="$1"
    echo "Testing R $r_version availability..."
    
    # Create temporary pixi.toml for testing
    cat > pixi-test.toml << EOF
[workspace]
channels = ["conda-forge"]
name = "test"
platforms = ["linux-64"]

[dependencies]
r-base = "~${r_version}.0"
EOF

    # Test if this R version can be installed
    if timeout 30s pixi install --manifest-path pixi-test.toml > /dev/null 2>&1; then
        echo "✓ R $r_version is available"
        rm -f pixi-test.toml
        rm -rf .pixi-test 2>/dev/null || true
        return 0
    else
        echo "✗ R $r_version not available"
        rm -f pixi-test.toml
        rm -rf .pixi-test 2>/dev/null || true
        return 1
    fi
}

# Function to build for specific R version
build_for_r_version() {
    local r_version="$1"
    echo ""
    echo "=== Building SAIGEQTL for R $r_version ==="
    
    # Create version-specific pixi.toml
    local pixi_file="pixi-r${r_version//.}.toml"
    
    cat > "$pixi_file" << EOF
[workspace]
channels = ["conda-forge", "bioconda"]
name = "SAIGEQTL-R${r_version//.}"
platforms = ["linux-64"]

[dependencies]
# System dependencies
boost-cpp = "*"
cget = "*"
openblas = "*"
zlib = "*"
zstd = "*"
superlu = "*"
zlib-ng = "*"
nomkl = "*"
llvm-openmp = "*"
cxx-compiler = "*"
make = "*"
libblas = { version = "*", build = "*openblas" }
liblapack = { version = "*", build = "*openblas" }

# R and packages
r-base = "~${r_version}.0"
r-devtools = "*"
r-rcpp = "*"
r-rcpparmadillo = "*"
r-rcppeigen = "*"
r-rcppparallel = "*"
r-rcppnumerical = "*"
r-bh = "*"
"r-data.table" = "*"
r-dplyr = "*"
r-matrix = "*"
r-nlme = "*"
r-mass = "*"
r-optparse = "*"
r-rsqlite = "*"
r-remotes = "*"
htslib = ">=1.22.1,<2"

[tasks]
install-deps = "Rscript -e 'if(!require(remotes)) install.packages(\"remotes\"); remotes::install_github(c(\"cysouw/qlcMatrix\", \"leeshawn/MetaSKAT\", \"barkasn/fastSave\"))'"
build = "Rscript -e 'if(!require(pkgbuild)) install.packages(\"pkgbuild\"); pkg <- pkgbuild::build(binary=FALSE, dest_path=\"${OUTPUT_DIR}\"); cat(\"Built:\", pkg, \"\n\")'"
EOF

    # Install environment
    echo "Installing R $r_version environment..."
    if ! CONDA_OVERRIDE_GLIBC=2.28 pixi install --manifest-path "$pixi_file"; then
        echo "✗ Failed to install R $r_version environment"
        rm -f "$pixi_file"
        return 1
    fi
    
    # Install GitHub dependencies
    echo "Installing GitHub dependencies..."
    if ! CONDA_OVERRIDE_GLIBC=2.28 pixi run --manifest-path "$pixi_file" install-deps; then
        echo "! GitHub dependencies failed (non-critical)"
    fi
    
    # Build package
    echo "Building package..."
    if CONDA_OVERRIDE_GLIBC=2.28 pixi run --manifest-path "$pixi_file" build; then
        echo "✓ Build successful for R $r_version"
        
        # Rename package with R version info
        latest_pkg=$(ls -t "${OUTPUT_DIR}"/SAIGEQTL_*.tar.gz 2>/dev/null | head -1)
        if [ -n "$latest_pkg" ]; then
            new_name="${OUTPUT_DIR}/SAIGEQTL_0.3.4_R-${r_version}_linux-x86_64_source.tar.gz"
            mv "$latest_pkg" "$new_name"
            echo "Package renamed to: $(basename "$new_name")"
            
            # Quick test
            echo "Testing package..."
            CONDA_OVERRIDE_GLIBC=2.28 pixi run --manifest-path "$pixi_file" \
                Rscript -e "
                temp_lib <- tempfile('test_')
                dir.create(temp_lib)
                install.packages('$new_name', lib=temp_lib, repos=NULL, type='source')
                library(SAIGEQTL, lib.loc=temp_lib)
                cat('✓ R', '$r_version', 'package test successful\\n')
                "
        fi
    else
        echo "✗ Build failed for R $r_version"
        rm -f "$pixi_file"
        return 1
    fi
    
    # Cleanup
    rm -f "$pixi_file"
    return 0
}

# Main execution
echo "Checking available R versions in conda-forge..."

# Test which R versions are available
available_versions=()
for version in "${R_VERSIONS[@]}"; do
    if test_r_version "$version"; then
        available_versions+=("$version")
    fi
done

if [ ${#available_versions[@]} -eq 0 ]; then
    echo "No R versions available. Exiting."
    exit 1
fi

echo ""
echo "Available R versions: ${available_versions[*]}"
echo ""

# Build for each available version
successful_builds=0
for r_version in "${available_versions[@]}"; do
    if build_for_r_version "$r_version"; then
        ((successful_builds++))
    fi
    echo ""
done

# Summary
echo "=== BUILD SUMMARY ==="
echo "Successful builds: $successful_builds/${#available_versions[@]}"
echo "Available versions: ${available_versions[*]}"
echo "Output directory: $OUTPUT_DIR"

if [ $successful_builds -gt 0 ]; then
    echo ""
    echo "Created packages:"
    ls -la "$OUTPUT_DIR"/*.tar.gz 2>/dev/null || echo "No packages found"
    
    echo ""
    echo "Test with:"
    echo "  install.packages('$OUTPUT_DIR/SAIGEQTL_*_R-X.X_*.tar.gz', repos=NULL, type='source')"
fi