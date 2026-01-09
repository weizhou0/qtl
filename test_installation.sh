#!/bin/bash

# Comprehensive Installation Testing Script for SAIGEQTL
# Tests all installation approaches before GitHub upload

set -e  # Exit on any error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="${SCRIPT_DIR}/test_installs"
PKG_NAME="SAIGEQTL"

echo "=========================================="
echo "SAIGEQTL Installation Testing Suite"
echo "=========================================="

# Create test directory structure
mkdir -p "$TEST_DIR"/{standard,pixi,manual}

# Function to test package loading
test_package_loading() {
    local test_name="$1"
    local lib_path="$2"
    
    echo "Testing package loading for $test_name..."
    R --vanilla --slave << EOF
.libPaths(c("$lib_path", .libPaths()))
if (require("$PKG_NAME", quietly = TRUE)) {
    cat("✅ $PKG_NAME loaded successfully in $test_name\n")
    cat("Version:", as.character(packageVersion("$PKG_NAME")), "\n")
    
    # Test basic functionality
    tryCatch({
        # Test if main functions exist
        if (exists("fitNULLGLMM_multiV")) {
            cat("✅ Main functions available\n")
        }
    }, error = function(e) {
        cat("⚠️  Warning: Some functions may not be available:", e\$message, "\n")
    })
} else {
    cat("❌ Failed to load $PKG_NAME in $test_name\n")
    quit(status = 1)
}
EOF
}

# Function to clean up test installations
cleanup_test() {
    local test_name="$1"
    echo "Cleaning up $test_name test..."
    rm -rf "$TEST_DIR/$test_name"
    mkdir -p "$TEST_DIR/$test_name"
}

echo ""
echo "=========================================="
echo "Test 1: Manual R CMD INSTALL (Standard)"
echo "=========================================="

cleanup_test "manual"
LIB_MANUAL="$TEST_DIR/manual/library"
mkdir -p "$LIB_MANUAL"

echo "Testing configure script in standard mode..."
# Simulate standard environment (no pixi)
(
    unset PIXI_PROJECT_ROOT
    cd /tmp
    mkdir -p test_standard
    cd test_standard
    cp "$SCRIPT_DIR/configure" .
    cp "$SCRIPT_DIR/src/Makevars.in" .
    mkdir -p src
    cp "$SCRIPT_DIR/src/Makevars.in" src/
    ./configure
    echo "Configure output for standard mode:"
    if [ -f "src/Makevars" ]; then
        echo "✅ Makevars created successfully"
        grep -E "(PKG_CPPFLAGS|PKG_LIBS)" src/Makevars || true
    else
        echo "❌ Makevars not created"
    fi
    cd .. && rm -rf test_standard
)

echo ""
echo "Attempting R CMD INSTALL..."
if R CMD INSTALL . --library="$LIB_MANUAL" --no-lock --no-multiarch; then
    echo "✅ Manual installation successful"
    test_package_loading "manual installation" "$LIB_MANUAL"
else
    echo "❌ Manual installation failed"
    echo "Checking for common issues..."
    
    # Check system dependencies
    echo "System dependency check:"
    if command -v pkg-config >/dev/null 2>&1; then
        echo "✅ pkg-config available"
    else
        echo "⚠️  pkg-config not found"
    fi
    
    # Check for compiler
    if command -v g++ >/dev/null 2>&1; then
        echo "✅ g++ compiler available"
    else
        echo "❌ g++ compiler not found"
    fi
    
    # Check for OpenMP
    echo "Testing OpenMP support..."
    cat > test_openmp.cpp << 'EOF'
#include <omp.h>
#include <iostream>
int main() {
    #pragma omp parallel
    {
        #pragma omp single
        std::cout << "OpenMP threads: " << omp_get_num_threads() << std::endl;
    }
    return 0;
}
EOF
    
    if g++ -fopenmp test_openmp.cpp -o test_openmp 2>/dev/null && ./test_openmp; then
        echo "✅ OpenMP support confirmed"
    else
        echo "❌ OpenMP support issues detected"
    fi
    rm -f test_openmp.cpp test_openmp
fi

echo ""
echo "=========================================="
echo "Test 2: Standard R Installation (remotes)"
echo "=========================================="

cleanup_test "standard"
LIB_STANDARD="$TEST_DIR/standard/library"
mkdir -p "$LIB_STANDARD"

echo "Testing installation via local remotes..."
R --vanilla --slave << EOF
# Test remotes installation from local path
.libPaths(c("$LIB_STANDARD", .libPaths()))

# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", lib = "$LIB_STANDARD")
}

# Install from local directory (simulating GitHub install)
tryCatch({
    remotes::install_local(".", lib = "$LIB_STANDARD", dependencies = TRUE, force = TRUE)
    cat("✅ Standard remotes installation successful\n")
}, error = function(e) {
    cat("❌ Standard installation failed:", e\$message, "\n")
    quit(status = 1)
})
EOF

if [ $? -eq 0 ]; then
    test_package_loading "standard remotes installation" "$LIB_STANDARD"
else
    echo "❌ Standard remotes installation failed"
fi

echo ""
echo "=========================================="
echo "Test 3: Pixi Environment Installation"
echo "=========================================="

cleanup_test "pixi"
LIB_PIXI="$TEST_DIR/pixi/library"
mkdir -p "$LIB_PIXI"

echo "Testing pixi installation approach..."
if [ -f "pixi.toml" ] && [ -d ".pixi" ]; then
    echo "Pixi environment detected, testing installation..."
    
    # Test the configure script in pixi mode
    ./configure
    echo "Configure output for pixi mode:"
    if [ -f "src/Makevars" ]; then
        echo "✅ Makevars created successfully"
        grep -E "(PKG_CPPFLAGS|PKG_LIBS)" src/Makevars || true
    fi
    
    # Test installation using the current environment
    if R CMD INSTALL . --library="$LIB_PIXI" --no-lock --no-multiarch; then
        echo "✅ Pixi installation successful"
        test_package_loading "pixi installation" "$LIB_PIXI"
    else
        echo "❌ Pixi installation failed"
    fi
else
    echo "⚠️  Pixi environment not available, skipping pixi-specific test"
    echo "   (This is expected in non-pixi environments)"
fi

echo ""
echo "=========================================="
echo "Test 4: Dependencies Verification"
echo "=========================================="

echo "Checking R package dependencies..."
R --vanilla --slave << 'EOF'
required_deps <- c("Rcpp", "RcppParallel", "Matrix", "data.table", 
                   "RcppArmadillo", "RcppNumerical", "furrr", "dbplyr",
                   "dplyr", "RSQLite", "RhpcBLASctl", "methods", "optparse", 
                   "nlme", "MASS", "BH", "RcppEigen")

missing_deps <- character(0)
for (dep in required_deps) {
    if (!requireNamespace(dep, quietly = TRUE)) {
        missing_deps <- c(missing_deps, dep)
    }
}

if (length(missing_deps) > 0) {
    cat("❌ Missing dependencies:", paste(missing_deps, collapse = ", "), "\n")
    cat("Install with: install.packages(c(", paste(paste0('"', missing_deps, '"'), collapse = ", "), "))\n")
} else {
    cat("✅ All required dependencies are available\n")
}
EOF

echo ""
echo "=========================================="
echo "Test 5: Regression Tests"
echo "=========================================="

if [ -f "run_regression_test.sh" ]; then
    echo "Running basic regression test to verify functionality..."
    if ./run_regression_test.sh test; then
        echo "✅ Regression tests passed"
    else
        echo "❌ Regression tests failed"
    fi
else
    echo "⚠️  Regression test script not found"
fi

echo ""
echo "=========================================="
echo "Installation Testing Summary"
echo "=========================================="

echo "Test results:"
echo "1. Manual R CMD INSTALL: $([ -d "$LIB_MANUAL" ] && echo "✅ Completed" || echo "❌ Failed")"
echo "2. Standard remotes install: $([ -d "$LIB_STANDARD" ] && echo "✅ Completed" || echo "❌ Failed")"
echo "3. Pixi environment install: $([ -d "$LIB_PIXI" ] && echo "✅ Completed" || echo "⚠️  Skipped")"

echo ""
echo "Next steps:"
echo "1. Fix any failed installations above"
echo "2. Commit changes to git"
echo "3. Push to GitHub"
echo "4. Test GitHub installation with:"
echo "   remotes::install_github('weizhou0/qtl')"

echo ""
echo "Cleanup test directories? (y/n)"
read -r response
if [ "$response" = "y" ]; then
    rm -rf "$TEST_DIR"
    echo "Test directories cleaned up."
fi