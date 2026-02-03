# SAIGEQTL Binary Installation

## Binary Information
- File: SAIGEQTL_0.3.4_R-4.4_linux-x86_64.tgz
- Built on: 2026-01-10 02:27:34.907702
- R Version: 4 4.3
- Platform: linux-x86_64

## Installation Instructions

### Method 1: Direct Installation
```r
install.packages("binaries/SAIGEQTL_0.3.4_R-4.4_linux-x86_64.tgz", repos = NULL, type = "binary")
```

### Method 2: Install to Custom Library
```r
# Create custom library directory
lib_dir <- "/path/to/custom/library"
dir.create(lib_dir, recursive = TRUE)

# Install binary
install.packages("binaries/SAIGEQTL_0.3.4_R-4.4_linux-x86_64.tgz", lib = lib_dir, repos = NULL, type = "binary")

# Load package
library(SAIGEQTL, lib.loc = lib_dir)
```

### Method 3: Test Installation
```r
# Test in temporary location
temp_lib <- tempfile("saigeqtl_test_")
dir.create(temp_lib)
install.packages("binaries/SAIGEQTL_0.3.4_R-4.4_linux-x86_64.tgz", lib = temp_lib, repos = NULL, type = "binary")
library(SAIGEQTL, lib.loc = temp_lib)

# Quick test
if (exists("fitNULLGLMM_multiV")) {
  cat("SAIGEQTL installed successfully!\n")
}
```

## Verification
After installation, verify with:
```r
library(SAIGEQTL)
packageVersion("SAIGEQTL")
```

