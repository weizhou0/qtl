# Using Custom Library Paths with SAIGEQTL

This guide explains how to install SAIGEQTL to a custom location and use it with the `--library` option in analysis scripts.

## Why Use Custom Library Paths?

- **Shared systems**: Install in personal directory when you don't have admin access
- **HPC clusters**: Install in user space or project-specific directories
- **Multiple versions**: Keep different versions in separate locations
- **Isolation**: Separate package environments for different projects

## Installing to Custom Location

### Method 1: Binary Installation (Recommended)

```bash
# Install to custom directory
Rscript scripts/install_binary.R /path/to/your/R/library

# Example: Install to home directory
Rscript scripts/install_binary.R ~/R-packages/

# Example: Install to project directory  
Rscript scripts/install_binary.R /project/myproject/R-packages/
```

### Method 2: Source Installation

```r
# Install to custom location
install.packages("remotes")
remotes::install_github("weizhou0/qtl", lib = "/path/to/your/R/library")
```

### Method 3: Pixi Installation

```bash
# Set custom library path environment variable
export SAIGEQTL_LIB_PATH="/path/to/your/R/library"
mkdir -p "$SAIGEQTL_LIB_PATH"
pixi run install-standard
```

## Using SAIGEQTL with Custom Library Path

After installing to a custom location, you need to specify the library path when using SAIGEQTL.

### Command Line Tools

**Step 1: Fit NULL GLMM**
```bash
step1_fitNULLGLMM_qtl.R \
  --library='/path/to/your/R/library' \
  --plinkFile=input/genotype \
  --phenoFile=input/phenotype.txt \
  --outputPrefix=output/step1_results
```

**Step 2: Association Testing**
```bash
step2_tests_qtl.R \
  --library='/path/to/your/R/library' \
  --vcfFile=input/genotype.vcf.gz \
  --nullModelFile=output/step1_results.rda \
  --outputFile=output/step2_results.txt
```

**Step 3: Gene-level P-values**
```bash
step3_gene_pvalue_qtl.R \
  --library='/path/to/your/R/library' \
  --inputFile=output/step2_results.txt \
  --outputFile=output/step3_gene_pvalues.txt
```

### R Session

**Load SAIGEQTL from custom location:**
```r
# Load package from custom path
library(SAIGEQTL, lib.loc = "/path/to/your/R/library")

# Or add to library paths permanently
.libPaths("/path/to/your/R/library")
library(SAIGEQTL)

# Verify installation
packageVersion("SAIGEQTL")
```

### Make Custom Path Permanent

**Option 1: Set R environment variable**
```bash
# Add to ~/.bashrc or ~/.profile
export R_LIBS_USER="/path/to/your/R/library:$R_LIBS_USER"
```

**Option 2: Create .Renviron file**
```bash
# Add to ~/.Renviron
echo 'R_LIBS_USER="/path/to/your/R/library"' >> ~/.Renviron
```

**Option 3: Modify R profile**
```r
# Add to ~/.Rprofile
.libPaths(c("/path/to/your/R/library", .libPaths()))
```

## Examples for Common Scenarios

### HPC Cluster Installation

```bash
# 1. Install to user directory
mkdir -p $HOME/R-packages
Rscript scripts/install_binary.R $HOME/R-packages

# 2. Add to job scripts
echo 'export R_LIBS_USER="$HOME/R-packages:$R_LIBS_USER"' >> ~/.bashrc
source ~/.bashrc

# 3. Use in analysis (no --library needed after step 2)
step1_fitNULLGLMM_qtl.R --plinkFile=data/geno --phenoFile=data/pheno.txt
```

### Project-Specific Installation

```bash
# 1. Create project R library
PROJECT_DIR="/projects/eqtl_study"
R_LIB="$PROJECT_DIR/R-packages"
mkdir -p "$R_LIB"

# 2. Install SAIGEQTL
Rscript scripts/install_binary.R "$R_LIB"

# 3. Create wrapper script for convenience
cat > "$PROJECT_DIR/run_saigeqtl.sh" << 'EOF'
#!/bin/bash
SCRIPT_NAME=$(basename "$1")
shift  # Remove first argument (script name)
exec $SCRIPT_NAME --library="/projects/eqtl_study/R-packages" "$@"
EOF
chmod +x "$PROJECT_DIR/run_saigeqtl.sh"

# 4. Use wrapper
./run_saigeqtl.sh step1_fitNULLGLMM_qtl.R --plinkFile=data/geno
```

### Multiple Versions

```bash
# Install different versions to separate directories
Rscript scripts/install_binary.R ~/R-packages/saigeqtl-v0.3.4
# Later version:
# Rscript scripts/install_binary.R ~/R-packages/saigeqtl-v0.3.5

# Use specific version
step1_fitNULLGLMM_qtl.R \
  --library='~/R-packages/saigeqtl-v0.3.4' \
  --plinkFile=data/geno
```

## Troubleshooting

### Common Issues

**1. Permission denied when creating directory**
```bash
# Solution: Use directory you have write access to
Rscript scripts/install_binary.R ~/my-R-packages/
```

**2. Package not found after installation**
```bash
# Check if package was installed correctly
R -e "library(SAIGEQTL, lib.loc='/path/to/your/R/library'); packageVersion('SAIGEQTL')"
```

**3. Command line tools can't find package**
```bash
# Make sure to use --library option
step1_fitNULLGLMM_qtl.R --library='/path/to/your/R/library' --help
```

**4. Forgot where you installed it**
```bash
# Find SAIGEQTL installations
find ~ -name "SAIGEQTL" -type d 2>/dev/null
```

### Verification

**Check installation:**
```bash
# Test with custom library path
R -e "library(SAIGEQTL, lib.loc='/path/to/your/R/library'); cat('✓ SAIGEQTL version:', as.character(packageVersion('SAIGEQTL')), '\n')"
```

**Test command line tools:**
```bash
step1_fitNULLGLMM_qtl.R --library='/path/to/your/R/library' --help
```

## Best Practices

1. **Document your library path**: Save the installation path for future reference
2. **Use absolute paths**: Avoid relative paths that might change
3. **Set environment variables**: For frequently used custom paths
4. **Test installation**: Always verify with `--help` or simple R commands
5. **Keep track of versions**: Use versioned directories for multiple installations

## Environment Variables Summary

| Variable | Purpose | Example |
|----------|---------|---------|
| `R_LIBS_USER` | Additional library paths for R | `export R_LIBS_USER="/custom/path:$R_LIBS_USER"` |
| `SAIGEQTL_LIB_PATH` | Pixi installation target | `export SAIGEQTL_LIB_PATH="/custom/path"` |
| `R_LIBS` | System-wide library paths | `export R_LIBS="/custom/path"` |

Choose the method that best fits your environment and workflow!