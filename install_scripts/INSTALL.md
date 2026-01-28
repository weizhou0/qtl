# SAIGEQTL Installation Guide

## Prerequisites

**SAIGEQTL uses Miniforge** (free, open-source) because it depends on specialized C++ libraries (`savvy` for VCF/BCF file handling) that are only available via conda-forge.

**Important:** Miniforge is NOT Anaconda - it's completely free, open-source, and works on HPC clusters. No sudo/root access required.

### Install Miniforge (automatic or manual)

The installation scripts will **automatically install Miniforge** if not found. You can also install it manually:

**Linux:**
```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b
source ~/miniforge3/bin/activate
```

**macOS (Intel):**
```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh
bash Miniforge3-MacOSX-x86_64.sh -b
source ~/miniforge3/bin/activate
```

**macOS (Apple Silicon/M1/M2/M3):**
```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3-MacOSX-arm64.sh -b
source ~/miniforge3/bin/activate
```

---

## Quick Installation

### One-Command Install

**Linux or macOS:**
```bash
# Download and run the installation script
curl -fsSL https://raw.githubusercontent.com/weizhou0/SAIGEQTL/main/install_scripts/install_linux.sh | bash

# Or for macOS specifically:
curl -fsSL https://raw.githubusercontent.com/weizhou0/SAIGEQTL/main/install_scripts/install_macOS.sh | bash
```

### Or Clone and Run Locally

```bash
git clone https://github.com/weizhou0/SAIGEQTL.git
cd SAIGEQTL/install_scripts

# Linux
chmod +x install_linux.sh && ./install_linux.sh

# macOS
chmod +x install_macOS.sh && ./install_macOS.sh
```

---

## Using SAIGEQTL

After installation, **always activate the environment first**:

```bash
conda activate saigeqtl_env
R
```

```r
library(SAIGEQTL)
```

---

## Manual Installation (Advanced)

If you prefer to set up the environment manually:

```bash
# Create conda environment with all dependencies
conda create -n saigeqtl_env -y \
    -c conda-forge -c bioconda \
    r-base=4.3 \
    r-devtools r-rcpp r-rcpparmadillo r-rcppeigen r-rcppparallel \
    r-rcppnumerical r-bh r-matrix r-data.table r-dplyr r-dbplyr \
    r-rsqlite r-optparse r-skat r-spatest r-metaskat r-qlcmatrix \
    r-rhpcblasctl r-furrr r-nlme r-mass \
    boost-cpp openblas superlu zlib zstd savvy compilers

# Activate and install SAIGEQTL
conda activate saigeqtl_env

# Configure R to find conda libraries
mkdir -p ~/.R
echo "CPPFLAGS += -I\${CONDA_PREFIX}/include" >> ~/.R/Makevars
echo "LDFLAGS += -L\${CONDA_PREFIX}/lib" >> ~/.R/Makevars

# Install SAIGEQTL from R
R -e 'remotes::install_github("weizhou0/SAIGEQTL")'
```

---

## System Requirements

| Component | Requirement |
|-----------|-------------|
| **OS** | Linux (x86_64) or macOS (Intel/Apple Silicon) |
| **RAM** | 8GB minimum, 32GB+ recommended |
| **Disk** | ~5GB for conda environment |
| **Conda** | Miniforge, Miniconda, or Anaconda |

### Key Dependencies (handled by conda)

| Library | Purpose |
|---------|---------|
| R ≥ 4.1 | Statistical computing |
| savvy | VCF/BCF file handling |
| Boost | C++ utilities |
| OpenBLAS | Linear algebra |
| SuperLU | Sparse linear solver |
| zstd | Compression |
| OpenMP | Parallelization |

---

## Troubleshooting

### "conda: command not found"
Initialize conda in your shell:
```bash
~/miniforge3/bin/conda init bash  # or zsh
source ~/.bashrc  # or ~/.zshrc
```

### "environment saigeqtl_env already exists"
Remove and recreate:
```bash
conda env remove -n saigeqtl_env
# Then run the install script again
```

### Compilation errors on macOS
Ensure Xcode command line tools are installed:
```bash
xcode-select --install
```

### Package loading errors
Make sure the environment is activated:
```bash
conda activate saigeqtl_env
R -e "library(SAIGEQTL)"
```

---

## Getting Help

- Documentation: https://weizhou0.github.io/SAIGE-QTL-doc/
- Issues: https://github.com/weizhou0/SAIGEQTL/issues
