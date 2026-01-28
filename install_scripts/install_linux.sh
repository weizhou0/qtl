#!/bin/bash
# SAIGEQTL Installation Script for Linux
#
# This script uses Miniforge (free, open-source) to install SAIGEQTL.
# Miniforge is required because the package depends on specialized C++
# libraries (savvy, superlu) only available via conda-forge.
#
# Miniforge is NOT the same as Anaconda - it's completely free and open-source.
# No sudo/root access required.

set -e  # Exit on error

echo "=============================================="
echo "SAIGEQTL Installation Script for Linux"
echo "=============================================="
echo ""
echo "This script uses Miniforge (free, open-source) to install"
echo "SAIGEQTL with all required dependencies. No sudo required."
echo ""

# Check for conda or mamba (provided by Miniforge)
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo "Found mamba - using for faster installation"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    echo "Found conda"
else
    echo "Miniforge not found. Installing now..."
    echo ""

    # Download and install Miniforge
    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
    rm Miniforge3-Linux-x86_64.sh

    # Initialize
    eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
    CONDA_CMD="mamba"
    echo "Miniforge installed successfully!"
fi

# Environment name
ENV_NAME="saigeqtl_env"

echo ""
echo "[1/4] Creating conda environment with system libraries..."
echo "----------------------------------------------"

# Check if environment already exists
if $CONDA_CMD env list | grep -q "^$ENV_NAME "; then
    echo "Environment '$ENV_NAME' already exists."
    read -p "Remove and recreate? [y/N] " response
    if [[ "$response" =~ ^[Yy]$ ]]; then
        $CONDA_CMD env remove -n $ENV_NAME -y
    else
        echo "Using existing environment..."
    fi
fi

# Create environment with all dependencies including savvy
$CONDA_CMD create -n $ENV_NAME -y \
    -c conda-forge -c bioconda \
    r-base=4.3 \
    r-devtools \
    r-rcpp \
    r-rcpparmadillo \
    r-rcppeigen \
    r-rcppparallel \
    r-rcppnumerical \
    r-bh \
    r-matrix \
    r-data.table \
    r-dplyr \
    r-dbplyr \
    r-rsqlite \
    r-optparse \
    r-skat \
    r-spatest \
    r-metaskat \
    r-qlcmatrix \
    r-rhpcblasctl \
    r-furrr \
    r-nlme \
    r-mass \
    boost-cpp \
    openblas \
    superlu \
    zlib \
    zstd \
    savvy \
    compilers

echo ""
echo "[2/4] Activating environment..."
echo "----------------------------------------------"

# Get conda base for activation
CONDA_BASE=$($CONDA_CMD info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
$CONDA_CMD activate $ENV_NAME

echo "Environment activated: $CONDA_PREFIX"

# Set up R Makevars for this environment
mkdir -p ~/.R
cat > ~/.R/Makevars << EOF
# Linux conda environment configuration for SAIGEQTL
CPPFLAGS += -I${CONDA_PREFIX}/include
LDFLAGS += -L${CONDA_PREFIX}/lib -Wl,-rpath,${CONDA_PREFIX}/lib
EOF

echo "Created ~/.R/Makevars pointing to conda environment"

echo ""
echo "[3/4] Installing additional R packages..."
echo "----------------------------------------------"

R --vanilla << 'EOF'
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

# Install GitHub packages
cat("Installing GitHub packages...\n")

if (!requireNamespace("fastSave", quietly = TRUE)) {
    remotes::install_github("barkasn/fastSave", upgrade = "never")
}

# Verify SPAtest version
if (packageVersion("SPAtest") != "3.1.2") {
    cat("Installing SPAtest version 3.1.2...\n")
    remotes::install_version("SPAtest", version = "3.1.2", upgrade = "never")
}

cat("GitHub packages installed!\n")
EOF

echo ""
echo "[4/4] Installing SAIGEQTL..."
echo "----------------------------------------------"

R --vanilla << 'EOF'
options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("Installing SAIGEQTL from GitHub...\n")
remotes::install_github("weizhou0/SAIGEQTL", dependencies = TRUE, upgrade = "never")

# Verify installation
cat("\n--- Verifying installation ---\n")
if (requireNamespace("SAIGEQTL", quietly = TRUE)) {
    library(SAIGEQTL)
    cat("SUCCESS: SAIGEQTL installed and loaded successfully!\n")
    cat(paste0("Version: ", packageVersion("SAIGEQTL"), "\n"))
} else {
    stop("FAILED: SAIGEQTL installation failed")
}
EOF

echo ""
echo "=============================================="
echo "Installation complete!"
echo "=============================================="
echo ""
echo "IMPORTANT: Always activate the environment before using SAIGEQTL:"
echo ""
echo "  conda activate $ENV_NAME"
echo "  R"
echo "  > library(SAIGEQTL)"
echo ""
echo "For HPC job scripts, add:"
echo "  source ~/.bashrc"
echo "  conda activate $ENV_NAME"
echo ""
