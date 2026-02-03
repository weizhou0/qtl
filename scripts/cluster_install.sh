#!/bin/bash

# SAIGEQTL Cluster Installation Script
# Designed for HPC environments that may lack proper compilers

set -e

echo "=== SAIGEQTL Cluster Installation ==="
echo "Designed for HPC systems without C++ compilers"
echo ""

# Function to install pixi
install_pixi() {
    if ! command -v pixi &> /dev/null; then
        echo "Installing pixi..."
        curl -fsSL https://pixi.sh/install.sh | bash
        export PATH="$HOME/.pixi/bin:$PATH"
        echo 'export PATH="$HOME/.pixi/bin:$PATH"' >> ~/.bashrc
    else
        echo "Pixi already installed"
    fi
}

# Function to setup pixi project
setup_pixi_project() {
    echo "Setting up SAIGEQTL pixi project..."
    
    # Create temporary directory
    TEMP_DIR=$(mktemp -d)
    cd "$TEMP_DIR"
    
    # Clone the repository
    git clone https://github.com/weizhou0/SAIGEQTL.git
    cd SAIGEQTL
    
    # Install dependencies using pixi
    pixi install
    
    # Install SAIGEQTL using pixi R environment
    pixi run install-github-deps
    pixi run R CMD INSTALL . --library="$HOME/R_libraries/SAIGEQTL"
    
    # Create activation script
    cat > "$HOME/activate_saigeqtl.sh" << 'EOF'
#!/bin/bash
# SAIGEQTL Environment Activation Script

export R_LIBS_USER="$HOME/R_libraries/SAIGEQTL"
export PATH="$HOME/.pixi/bin:$PATH"

echo "SAIGEQTL environment activated"
echo "Start R and run: library(SAIGEQTL)"
EOF
    
    chmod +x "$HOME/activate_saigeqtl.sh"
    
    echo ""
    echo "Installation complete!"
    echo "To use SAIGEQTL:"
    echo "  source ~/activate_saigeqtl.sh"
    echo "  R"
    echo "  > library(SAIGEQTL)"
    echo ""
    
    # Cleanup
    cd /
    rm -rf "$TEMP_DIR"
}

# Function for conda-based installation (alternative)
install_via_conda() {
    echo "Setting up conda environment for SAIGEQTL..."
    
    # Check if conda/mamba available
    if command -v mamba &> /dev/null; then
        CONDA_CMD="mamba"
    elif command -v conda &> /dev/null; then
        CONDA_CMD="conda"
    else
        echo "Neither conda nor mamba found. Installing miniconda..."
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh -b -p "$HOME/miniconda3"
        export PATH="$HOME/miniconda3/bin:$PATH"
        echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
        CONDA_CMD="conda"
    fi
    
    # Create environment
    $CONDA_CMD create -n saigeqtl -y r-base r-devtools r-rcpp r-rcpparmadillo r-matrix r-data.table gcc_linux-64 gxx_linux-64
    
    # Activate environment and install
    source activate saigeqtl
    R -e "devtools::install_github('weizhou0/qtl')"
    
    echo ""
    echo "Installation complete!"
    echo "To use SAIGEQTL:"
    echo "  conda activate saigeqtl"
    echo "  R"
    echo "  > library(SAIGEQTL)"
    echo ""
}

# Main installation logic
main() {
    echo "Detecting system capabilities..."
    
    # Check for existing R installation
    if ! command -v R &> /dev/null; then
        echo "R not found. Please install R first or use a module system:"
        echo "  module load R"
        exit 1
    fi
    
    # Check R version
    R_VERSION=$(R --version | head -n1 | grep -oE '[0-9]+\.[0-9]+' | head -n1)
    if [[ $(echo "$R_VERSION < 3.5" | bc -l) -eq 1 ]]; then
        echo "R version $R_VERSION found. Need R >= 3.5.0"
        exit 1
    fi
    echo "R version $R_VERSION ✓"
    
    # Check for compiler
    if command -v gcc &> /dev/null && command -v g++ &> /dev/null; then
        echo "Compilers found. You may be able to use standard installation:"
        echo "  R -e \"remotes::install_github('weizhou0/qtl')\""
        echo ""
        read -p "Continue with pixi installation anyway? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            exit 0
        fi
    fi
    
    # Choose installation method
    echo "Available installation methods:"
    echo "1) Pixi-based installation (recommended)"
    echo "2) Conda-based installation"
    echo "3) Exit and try manual installation"
    echo ""
    read -p "Choose method (1-3): " -n 1 -r method
    echo ""
    
    case $method in
        1)
            install_pixi
            setup_pixi_project
            ;;
        2)
            install_via_conda
            ;;
        3)
            echo "Manual installation instructions:"
            echo "1. Load required modules (if available):"
            echo "   module load gcc R"
            echo "2. Install in R:"
            echo "   install.packages('remotes')"
            echo "   remotes::install_github('weizhou0/qtl')"
            ;;
        *)
            echo "Invalid choice"
            exit 1
            ;;
    esac
}

# Run main function
main "$@"