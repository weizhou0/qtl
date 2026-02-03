#!/usr/bin/env Rscript

#' SAIGEQTL Smart Installation Script
#' 
#' This script provides multiple installation methods depending on the system:
#' 1. Standard remotes::install_github (if compiler available)
#' 2. Pixi-based installation (if pixi available)
#' 3. Pre-compiled binary download (fallback)
#' 4. Docker-based installation (last resort)

cat("=== SAIGEQTL Smart Installer ===\n\n")

# Function to check system capabilities
check_system <- function() {
  results <- list()
  
  # Check for R version
  r_version <- R.Version()
  results$r_version <- paste(r_version$major, r_version$minor, sep=".")
  results$r_ok <- as.numeric(results$r_version) >= 3.5
  
  # Check for compiler
  results$has_gcc <- system("which gcc", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  results$has_gpp <- system("which g++", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  results$has_clang <- system("which clang", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  results$compiler_ok <- results$has_gcc || results$has_gpp || results$has_clang
  
  # Check for pixi
  results$has_pixi <- system("which pixi", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  
  # Check for essential system libraries
  results$has_blas <- system("ldconfig -p | grep -i blas", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  results$has_lapack <- system("ldconfig -p | grep -i lapack", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  
  # Check for remotes package
  results$has_remotes <- requireNamespace("remotes", quietly = TRUE)
  
  # Check for jsonlite (needed for GitHub API)
  results$has_jsonlite <- requireNamespace("jsonlite", quietly = TRUE)
  
  return(results)
}

# Function to install via standard remotes
install_via_remotes <- function() {
  cat("Installing via remotes::install_github...\n")
  
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org/")
  }
  
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    install.packages("jsonlite", repos = "https://cloud.r-project.org/")
  }
  
  tryCatch({
    remotes::install_github("weizhou0/qtl", dependencies = TRUE, upgrade = "never")
    return(TRUE)
  }, error = function(e) {
    cat("Standard installation failed:", e$message, "\n")
    return(FALSE)
  })
}

# Function to install via pixi environment
install_via_pixi <- function() {
  cat("Installing via pixi environment...\n")
  
  # Check if we're already in a pixi environment
  in_pixi <- Sys.getenv("PIXI_PROJECT_ROOT") != "" || dir.exists(".pixi")
  
  if (!in_pixi) {
    cat("Setting up pixi environment...\n")
    system("curl -fsSL https://pixi.sh/install.sh | bash")
    Sys.setenv(PATH = paste(Sys.getenv("PATH"), "~/.pixi/bin", sep = ":"))
  }
  
  # Install using pixi
  result <- system("pixi run R -e \"remotes::install_github('weizhou0/qtl')\"")
  return(result == 0)
}

# Function to download pre-compiled binary
install_precompiled <- function() {
  cat("Attempting to download pre-compiled binary...\n")
  
  # Detect platform and R version
  r_version <- R.Version()
  r_major_minor <- paste(r_version$major, strsplit(r_version$minor, "\\.")[[1]][1], sep = ".")
  
  # Map R platform to our naming convention
  platform_map <- function() {
    if (Sys.info()["sysname"] == "Linux") {
      return("linux-x86_64")
    } else if (Sys.info()["sysname"] == "Darwin") {
      # Check for Apple Silicon
      arch <- system("uname -m", intern = TRUE)
      if (arch == "arm64") {
        return("darwin-arm64")
      } else {
        return("darwin-x86_64")
      }
    } else if (Sys.info()["sysname"] == "Windows") {
      return("windows-x86_64")
    } else {
      return(R.Version()$platform)
    }
  }
  
  platform <- platform_map()
  
  # Determine file extension
  file_ext <- if (Sys.info()["sysname"] == "Windows") ".zip" else ".tgz"
  
  # Try different R version suffixes
  r_suffixes <- c(
    r_major_minor,
    paste0(r_major_minor, "-arm64"),  # For Apple Silicon
    "4.3", "4.2", "4.1"  # Fallback versions
  )
  
  for (r_suffix in r_suffixes) {
    binary_name <- sprintf("SAIGEQTL_*_R-%s_%s%s", r_suffix, platform, file_ext)
    
    # Get latest release info from GitHub API
    tryCatch({
      cat(sprintf("  Trying R %s for %s...\n", r_suffix, platform))
      
      # Get release assets
      api_url <- "https://api.github.com/repos/weizhou0/qtl/releases/latest"
      release_info <- jsonlite::fromJSON(api_url)
      
      # Find matching asset
      matching_assets <- grep(sprintf("R-%s.*%s", r_suffix, platform), 
                             release_info$assets$name, value = TRUE)
      
      if (length(matching_assets) > 0) {
        asset_name <- matching_assets[1]
        download_url <- release_info$assets$browser_download_url[
          release_info$assets$name == asset_name
        ]
        
        cat(sprintf("  Found binary: %s\n", asset_name))
        
        # Download and install
        temp_file <- tempfile(fileext = file_ext)
        download.file(download_url, temp_file, mode = "wb", quiet = TRUE)
        
        install.packages(temp_file, repos = NULL, type = "binary")
        cat(sprintf("  ✓ Successfully installed binary for R %s\n", r_suffix))
        return(TRUE)
      }
    }, error = function(e) {
      # Continue to next version
    })
  }
  
  cat("  ✗ No compatible pre-compiled binary found\n")
  return(FALSE)
}

# Function to provide Docker installation
suggest_docker <- function() {
  cat("\n=== Docker Installation Option ===\n")
  cat("If other methods fail, you can use our Docker image:\n\n")
  cat("docker pull weizhou0/saigeqtl\n")
  cat("docker run -it --rm -v $(pwd):/data weizhou0/saigeqtl R\n\n")
  cat("Then in R:\n")
  cat("library(SAIGEQTL)\n\n")
}

# Main installation logic
main <- function() {
  cat("Checking system configuration...\n")
  sys_check <- check_system()
  
  cat(sprintf("  R version: %s %s\n", sys_check$r_version, 
              ifelse(sys_check$r_ok, "✓", "✗ (need >= 3.5)")))
  cat(sprintf("  Compiler: %s\n", ifelse(sys_check$compiler_ok, "✓", "✗")))
  cat(sprintf("  Pixi: %s\n", ifelse(sys_check$has_pixi, "✓", "✗")))
  cat(sprintf("  BLAS/LAPACK: %s\n", ifelse(sys_check$has_blas && sys_check$has_lapack, "✓", "✗")))
  
  if (!sys_check$r_ok) {
    stop("R version >= 3.5.0 required")
  }
  
  # Strategy 1: Standard installation if compiler available
  if (sys_check$compiler_ok && sys_check$has_blas && sys_check$has_lapack) {
    cat("\nAttempting standard installation...\n")
    if (install_via_remotes()) {
      cat("\n✓ Installation successful via remotes::install_github!\n")
      return(invisible(TRUE))
    }
  }
  
  # Strategy 2: Pixi installation if available
  if (sys_check$has_pixi) {
    cat("\nAttempting pixi-based installation...\n")
    if (install_via_pixi()) {
      cat("\n✓ Installation successful via pixi!\n")
      return(invisible(TRUE))
    }
  }
  
  # Strategy 3: Pre-compiled binary
  cat("\nAttempting pre-compiled binary installation...\n")
  if (install_precompiled()) {
    cat("\n✓ Installation successful via pre-compiled binary!\n")
    return(invisible(TRUE))
  }
  
  # Strategy 4: Suggest Docker
  cat("\n✗ All installation methods failed\n")
  suggest_docker()
  
  cat("\n=== Manual Installation Instructions ===\n")
  cat("1. Install system dependencies:\n")
  cat("   # Ubuntu/Debian:\n")
  cat("   sudo apt-get install build-essential libopenblas-dev liblapack-dev\n\n")
  cat("   # CentOS/RHEL:\n")
  cat("   sudo yum groupinstall 'Development Tools'\n")
  cat("   sudo yum install openblas-devel lapack-devel\n\n")
  cat("   # macOS:\n")
  cat("   xcode-select --install\n\n")
  cat("2. Then try again:\n")
  cat("   remotes::install_github('weizhou0/qtl')\n\n")
  
  return(invisible(FALSE))
}

# Run if called directly
if (!interactive()) {
  main()
}