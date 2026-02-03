#!/usr/bin/env Rscript

#' SAIGEQTL Binary Installer
#' 
#' Simple script to download and install pre-compiled SAIGEQTL binaries
#' 
#' Usage:
#'   Rscript install_binary.R [custom_lib_path]
#'   
#' Examples:
#'   Rscript install_binary.R                      # Install to default location
#'   Rscript install_binary.R /path/to/R/library   # Install to custom location

cat("=== SAIGEQTL Binary Installer ===\n\n")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
custom_lib_path <- if (length(args) > 0) args[1] else NULL

if (!is.null(custom_lib_path)) {
  cat(sprintf("Custom library path: %s\n", custom_lib_path))
  if (!dir.exists(custom_lib_path)) {
    cat("Creating custom library directory...\n")
    dir.create(custom_lib_path, recursive = TRUE)
  }
}

# Install required packages if needed
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  cat("Installing jsonlite package...\n")
  install.packages("jsonlite", repos = "https://cloud.r-project.org/")
}

# Function to detect platform
detect_platform <- function() {
  if (Sys.info()["sysname"] == "Linux") {
    return("linux-x86_64")
  } else if (Sys.info()["sysname"] == "Darwin") {
    arch <- system("uname -m", intern = TRUE, ignore.stderr = TRUE)
    if (length(arch) > 0 && arch == "arm64") {
      return("darwin-arm64")
    } else {
      return("darwin-x86_64")
    }
  } else if (Sys.info()["sysname"] == "Windows") {
    return("windows-x86_64")
  } else {
    return("unknown")
  }
}

# Function to get R version
get_r_version <- function() {
  r_version <- R.Version()
  return(paste(r_version$major, strsplit(r_version$minor, "\\.")[[1]][1], sep = "."))
}

# Main installation function
install_saigeqtl_binary <- function() {
  platform <- detect_platform()
  r_version <- get_r_version()
  
  cat(sprintf("Platform: %s\n", platform))
  cat(sprintf("R Version: %s\n", r_version))
  
  if (platform == "unknown") {
    stop("Unsupported platform. Please use source installation.")
  }
  
  # File extension based on platform
  file_ext <- if (Sys.info()["sysname"] == "Windows") ".zip" else ".tgz"
  
  # Try different R versions in order of preference
  r_versions_to_try <- c(r_version, "4.3", "4.2", "4.1")
  if (platform == "darwin-arm64") {
    r_versions_to_try <- c(paste0(r_version, "-arm64"), r_versions_to_try)
  }
  
  success <- FALSE
  
  for (try_version in r_versions_to_try) {
    cat(sprintf("\nTrying R %s...\n", try_version))
    
    tryCatch({
      # Get latest release from GitHub API
      api_url <- "https://api.github.com/repos/weizhou0/qtl/releases/latest"
      cat("Fetching release information...\n")
      
      release_info <- jsonlite::fromJSON(api_url)
      
      # Find matching binary
      pattern <- sprintf("R-%s.*%s.*%s$", try_version, platform, gsub("\\.", "\\\\.", file_ext))
      matching_assets <- grep(pattern, release_info$assets$name, value = TRUE)
      
      if (length(matching_assets) > 0) {
        asset_name <- matching_assets[1]
        download_url <- release_info$assets$browser_download_url[
          release_info$assets$name == asset_name
        ]
        
        cat(sprintf("Found: %s\n", asset_name))
        cat("Downloading...\n")
        
        # Download binary
        temp_file <- tempfile(fileext = file_ext)
        download.file(download_url, temp_file, mode = "wb", quiet = FALSE)
        
        cat("Installing...\n")
        if (!is.null(custom_lib_path)) {
          install.packages(temp_file, repos = NULL, type = "binary", lib = custom_lib_path)
        } else {
          install.packages(temp_file, repos = NULL, type = "binary")
        }
        
        # Test installation
        cat("Testing installation...\n")
        if (!is.null(custom_lib_path)) {
          library(SAIGEQTL, lib.loc = custom_lib_path)
        } else {
          library(SAIGEQTL)
        }
        
        cat(sprintf("\n✓ Successfully installed SAIGEQTL binary (R %s)!\n", try_version))
        if (!is.null(custom_lib_path)) {
          cat(sprintf("\nIMPORTANT: SAIGEQTL installed to custom path: %s\n", custom_lib_path))
          cat("To use SAIGEQTL commands, specify the library path:\n")
          cat(sprintf("  step1_fitNULLGLMM_qtl.R --library='%s' [other_options]\n", custom_lib_path))
          cat(sprintf("  step2_tests_qtl.R --library='%s' [other_options]\n", custom_lib_path))
          cat(sprintf("  step3_gene_pvalue_qtl.R --library='%s' [other_options]\n", custom_lib_path))
          cat("\nIn R, load with:\n")
          cat(sprintf("  library(SAIGEQTL, lib.loc='%s')\n", custom_lib_path))
        }
        success <- TRUE
        break
        
      } else {
        cat(sprintf("No binary found for R %s on %s\n", try_version, platform))
      }
      
    }, error = function(e) {
      cat(sprintf("Failed to install R %s binary: %s\n", try_version, e$message))
    })
  }
  
  if (!success) {
    cat("\n✗ No compatible binary found.\n")
    cat("\nFallback options:\n")
    cat("1. Source installation: remotes::install_github('weizhou0/qtl')\n")
    cat("2. Smart installer: source('https://raw.githubusercontent.com/weizhou0/qtl/main/scripts/install.R')\n")
    cat("3. Docker: docker run -it weizhou0/saigeqtl\n")
    stop("Binary installation failed")
  }
}

# Run if called directly
if (!interactive()) {
  install_saigeqtl_binary()
}