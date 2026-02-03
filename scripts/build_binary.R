#!/usr/bin/env Rscript

#' Manual Binary Building Script for SAIGEQTL
#' 
#' Use this to test binary creation before pushing to GitHub

cat("=== SAIGEQTL Manual Binary Builder ===\n")
cat("Testing binary creation locally\n\n")

# Function to get package info
get_package_info <- function() {
  if (!file.exists("DESCRIPTION")) {
    stop("DESCRIPTION file not found. Run this from package root directory.")
  }
  
  desc <- read.dcf("DESCRIPTION")
  list(
    name = as.character(desc[,"Package"]),
    version = as.character(desc[,"Version"])
  )
}

# Function to detect current platform
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

# Function to setup build environment
setup_build_environment <- function() {
  cat("Setting up build environment...\n")
  
  # Install required packages
  required_packages <- c("devtools", "pkgbuild", "remotes")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("  Installing %s...\n", pkg))
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
  }
  
  # Run configure script if it exists
  if (file.exists("configure")) {
    cat("  Running configure script...\n")
    if (Sys.info()["sysname"] == "Windows") {
      system("configure.win")
    } else {
      system("bash configure")
    }
  }
  
  cat("  Build environment ready!\n")
}

# Function to install dependencies
install_dependencies <- function() {
  cat("Installing package dependencies...\n")
  
  tryCatch({
    remotes::install_deps(dependencies = TRUE, upgrade = "never")
    cat("  Dependencies installed successfully!\n")
  }, error = function(e) {
    cat("  Warning: Some dependencies failed to install:\n")
    cat(sprintf("  %s\n", e$message))
    cat("  Continuing with available dependencies...\n")
  })
}

# Function to build binary
build_binary <- function(output_dir = ".") {
  cat("Building binary package...\n")
  
  # Get package info
  pkg_info <- get_package_info()
  r_version <- R.Version()
  r_major_minor <- paste(r_version$major, strsplit(r_version$minor, "\\.")[[1]][1], sep = ".")
  platform <- detect_platform()
  
  cat(sprintf("  Package: %s v%s\n", pkg_info$name, pkg_info$version))
  cat(sprintf("  R Version: %s\n", r_major_minor))
  cat(sprintf("  Platform: %s\n", platform))
  
  # Build both source and binary
  source_path <- NULL
  binary_path <- NULL
  
  tryCatch({
    # Build source package first
    cat("  Building source package...\n")
    source_path <- pkgbuild::build(binary = FALSE, dest_path = output_dir)
    cat(sprintf("  Source built: %s\n", source_path))
    
    # Build binary package (compiled)
    cat("  Building binary package...\n")
    binary_path <- pkgbuild::build(binary = TRUE, dest_path = output_dir)
    cat(sprintf("  Binary built: %s\n", binary_path))
    
    # For Linux, check if we actually got a compiled binary
    if (Sys.info()["sysname"] == "Linux") {
      # Check if binary contains compiled code
      temp_dir <- tempfile()
      dir.create(temp_dir)
      system(sprintf("cd %s && tar -xzf %s", temp_dir, shQuote(binary_path)))
      
      # Look for compiled shared objects
      so_files <- list.files(temp_dir, pattern = "\\.so$", recursive = TRUE)
      if (length(so_files) > 0) {
        cat("  ✓ Binary contains compiled shared libraries\n")
        
        # Rename binary with platform info
        file_ext <- ".tgz"
        new_name <- sprintf("%s_%s_R-%s_%s_binary%s", 
                            pkg_info$name, pkg_info$version, 
                            r_major_minor, platform, file_ext)
        
        new_binary_path <- file.path(output_dir, new_name)
        file.rename(binary_path, new_binary_path)
        cat(sprintf("  Renamed binary to: %s\n", new_name))
        binary_path <- new_binary_path
      } else {
        cat("  ! Binary is same as source (no compiled code)\n")
        # Remove duplicate and keep only source
        file.remove(binary_path)
        binary_path <- NULL
      }
      
      unlink(temp_dir, recursive = TRUE)
    }
    
    # Rename source with platform info
    source_ext <- ".tar.gz"
    new_source_name <- sprintf("%s_%s_R-%s_%s_source%s", 
                        pkg_info$name, pkg_info$version, 
                        r_major_minor, platform, source_ext)
    
    new_source_path <- file.path(output_dir, new_source_name)
    file.rename(source_path, new_source_path)
    cat(sprintf("  Renamed source to: %s\n", new_source_name))
    
    return(list(source = new_source_path, binary = binary_path))
    
  }, error = function(e) {
    cat(sprintf("  Build failed: %s\n", e$message))
    return(NULL)
  })
}

# Function to test binary
test_binary <- function(package_paths, test_lib_dir = NULL) {
  cat("Testing package installation...\n")
  
  if (is.null(package_paths)) {
    cat("  No packages to test\n")
    return(FALSE)
  }
  
  # Handle both single path and list of paths
  if (is.character(package_paths)) {
    # Single package path
    test_paths <- package_paths
  } else {
    # List with source and binary
    test_paths <- c(package_paths$source, package_paths$binary)
    test_paths <- test_paths[!is.null(test_paths)]
  }
  
  if (length(test_paths) == 0) {
    cat("  No valid packages to test\n")
    return(FALSE)
  }
  
  # Create temporary library if not specified
  if (is.null(test_lib_dir)) {
    test_lib_dir <- tempfile("test_lib_")
    dir.create(test_lib_dir, recursive = TRUE)
    on.exit(unlink(test_lib_dir, recursive = TRUE), add = TRUE)
  } else {
    dir.create(test_lib_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  success <- FALSE
  
  for (pkg_path in test_paths) {
    if (!file.exists(pkg_path)) next
    
    cat(sprintf("  Testing: %s\n", basename(pkg_path)))
    
    tryCatch({
      # Determine installation type based on filename
      pkg_type <- if (grepl("_binary\\.", pkg_path)) {
        "source"  # Even "binary" packages on Linux are installed as source
      } else {
        "source"
      }
      
      # Install package to test location
      cat(sprintf("    Installing to: %s\n", test_lib_dir))
      install.packages(pkg_path, lib = test_lib_dir, repos = NULL, type = pkg_type)
      
      # Test loading
      cat("    Testing package loading...\n")
      library(SAIGEQTL, lib.loc = test_lib_dir)
      
      # Test basic functionality
      cat("    Testing basic functions...\n")
      if (exists("fitNULLGLMM_multiV")) {
        cat("      ✓ Main function fitNULLGLMM_multiV available\n")
      }
      if (exists("SAIGE_SPATest")) {
        cat("      ✓ Test function SAIGE_SPATest available\n")
      }
      
      # Check if package has compiled code
      pkg_lib_path <- file.path(test_lib_dir, "SAIGEQTL")
      so_files <- list.files(pkg_lib_path, pattern = "\\.so$", recursive = TRUE)
      if (length(so_files) > 0) {
        cat("      ✓ Package contains compiled shared libraries\n")
      }
      
      # Unload for clean test
      try(detach("package:SAIGEQTL", unload = TRUE), silent = TRUE)
      
      cat("    ✓ Package test successful!\n")
      success <- TRUE
      break  # Success, no need to test other packages
      
    }, error = function(e) {
      cat(sprintf("    ✗ Package test failed: %s\n", e$message))
    })
  }
  
  return(success)
}

# Function to create installation instructions
create_install_instructions <- function(binary_path, output_dir) {
  if (is.null(binary_path)) return()
  
  instructions_file <- file.path(output_dir, "BINARY_INSTALL.md")
  
  instructions <- sprintf("# SAIGEQTL Binary Installation

## Binary Information
- File: %s
- Built on: %s
- R Version: %s
- Platform: %s

## Installation Instructions

### Method 1: Direct Installation
```r
install.packages(\"%s\", repos = NULL, type = \"binary\")
```

### Method 2: Install to Custom Library
```r
# Create custom library directory
lib_dir <- \"/path/to/custom/library\"
dir.create(lib_dir, recursive = TRUE)

# Install binary
install.packages(\"%s\", lib = lib_dir, repos = NULL, type = \"binary\")

# Load package
library(SAIGEQTL, lib.loc = lib_dir)
```

### Method 3: Test Installation
```r
# Test in temporary location
temp_lib <- tempfile(\"saigeqtl_test_\")
dir.create(temp_lib)
install.packages(\"%s\", lib = temp_lib, repos = NULL, type = \"binary\")
library(SAIGEQTL, lib.loc = temp_lib)

# Quick test
if (exists(\"fitNULLGLMM_multiV\")) {
  cat(\"SAIGEQTL installed successfully!\\n\")
}
```

## Verification
After installation, verify with:
```r
library(SAIGEQTL)
packageVersion(\"SAIGEQTL\")
```
",
    basename(binary_path),
    Sys.time(),
    paste(R.Version()$major, R.Version()$minor),
    detect_platform(),
    binary_path,
    binary_path,
    binary_path
  )
  
  writeLines(instructions, instructions_file)
  cat(sprintf("  Installation instructions written to: %s\n", instructions_file))
}

# Main function
main <- function(output_dir = "binaries", test_installation = TRUE) {
  cat("Starting binary build process...\n\n")
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Step 1: Setup
  setup_build_environment()
  cat("\n")
  
  # Step 2: Install dependencies
  install_dependencies()
  cat("\n")
  
  # Step 3: Build packages
  build_result <- build_binary(output_dir)
  cat("\n")
  
  # Step 4: Test packages (optional)
  if (test_installation && !is.null(build_result)) {
    test_success <- test_binary(build_result)
    cat("\n")
    
    if (!test_success) {
      cat("⚠️  Package test failed, but packages were created\n")
    }
  }
  
  # Step 5: Create instructions
  create_install_instructions(build_result, output_dir)
  cat("\n")
  
  # Summary
  if (!is.null(build_result)) {
    cat("=== BUILD SUMMARY ===\n")
    
    if (is.character(build_result)) {
      # Single package
      cat(sprintf("✓ Package created: %s\n", basename(build_result)))
      cat(sprintf("✓ Size: %.2f MB\n", file.info(build_result)$size / 1024 / 1024))
    } else {
      # List of packages
      if (!is.null(build_result$source)) {
        cat(sprintf("✓ Source package: %s\n", basename(build_result$source)))
        cat(sprintf("  Size: %.2f MB\n", file.info(build_result$source)$size / 1024 / 1024))
      }
      if (!is.null(build_result$binary)) {
        cat(sprintf("✓ Binary package: %s\n", basename(build_result$binary)))
        cat(sprintf("  Size: %.2f MB\n", file.info(build_result$binary)$size / 1024 / 1024))
      }
    }
    
    cat(sprintf("✓ Location: %s\n", output_dir))
    cat("\nNext steps:\n")
    cat("1. Test the packages on different systems using 'source' type\n")
    cat("2. If successful, push to GitHub to trigger automated builds\n")
    cat("3. Create a release to publish packages\n")
  } else {
    cat("=== BUILD FAILED ===\n")
    cat("Check the error messages above and fix issues before retrying\n")
  }
}

# Command line interface
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    main()
  } else if (args[1] == "--help") {
    cat("Usage: Rscript build_binary.R [output_dir] [--no-test]\n")
    cat("Options:\n")
    cat("  output_dir    Directory for binary output (default: binaries)\n")
    cat("  --no-test     Skip binary testing\n")
    cat("  --help        Show this help\n")
  } else {
    output_dir <- if (length(args) >= 1 && args[1] != "--no-test") args[1] else "binaries"
    test_install <- !("--no-test" %in% args)
    main(output_dir, test_install)
  }
} else {
  # Interactive mode - just run with defaults
  cat("Run main() to start binary building, or main(\"custom_dir\", FALSE) to skip testing\n")
}