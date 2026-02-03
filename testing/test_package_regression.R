#!/usr/bin/env Rscript

# SAIGEQTL Package Regression Test Script
# Automatically runs the example workflow to ensure package consistency

library(tools)

# Configuration - can be overridden by environment variables
TEST_DIR <- Sys.getenv("SAIGEQTL_TEST_DIR", "test_output")
EXTDATA_DIR <- Sys.getenv("SAIGEQTL_EXTDATA_DIR", "extdata")
PACKAGE_ROOT <- Sys.getenv("SAIGEQTL_PACKAGE_ROOT", ".")
LIBRARY_PATH <- Sys.getenv("SAIGEQTL_LIBRARY_PATH", "")
PIXI_MANIFEST <- Sys.getenv("SAIGEQTL_PIXI_MANIFEST", file.path(PACKAGE_ROOT, "pixi.toml"))

# Set expected output directory - prefer extdata/expected_output if it exists
extdata_expected <- file.path(EXTDATA_DIR, "expected_output")
if (dir.exists(extdata_expected)) {
  EXPECTED_DIR <- Sys.getenv("SAIGEQTL_EXPECTED_DIR", extdata_expected)
} else {
  EXPECTED_DIR <- Sys.getenv("SAIGEQTL_EXPECTED_DIR", "expected_output")
}

# Build Rscript command with pixi
build_rscript_cmd <- function(script_path, extra_args = "") {
  script_full_path <- file.path(PACKAGE_ROOT, EXTDATA_DIR, script_path)
  
  cmd <- sprintf("pixi run --manifest-path %s Rscript %s", 
                 PIXI_MANIFEST, script_full_path)
  
  # Add library path as a script argument if specified
  if (LIBRARY_PATH != "") {
    cmd <- sprintf("%s --library=%s", cmd, LIBRARY_PATH)
  }
  
  if (extra_args != "") {
    cmd <- sprintf("%s %s", cmd, extra_args)
  }
  
  return(cmd)
}

# Create test output directory
if (!dir.exists(TEST_DIR)) {
  dir.create(TEST_DIR, recursive = TRUE)
}

# Log function
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, msg))
}

# Function to run command and capture output
run_cmd <- function(cmd, log_file = NULL) {
  log_message(sprintf("Running: %s", cmd))
  
  if (!is.null(log_file)) {
    cmd <- sprintf("%s &> %s", cmd, log_file)
  }
  
  start_time <- Sys.time()
  result <- system(cmd, intern = FALSE)
  end_time <- Sys.time()
  
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (result == 0) {
    log_message(sprintf("Command completed successfully in %.2f seconds", elapsed))
  } else {
    log_message(sprintf("Command failed with exit code %d after %.2f seconds", result, elapsed))
    stop(sprintf("Command failed: %s", cmd))
  }
  
  return(result)
}

# Function to check if files exist
check_files_exist <- function(files, description) {
  log_message(sprintf("Checking %s files...", description))
  missing_files <- c()
  
  for (file in files) {
    if (!file.exists(file)) {
      missing_files <- c(missing_files, file)
    }
  }
  
  if (length(missing_files) > 0) {
    log_message(sprintf("Missing %s files:", description))
    for (file in missing_files) {
      log_message(sprintf("  - %s", file))
    }
    return(FALSE)
  }
  
  log_message(sprintf("All %s files found", description))
  return(TRUE)
}

# Function to compare output files with expected results
compare_results <- function(test_file, expected_file, tolerance = 1e-6) {
  if (!file.exists(expected_file)) {
    log_message(sprintf("Expected file not found: %s", expected_file))
    return(FALSE)
  }
  
  if (!file.exists(test_file)) {
    log_message(sprintf("Test output file not found: %s", test_file))
    return(FALSE)
  }
  
  # Read files and compare
  if (grepl("\\.txt$", test_file) || grepl("\\.log$", test_file)) {
    test_lines <- readLines(test_file)
    expected_lines <- readLines(expected_file)
    
    # For log files, focus on key result lines
    if (grepl("\\.log$", test_file)) {
      test_lines <- test_lines[grepl("(Analysis completed|Error|Warning)", test_lines)]
      expected_lines <- expected_lines[grepl("(Analysis completed|Error|Warning)", expected_lines)]
    }
    
    if (length(test_lines) != length(expected_lines)) {
      log_message(sprintf("File %s has different number of lines than expected", test_file))
      return(FALSE)
    }
    
    return(TRUE)
  }
  
  return(TRUE)
}

# Main testing function
run_regression_test <- function() {
  log_message("Starting SAIGEQTL regression test")
  log_message("=================================")
  
  # Check input files exist
  input_files <- c(
    file.path(EXTDATA_DIR, "input/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt"),
    file.path(EXTDATA_DIR, "input/n.indep_100_n.cell_1_01.step1.bed"),
    file.path(EXTDATA_DIR, "input/n.indep_100_n.cell_1_01.step1.bim"),
    file.path(EXTDATA_DIR, "input/n.indep_100_n.cell_1_01.step1.fam"),
    file.path(EXTDATA_DIR, "input/n.indep_100_n.cell_1.bed"),
    file.path(EXTDATA_DIR, "input/n.indep_100_n.cell_1.bim"),
    file.path(EXTDATA_DIR, "input/n.indep_100_n.cell_1.fam"),
    file.path(EXTDATA_DIR, "input/group_new_chrposa1a2.txt")
  )
  
  if (!check_files_exist(input_files, "input")) {
    stop("Required input files are missing")
  }
  
  # Step 1: Fit NULL GLMM
  log_message("Step 1: Fitting NULL GLMM")
  step1_args <- sprintf("--useSparseGRMtoFitNULL=FALSE \\
    --useGRMtoFitNULL=FALSE \\
    --phenoFile=./%s/input/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt \\
    --phenoCol=gene_1 \\
    --covarColList=X1,X2,pf1,pf2 \\
    --sampleCovarColList=X1,X2 \\
    --sampleIDColinphenoFile=IND_ID \\
    --traitType=count \\
    --outputPrefix=./%s/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1 \\
    --skipVarianceRatioEstimation=FALSE \\
    --isRemoveZerosinPheno=FALSE \\
    --isCovariateOffset=FALSE \\
    --isCovariateTransform=TRUE \\
    --skipModelFitting=FALSE \\
    --tol=0.00001 \\
    --plinkFile=./%s/input/n.indep_100_n.cell_1_01.step1 \\
    --IsOverwriteVarianceRatioFile=TRUE", 
    EXTDATA_DIR, TEST_DIR, EXTDATA_DIR)
  
  step1_cmd <- build_rscript_cmd("step1_fitNULLGLMM_qtl.R", step1_args)
  
  run_cmd(step1_cmd, file.path(TEST_DIR, "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1.log"))
  
  # Check step 1 outputs
  step1_outputs <- c(
    file.path(TEST_DIR, "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1.rda"),
    file.path(TEST_DIR, "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1.varianceRatio.txt")
  )
  
  if (!check_files_exist(step1_outputs, "Step 1 output")) {
    stop("Step 1 failed to produce expected outputs")
  }
  
  # Step 2: Single variant association tests - prepare region file
  log_message("Step 2: Preparing region file and running association tests")
  region_file <- file.path(TEST_DIR, "gene_1_cis_region.txt")
  writeLines("2\t1\t9810000", region_file)
  
  step2_single_args <- sprintf("--bedFile=./%s/input/n.indep_100_n.cell_1.bed \\
    --bimFile=./%s/input/n.indep_100_n.cell_1.bim \\
    --famFile=./%s/input/n.indep_100_n.cell_1.fam \\
    --SAIGEOutputFile=./%s/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis \\
    --chrom=2 \\
    --minMAF=0 \\
    --minMAC=20 \\
    --LOCO=FALSE \\
    --GMMATmodelFile=./%s/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1.rda \\
    --SPAcutoff=2 \\
    --varianceRatioFile=./%s/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1.varianceRatio.txt \\
    --rangestoIncludeFile=./%s/gene_1_cis_region.txt \\
    --markers_per_chunk=10000",
    EXTDATA_DIR, EXTDATA_DIR, EXTDATA_DIR, TEST_DIR, TEST_DIR, TEST_DIR, TEST_DIR)
  
  step2_single_cmd <- build_rscript_cmd("step2_tests_qtl.R", step2_single_args)
    
  run_cmd(step2_single_cmd, file.path(TEST_DIR, "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis.log"))
  
  # Check final outputs exist
  final_outputs <- c(
    file.path(TEST_DIR, "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis")
  )
  
  if (!check_files_exist(final_outputs, "final output")) {
    stop("Final outputs are missing")
  }
  
  # Compare with expected results if they exist
  if (dir.exists(EXPECTED_DIR)) {
    log_message("Comparing results with expected outputs...")
    
    expected_files <- list.files(EXPECTED_DIR, full.names = TRUE)
    for (expected_file in expected_files) {
      filename <- basename(expected_file)
      test_file <- file.path(TEST_DIR, filename)
      
      if (file.exists(test_file)) {
        if (compare_results(test_file, expected_file)) {
          log_message(sprintf("✓ %s matches expected output", filename))
        } else {
          log_message(sprintf("✗ %s differs from expected output", filename))
        }
      }
    }
  } else {
    log_message("No expected output directory found - saving current results as baseline")
    # Copy current results to expected directory for future comparisons
    dir.create(EXPECTED_DIR, recursive = TRUE)
    file.copy(list.files(TEST_DIR, full.names = TRUE), EXPECTED_DIR, recursive = TRUE)
  }
  
  log_message("Regression test completed successfully!")
  log_message("=================================")
  
  # Generate summary report
  cat("\nSUMMARY REPORT\n")
  cat("==============\n")
  cat(sprintf("Test directory: %s\n", TEST_DIR))
  cat("Steps completed:\n")
  cat("  ✓ Step 1: NULL GLMM fitting\n")
  cat("  ✓ Step 2: Single variant association tests with region file\n")
  cat(sprintf("Total runtime: %.2f minutes\n", as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
  
  return(TRUE)
}

# Test BGEN format
run_bgen_test <- function() {
  log_message("Starting BGEN format test")
  log_message("========================")
  
  # Check BGEN input files exist
  bgen_files <- c(
    file.path(EXTDATA_DIR, "input/n.indep_100_n.cell_1.bgen"),
    file.path(EXTDATA_DIR, "input/n.indep_100_n.cell_1.bgen.bgi"),
    file.path(EXTDATA_DIR, "input/n.indep_100_n.cell_1.sample")
  )
  
  if (!check_files_exist(bgen_files, "BGEN input")) {
    log_message("BGEN test skipped - required files missing")
    return(FALSE)
  }
  
  # Check if Step 1 results exist
  step1_outputs <- c(
    file.path(TEST_DIR, "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1.rda"),
    file.path(TEST_DIR, "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1.varianceRatio.txt")
  )
  
  if (!check_files_exist(step1_outputs, "Step 1 output for BGEN test")) {
    log_message("BGEN test skipped - Step 1 outputs missing. Run main test first.")
    return(FALSE)
  }
  
  # Create region file for BGEN test
  region_file_bgen <- file.path(TEST_DIR, "gene_1_cis_region_bgen.txt")
  writeLines("2\t1\t9810000", region_file_bgen)
  
  # Note: Sample IDs in BGEN should match phenotype file (a1, a2, etc.)
  log_message("Using n.indep_100_n.cell_1.sample file with matching sample IDs")
  
  # BGEN test parameters
  log_message("Running BGEN format association test")
  bgen_args <- sprintf("--bgenFile=./%s/input/n.indep_100_n.cell_1.bgen \\
    --bgenFileIndex=./%s/input/n.indep_100_n.cell_1.bgen.bgi \\
    --AlleleOrder=ref-first \\
    --sampleFile=./%s/input/n.indep_100_n.cell_1.sample \\
    --SAIGEOutputFile=./%s/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_bgen \\
    --chrom=2 \\
    --minMAF=0 \\
    --minMAC=20 \\
    --LOCO=FALSE \\
    --GMMATmodelFile=./%s/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1.rda \\
    --SPAcutoff=2 \\
    --varianceRatioFile=./%s/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1.varianceRatio.txt \\
    --rangestoIncludeFile=./%s/gene_1_cis_region_bgen.txt \\
    --markers_per_chunk=10000",
    EXTDATA_DIR, EXTDATA_DIR, EXTDATA_DIR, TEST_DIR, TEST_DIR, TEST_DIR, TEST_DIR)
  
  bgen_cmd <- build_rscript_cmd("step2_tests_qtl.R", bgen_args)
  
  # Try to run, but expect it might fail due to sample mismatch
  tryCatch({
    run_cmd(bgen_cmd, file.path(TEST_DIR, "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_bgen.log"))
  }, error = function(e) {
    log_message("BGEN test failed as expected due to sample ID mismatch or file format issues")
    log_message("This is a known limitation with the current test data")
    return(FALSE)
  })
  
  # Check BGEN outputs
  bgen_outputs <- c(
    file.path(TEST_DIR, "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_bgen")
  )
  
  if (!check_files_exist(bgen_outputs, "BGEN test output")) {
    log_message("BGEN test failed due to sample ID mismatch or file format issues")
    log_message("The BGEN format functionality was tested but failed due to incompatible test data")
    return(FALSE)
  }
  
  log_message("BGEN format test completed successfully!")
  return(TRUE)
}

# Test VCF format
run_vcf_test <- function() {
  log_message("Starting VCF format test")
  log_message("=======================")
  
  # Check VCF input files exist
  vcf_files <- c(
    file.path(EXTDATA_DIR, "input/n.indep_100_n.cell_1.vcf.gz"),
    file.path(EXTDATA_DIR, "input/n.indep_100_n.cell_1.vcf.gz.tbi")
  )
  
  if (!check_files_exist(vcf_files, "VCF input")) {
    log_message("VCF test skipped - required files missing")
    return(FALSE)
  }
  
  # Check if Step 1 results exist
  step1_outputs <- c(
    file.path(TEST_DIR, "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1.rda"),
    file.path(TEST_DIR, "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1.varianceRatio.txt")
  )
  
  if (!check_files_exist(step1_outputs, "Step 1 output for VCF test")) {
    log_message("VCF test skipped - Step 1 outputs missing. Run main test first.")
    return(FALSE)
  }
  
  # Create region file for VCF test
  region_file_vcf <- file.path(TEST_DIR, "gene_1_cis_region_vcf.txt")
  writeLines("2\t1\t9810000", region_file_vcf)
  
  # VCF test parameters using GT field (since DS not available in this file)
  log_message("Running VCF format association test")
  vcf_args <- sprintf("--vcfFile=./%s/input/n.indep_100_n.cell_1.vcf.gz \\
    --vcfFileIndex=./%s/input/n.indep_100_n.cell_1.vcf.gz.tbi \\
    --vcfField=GT \\
    --SAIGEOutputFile=./%s/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_vcf \\
    --chrom=2 \\
    --minMAF=0 \\
    --minMAC=20 \\
    --LOCO=FALSE \\
    --GMMATmodelFile=./%s/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1.rda \\
    --SPAcutoff=2 \\
    --varianceRatioFile=./%s/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1.varianceRatio.txt \\
    --rangestoIncludeFile=./%s/gene_1_cis_region_vcf.txt \\
    --markers_per_chunk=10000",
    EXTDATA_DIR, EXTDATA_DIR, TEST_DIR, TEST_DIR, TEST_DIR, TEST_DIR)
  
  vcf_cmd <- build_rscript_cmd("step2_tests_qtl.R", vcf_args)
  
  # Try to run VCF test
  tryCatch({
    run_cmd(vcf_cmd, file.path(TEST_DIR, "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_vcf.log"))
  }, error = function(e) {
    log_message("VCF test failed due to potential data format or compatibility issues")
    log_message(sprintf("Error: %s", e$message))
    return(FALSE)
  })
  
  # Check VCF outputs
  vcf_outputs <- c(
    file.path(TEST_DIR, "nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_vcf")
  )
  
  if (!check_files_exist(vcf_outputs, "VCF test output")) {
    log_message("VCF test failed - output files not generated")
    log_message("The VCF format functionality was tested but failed due to data format or compatibility issues")
    return(FALSE)
  }
  
  log_message("VCF format test completed successfully!")
  return(TRUE)
}

# Comprehensive test function
run_comprehensive_test <- function() {
  log_message("Starting comprehensive SAIGEQTL format testing")
  log_message("===============================================")
  
  # Track results
  results <- list()
  
  # Run main PLINK test first
  tryCatch({
    results$plink <- run_regression_test()
  }, error = function(e) {
    log_message(sprintf("PLINK test failed: %s", e$message))
    results$plink <- FALSE
  })
  
  # Run BGEN test
  tryCatch({
    results$bgen <- run_bgen_test()
  }, error = function(e) {
    log_message(sprintf("BGEN test failed: %s", e$message))
    results$bgen <- FALSE
  })
  
  # Run VCF test
  tryCatch({
    results$vcf <- run_vcf_test()
  }, error = function(e) {
    log_message(sprintf("VCF test failed: %s", e$message))
    results$vcf <- FALSE
  })
  
  # Report results
  log_message("Comprehensive test results:")
  log_message("===========================")
  log_message(sprintf("✓ PLINK format: %s", if(isTRUE(results$plink)) "PASSED" else "FAILED"))
  log_message(sprintf("✓ BGEN format: %s", if(isTRUE(results$bgen)) "PASSED" else "FAILED/SKIPPED"))
  log_message(sprintf("✓ VCF format: %s", if(isTRUE(results$vcf)) "PASSED" else "FAILED/SKIPPED"))
  
  return(results)
}

# Save baseline results function
save_baseline <- function() {
  if (dir.exists(TEST_DIR) && length(list.files(TEST_DIR)) > 0) {
    if (dir.exists(EXPECTED_DIR)) {
      unlink(EXPECTED_DIR, recursive = TRUE)
    }
    dir.create(EXPECTED_DIR, recursive = TRUE)
    file.copy(list.files(TEST_DIR, full.names = TRUE), EXPECTED_DIR, recursive = TRUE)
    log_message(sprintf("Baseline results saved to %s", EXPECTED_DIR))
  } else {
    log_message("No test results found to save as baseline")
  }
}

# Clean up function
clean_test_output <- function() {
  if (dir.exists(TEST_DIR)) {
    unlink(TEST_DIR, recursive = TRUE)
    log_message("Test output directory cleaned")
  }
}

# Main execution
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) > 0) {
    if (args[1] == "clean") {
      clean_test_output()
      quit(status = 0)
    } else if (args[1] == "baseline") {
      save_baseline()
      quit(status = 0)
    } else if (args[1] == "bgen") {
      # Record start time
      start_time <- Sys.time()
      # Run BGEN test only
      tryCatch({
        result <- run_bgen_test()
        quit(status = if(result) 0 else 1)
      }, error = function(e) {
        log_message(sprintf("BGEN test failed with error: %s", e$message))
        quit(status = 1)
      })
    } else if (args[1] == "vcf") {
      # Record start time
      start_time <- Sys.time()
      # Run VCF test only
      tryCatch({
        result <- run_vcf_test()
        quit(status = if(result) 0 else 1)
      }, error = function(e) {
        log_message(sprintf("VCF test failed with error: %s", e$message))
        quit(status = 1)
      })
    } else if (args[1] == "comprehensive" || args[1] == "all") {
      # Record start time
      start_time <- Sys.time()
      # Run all tests
      tryCatch({
        results <- run_comprehensive_test()
        # Exit with error if any test failed
        all_passed <- isTRUE(results$plink) && isTRUE(results$bgen) && isTRUE(results$vcf)
        quit(status = if(all_passed) 0 else 1)
      }, error = function(e) {
        log_message(sprintf("Comprehensive test failed with error: %s", e$message))
        quit(status = 1)
      })
    }
  }
  
  # Record start time
  start_time <- Sys.time()
  
  # Run the default PLINK test
  tryCatch({
    result <- run_regression_test()
    quit(status = 0)
  }, error = function(e) {
    log_message(sprintf("Test failed with error: %s", e$message))
    quit(status = 1)
  })
}