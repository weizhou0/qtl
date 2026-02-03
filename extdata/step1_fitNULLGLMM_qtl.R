#!/usr/bin/env -S pixi run --manifest-path ../pixi.toml Rscript

options(stringsAsFactors = F)

## load R libraries
require(optparse) # install.packages("optparse")

## set list of cmd line arguments
option_list <- list(
  make_option("--plinkFile",
    type = "character", default = "",
    help = "Path to plink file for creating the genetic relationship matrix (GRM). minMAFforGRM can be used to specify the minimum MAF and maxMissingRate can be used to specify the maximum missing rates  of markers in the plink file to be used for constructing GRM. Genetic markers are also randomly selected from the plink file to estimate the variance ratios"
  ),
  make_option("--bedFile",
    type = "character", default = "",
    help = "Path to bed file. If plinkFile is specified, 'plinkFile'.bed will be used"
  ),
  make_option("--bimFile",
    type = "character", default = "",
    help = "Path to bim file. If plinkFile is specified, 'plinkFile'.bim will be used"
  ),
  make_option("--famFile",
    type = "character", default = "",
    help = "Path to fam file. If plinkFile is specified, 'plinkFile'.fam will be used"
  ),
  make_option("--phenoFile",
    type = "character", default = "",
    help = "Required. Path to the phenotype file. The file can be either tab or space delimited. The phenotype file has a header and contains at least two columns. One column is for phentoype and the other column is for sample IDs. Additional columns can be included in the phenotype file for covariates in the null model. Please specify the names of the covariates using the argument covarColList and specify categorical covariates using the argument qCovarColList. All categorical covariates must also be included in covarColList."
  ),
  make_option("--phenoCol",
    type = "character", default = "",
    help = "Required. Column name for phenotype to be tested in the phenotype file, e.g CAD"
  ),
  make_option("--isRemoveZerosinPheno",
    type = "logical", default = FALSE,
    help = "Optional. Whether to remove zeros in the phenotype"
  ),
  make_option("--traitType", type = "character", default = "binary", help = "Required. binary or quantitative [default=binary]"),
  make_option("--invNormalize",
    type = "logical", default = FALSE,
    help = "Optional. Only for quantitative. Whether to perform the inverse normalization for the phenotype [default='FALSE']"
  ),
  make_option("--covarColList",
    type = "character", default = "",
    help = "List of covariates (comma separated)"
  ),
  make_option("--sampleCovarColList",
    type = "character", default = "",
    help = "List of covariates that are on sample level (comma separated)"
  ),
  make_option("--longlCol",
    type = "character", default = "",
    help = ""
  ),
  make_option("--qCovarColList",
    type = "character", default = "",
    help = "List of categorical covariates (comma separated). All categorical covariates must also be in covarColList"
  ),
  make_option("--sampleIDColinphenoFile",
    type = "character", default = "IID",
    help = "Required. Column name of sample IDs in the phenotype file, e.g. IID"
  ),
  make_option("--cellIDColinphenoFile",
    type = "character", default = "",
    help = "Column name of cell IDs in the phenotype file, e.g. barcode"
  ),
  make_option("--tol",
    type = "numeric", default = 0.02,
    help = "Optional. Tolerance for fitting the null GLMM to converge [default=0.02]."
  ),
  make_option("--maxiter",
    type = "integer", default = 20,
    help = "Optional. Maximum number of iterations used to fit the null GLMM [default=20]."
  ),
  make_option("--tolPCG",
    type = "numeric", default = 1e-5,
    help = "Optional. Tolerance for PCG to converge [default=1e-5]."
  ),
  make_option("--maxiterPCG",
    type = "integer", default = 500,
    help = "Optional. Maximum number of iterations for PCG [default=500]."
  ),
  make_option("--nThreads",
    type = "integer", default = 1,
    help = "Optional. Number of threads (CPUs) to use [default=1]."
  ),
  make_option("--SPAcutoff",
    type = "numeric", default = 2,
    help = "Optional. Cutoff for the deviation of score test statistics from mean in the unit of sd to perform SPA [default=2]."
  ),
  make_option("--numRandomMarkerforVarianceRatio",
    type = "integer", default = 30,
    help = "Optional. An integer greater than 0. Number of markers to be randomly selected for estimating the variance ratio. The number will be automatically added by 10 until the coefficient of variantion (CV) for the variance ratio estimate is below ratioCVcutoff [default=30]."
  ),
  make_option("--skipModelFitting",
    type = "logical", default = FALSE,
    help = "Optional. Whether to skip model fitting and only to estimate the variance ratio. If TRUE, the file outputPrefix.rda is required [default='FALSE']"
  ),
  make_option("--skipVarianceRatioEstimation",
    type = "logical", default = FALSE,
    help = "Optional. Whether to skip model fitting and only to estimate the variance ratio. If TRUE, the file outputPrefix.rda is required [default='FALSE']"
  ),
  make_option("--memoryChunk",
    type = "numeric", default = 2,
    help = "Optional. Size (Gb) for each memory chunk [default=2]"
  ),
  make_option("--tauInit",
    type = "character", default = "0,0",
    help = "Optional. Initial values for tau. [default=0,0]"
  ),
  make_option("--LOCO",
    type = "logical", default = TRUE,
    help = "Whether to apply the leave-one-chromosome-out (LOCO) approach when fitting the null model using the full GRM [default=TRUE]."
  ),
  make_option("--isLowMemLOCO",
    type = "logical", default = FALSE,
    help = "Whehter to output the model file by chromosome when LOCO=TRUE. If TRUE, the memory usage in Step 1 and Step 2 will be lower [default=FALSE]"
  ),
  make_option("--traceCVcutoff",
    type = "numeric", default = 0.0025,
    help = "Optional. Threshold for coefficient of variation (CV) for the trace estimator. Number of runs for trace estimation will be increased until the CV is below the threshold [default=0.0025]."
  ),
  make_option("--nrun",
    type = "numeric", default = 30,
    help = "Number of rums in trace estimation. [default=30]"
  ),
  make_option("--ratioCVcutoff",
    type = "numeric", default = 0.001,
    help = "Optional. Threshold for coefficient of variation (CV) for estimating the variance ratio. The number of randomly selected markers will be increased until the CV is below the threshold [default=0.001]"
  ),
  make_option("--outputPrefix",
    type = "character", default = "~/",
    help = "Required. Path and prefix of the output files [default='~/']"
  ),
  make_option("--outputPrefix_varRatio",
    type = "character", default = "",
    help = "Optional. Path and prefix of the output the variance ratio file. if not specified, it will be the same as the outputPrefix"
  ),
  make_option("--IsOverwriteVarianceRatioFile",
    type = "logical", default = FALSE,
    help = "Optional. Whether to overwrite the variance ratio file if the file exist.[default='FALSE']"
  ),
  make_option("--sparseGRMFile",
    type = "character", default = "",
    help = "Path to the pre-calculated sparse GRM file. If not specified and  IsSparseKin=TRUE, sparse GRM will be computed [default=NULL]"
  ),
  make_option("--sparseGRMSampleIDFile",
    type = "character", default = "",
    help = "Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to sample IDs in the sparse GRM [default=NULL]"
  ),
  make_option("--isCateVarianceRatio",
    type = "logical", default = FALSE,
    help = "Required. Whether to estimate variance ratio based on different MAC categories. If yes, variance ratio will be estiamted for multiple MAC categories corresponding to cateVarRatioMinMACVecExclude and cateVarRatioMaxMACVecInclude. Currently, if isCateVarianceRatio=TRUE, then LOCO=FALSE [default=FALSE]"
  ),
  make_option("--relatednessCutoff",
    type = "numeric", default = 0,
    help = "Optional. Threshold (minimum relatedness coefficient) to treat two samples as unrelated when the sparse GRM is used [default=0]"
  ),
  make_option("--cateVarRatioMinMACVecExclude",
    type = "character", default = "10,20.5",
    help = "Optional. vector of float. Lower bound for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. [default='10,20.5']"
  ),
  make_option("--cateVarRatioMaxMACVecInclude",
    type = "character", default = "20.5",
    help = "Optional. vector of float. Higher bound for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. [default='20.5']"
  ),
  make_option("--isCovariateTransform",
    type = "logical", default = TRUE,
    help = "Optional. Whether use qr transformation on covariates [default='TRUE']."
  ),
  make_option("--isDiagofKinSetAsOne",
    type = "logical", default = FALSE,
    help = "Optional. Whether to set the diagnal elements in GRM to be 1 [default='FALSE']."
  ),
  make_option("--useSparseGRMtoFitNULL",
    type = "logical", default = FALSE,
    help = "Optional. Whether to use sparse GRM to fit the null model [default='FALSE']."
  ),
  make_option("--useSparseGRMforVarRatio",
    type = "logical", default = FALSE,
    help = "Optional. Whether to use sparse GRM to estimate the variance Ratios. If TRUE, the variance ratios will be estimated using the full GRM (numerator) and the sparse GRM (denominator). By default, FALSE"
  ),
  make_option("--minMAFforGRM",
    type = "numeric", default = 0.01,
    help = "Optional. Minimum MAF of markers used for GRM"
  ),
  make_option("--maxMissingRateforGRM",
    type = "numeric", default = 0.15,
    help = "Optional. Maximum missing rate of markers used for GRM"
  ),
  make_option("--minCovariateCount",
    type = "numeric", default = -1,
    help = "Optional. Binary covariates with a count less than minCovariateCount will be excluded from the model to avoid convergence issues [default=-1] (no covariates will be excluded)."
  ),
  make_option("--includeNonautoMarkersforVarRatio",
    type = "logical", default = FALSE,
    help = "Optional. Whether to allow for non-autosomal markers for variance ratio. [default, 'FALSE']"
  ),
  make_option("--FemaleOnly",
    type = "logical", default = FALSE,
    help = "Optional. Whether to run Step 1 for females only [default=FALSE]. if TRUE, --sexCol and --FemaleCode need to be specified"
  ),
  make_option("--MaleOnly",
    type = "logical", default = FALSE,
    help = "Optional. Whether to run Step 1 for males only [default=FALSE]. if TRUE, --sexCol and --MaleCode need to be specified"
  ),
  make_option("--sexCol",
    type = "character", default = "",
    help = "Optional. Column name for sex in the phenotype file, e.g Sex"
  ),
  make_option("--FemaleCode",
    type = "character", default = "1",
    help = "Optional. Values in the column for sex in the phenotype file are used for females [default, '1']"
  ),
  make_option("--MaleCode",
    type = "character", default = "0",
    help = "Optional. Values in the column for sex in the phenotype file are used for males [default, '0']"
  ),
  make_option("--isCovariateOffset",
    type = "logical", default = TRUE,
    help = "Optional. Whether to estimate fixed effect coeffciets. [default, 'TRUE']"
  ),
  make_option("--SampleIDIncludeFile",
    type = "character", default = "",
    help = "Path to the file that contains one column for IDs of samples who will be include for null model fitting."
  ),
  make_option("--VmatFilelist",
    type = "character", default = "",
    help = "List of additional V (comma separated)"
  ),
  make_option("--VmatSampleFilelist",
    type = "character", default = "",
    help = "List of additional V (comma separated)"
  ),
  make_option("--useGRMtoFitNULL", type = "logical", default = TRUE, help = ""),
  make_option("--offsetCol",
    type = "character", default = NULL,
    help = "offset column"
  ),
  make_option("--varWeightsCol",
    type = "character", default = NULL,
    help = "variance weight column"
  ),
  make_option("--isStoreSigma",
    type = "logical", default = TRUE,
    help = "Optional. Whether to store the inv Sigma matrix. [default, 'TRUE']"
  ),
  make_option("--isShrinkModelOutput",
    type = "logical", default = TRUE,
    help = "Optional. Whether to remove unnecessary objects for step2 from the model output. [default, 'TRUE']"
  ),
  make_option("--isExportResiduals",
    type = "logical", default = FALSE,
    help = "Optional. Whether to export residual vector. [default, 'FALSE']"
  ),
  make_option("--varRatioBatchSize",
    type = "integer", default = 1,
    help = "Optional. Batch size for variance ratio estimation. Higher values use more memory but may be faster. Set to 1 for sequential processing (low memory), or >1 for batch processing. [default=1]"
  ),
  make_option("--solverMethod",
    type = "character", default = "auto",
    help = "Optional. Solver method for fitting null model: 'auto' (automatic selection based on data structure), 'pcg' (Preconditioned Conjugate Gradient), 'smw' (Sherman-Morrison-Woodbury). When 'auto': SMW for repeated cells without GRM, PCG when GRM is provided or no repeated cells. [default='auto']"
  ),
  make_option("--library",
    type = "character", default = "",
    help = "Optional. Path to custom R library directory where SAIGEQTL package is installed. If not specified, uses default R library paths."
  )
)

# Parse options to get library path
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

# Set custom library path if provided
if (!is.null(opt$library) && opt$library != "") {
  .libPaths(c(opt$library, .libPaths()))
  cat("Using custom library path:", opt$library, "\n")
}

if (!is.null(opt$library) && opt$library != "") {
  library(SAIGEQTL, lib.loc = opt$library)
} else {
  library(SAIGEQTL)
}

print(sessionInfo())

print(opt)

covars <- strsplit(opt$covarColList, ",")[[1]]
qcovars <- strsplit(opt$qCovarColList, ",")[[1]]
scovars <- strsplit(opt$sampleCovarColList, ",")[[1]]
convertoNumeric <- function(x, stringOutput) {
  y <- tryCatch(expr = as.numeric(x), warning = function(w) {
    return(NULL)
  })
  if (is.null(y)) {
    stop(stringOutput, " is not numeric\n")
  } else {
    cat(stringOutput, " is ", y, "\n")
  }
  return(y)
}

tauInit <- convertoNumeric(strsplit(opt$tauInit, ",")[[1]], "tauInit")
cateVarRatioMinMACVecExclude <- convertoNumeric(x = strsplit(opt$cateVarRatioMinMACVecExclude, ",")[[1]], "cateVarRatioMinMACVecExclude")
cateVarRatioMaxMACVecInclude <- convertoNumeric(x = strsplit(opt$cateVarRatioMaxMACVecInclude, ",")[[1]], "cateVarRatioMaxMACVecInclude")

BLASctl_installed <- require(RhpcBLASctl)
if (BLASctl_installed) {
  # Set number of threads for BLAS to 1, this step does not benefit from multithreading or multiprocessing
  original_num_threads <- blas_get_num_procs()
  blas_set_num_threads(1)
}


# Check if varRatioBatchSize parameter is supported in this version
varRatioBatchSize_supported <- "varRatioBatchSize" %in% names(formals(fitNULLGLMM_multiV))

# Debug: Print some key variables
cat("covars:", covars, "\n")
cat("qcovars:", qcovars, "\n") 
cat("scovars:", scovars, "\n")

# Prepare arguments list - remove empty/NULL values as we build it
args_list <- list()

# Add non-empty arguments only
if (!is.null(opt$plinkFile) && opt$plinkFile != "") args_list$plinkFile <- opt$plinkFile
if (!is.null(opt$bedFile) && opt$bedFile != "") args_list$bedFile <- opt$bedFile
if (!is.null(opt$bimFile) && opt$bimFile != "") args_list$bimFile <- opt$bimFile
if (!is.null(opt$famFile) && opt$famFile != "") args_list$famFile <- opt$famFile
args_list$useSparseGRMtoFitNULL <- opt$useSparseGRMtoFitNULL
if (!is.null(opt$sparseGRMFile) && opt$sparseGRMFile != "") args_list$sparseGRMFile <- opt$sparseGRMFile
if (!is.null(opt$sparseGRMSampleIDFile) && opt$sparseGRMSampleIDFile != "") args_list$sparseGRMSampleIDFile <- opt$sparseGRMSampleIDFile
args_list$phenoFile <- opt$phenoFile
args_list$phenoCol <- opt$phenoCol
args_list$isRemoveZerosinPheno <- opt$isRemoveZerosinPheno
args_list$sampleIDColinphenoFile <- opt$sampleIDColinphenoFile
if (!is.null(opt$cellIDColinphenoFile) && opt$cellIDColinphenoFile != "") args_list$cellIDColinphenoFile <- opt$cellIDColinphenoFile
args_list$traitType <- opt$traitType
args_list$outputPrefix <- opt$outputPrefix
args_list$isCovariateOffset <- opt$isCovariateOffset
args_list$nThreads <- opt$nThreads
args_list$useSparseGRMforVarRatio <- opt$useSparseGRMforVarRatio
args_list$invNormalize <- opt$invNormalize
if (!is.null(covars) && length(covars) > 0 && !all(covars == "")) args_list$covarColList <- covars
if (!is.null(qcovars) && length(qcovars) > 0 && !all(qcovars == "")) args_list$qCovarCol <- qcovars
args_list$tol <- opt$tol
args_list$maxiter <- opt$maxiter
args_list$tolPCG <- opt$tolPCG
args_list$maxiterPCG <- opt$maxiterPCG
args_list$SPAcutoff <- opt$SPAcutoff
args_list$numMarkersForVarRatio <- opt$numRandomMarkerforVarianceRatio
args_list$skipModelFitting <- opt$skipModelFitting
args_list$skipVarianceRatioEstimation <- opt$skipVarianceRatioEstimation
args_list$memoryChunk <- opt$memoryChunk
args_list$tauInit <- tauInit
args_list$LOCO <- opt$LOCO
args_list$isLowMemLOCO <- opt$isLowMemLOCO
args_list$traceCVcutoff <- opt$traceCVcutoff
args_list$nrun <- opt$nrun
args_list$ratioCVcutoff <- opt$ratioCVcutoff
if (!is.null(opt$outputPrefix_varRatio) && opt$outputPrefix_varRatio != "") args_list$outputPrefix_varRatio <- opt$outputPrefix_varRatio
args_list$IsOverwriteVarianceRatioFile <- opt$IsOverwriteVarianceRatioFile
args_list$relatednessCutoff <- opt$relatednessCutoff
args_list$isCateVarianceRatio <- opt$isCateVarianceRatio
args_list$cateVarRatioMinMACVecExclude <- cateVarRatioMinMACVecExclude
args_list$cateVarRatioMaxMACVecInclude <- cateVarRatioMaxMACVecInclude
args_list$isCovariateTransform <- opt$isCovariateTransform
args_list$isDiagofKinSetAsOne <- opt$isDiagofKinSetAsOne
args_list$minMAFforGRM <- opt$minMAFforGRM
args_list$maxMissingRateforGRM <- opt$maxMissingRateforGRM
args_list$minCovariateCount <- opt$minCovariateCount
args_list$includeNonautoMarkersforVarRatio <- opt$includeNonautoMarkersforVarRatio
if (!is.null(opt$sexCol) && opt$sexCol != "") args_list$sexCol <- opt$sexCol
if (!is.null(opt$FemaleCode) && opt$FemaleCode != "") args_list$FemaleCode <- opt$FemaleCode
args_list$FemaleOnly <- opt$FemaleOnly
if (!is.null(opt$MaleCode) && opt$MaleCode != "") args_list$MaleCode <- opt$MaleCode
args_list$MaleOnly <- opt$MaleOnly
if (!is.null(opt$SampleIDIncludeFile) && opt$SampleIDIncludeFile != "") args_list$SampleIDIncludeFile <- opt$SampleIDIncludeFile
if (!is.null(opt$VmatFilelist) && opt$VmatFilelist != "") args_list$VmatFilelist <- opt$VmatFilelist
if (!is.null(opt$VmatSampleFilelist) && opt$VmatSampleFilelist != "") args_list$VmatSampleFilelist <- opt$VmatSampleFilelist
if (!is.null(opt$longlCol) && opt$longlCol != "") args_list$longlCol <- opt$longlCol
args_list$useGRMtoFitNULL <- opt$useGRMtoFitNULL
if (!is.null(opt$offsetCol) && opt$offsetCol != "") args_list$offsetCol <- opt$offsetCol
if (!is.null(opt$varWeightsCol) && opt$varWeightsCol != "") args_list$varWeightsCol <- opt$varWeightsCol
if (!is.null(scovars) && length(scovars) > 0 && !all(scovars == "")) args_list$sampleCovarCol <- scovars
args_list$isStoreSigma <- opt$isStoreSigma
args_list$isShrinkModelOutput <- opt$isShrinkModelOutput
args_list$isExportResiduals <- opt$isExportResiduals

# Conditionally add varRatioBatchSize if supported
if (varRatioBatchSize_supported) {
  args_list$varRatioBatchSize <- opt$varRatioBatchSize
}

# Add solverMethod if function supports it
solverMethod_supported <- "solverMethod" %in% names(formals(fitNULLGLMM_multiV))
if (solverMethod_supported) {
  args_list$solverMethod <- opt$solverMethod
}

cat("Final number of arguments:", length(args_list), "\n")

# Function to remove all output files when convergence fails
remove_failed_output_files <- function(output_prefix, output_prefix_varRatio = "") {
  files_to_remove <- c(
    paste0(output_prefix, ".rda"),
    paste0(output_prefix, ".varianceRatio.txt")
  )
  
  # Add variance ratio file with alternative prefix if specified
  if (output_prefix_varRatio != "" && output_prefix_varRatio != output_prefix) {
    files_to_remove <- c(files_to_remove, paste0(output_prefix_varRatio, ".varianceRatio.txt"))
  }
  
  removed_files <- c()
  for (file_path in files_to_remove) {
    if (file.exists(file_path)) {
      tryCatch({
        file.remove(file_path)
        removed_files <- c(removed_files, file_path)
        cat(sprintf("Removed failed output file: %s\n", file_path))
      }, error = function(e) {
        cat(sprintf("Warning: Could not remove file %s: %s\n", file_path, e$message))
      })
    }
  }
  
  if (length(removed_files) > 0) {
    cat(sprintf("Cleaned up %d failed output files\n", length(removed_files)))
  }
  
  return(removed_files)
}

# Function to check model convergence status and write status file
check_convergence_and_write_status <- function(output_prefix, args_used, package_version, final_status = NULL) {
  status_file <- paste0(output_prefix, ".status.txt")
  
  # Initialize status information
  status_info <- list(
    timestamp = Sys.time(),
    package_version = package_version,
    convergence_status = "UNKNOWN",
    convergence_details = "",
    final_model_file = "",
    variance_ratio_file = "",
    arguments_used = args_used
  )
  
  # Check if model file exists and examine convergence
  rda_file <- paste0(output_prefix, ".rda")
  if (file.exists(rda_file)) {
    tryCatch({
      my_env <- new.env()
      load(rda_file, envir = my_env)
      modglmm <- my_env$modglmm
      
      if (!is.null(modglmm)) {
        # Check convergence flag
        model_converged <- modglmm$converged
        
        # Check variance component bounds
        theta_valid <- TRUE
        theta_sum <- sum(modglmm$theta[2:length(modglmm$theta)])
        if (theta_sum <= 0 || theta_sum > 1) {
          theta_valid <- FALSE
        }
        
        # Determine overall status
        if (model_converged && theta_valid) {
          status_info$convergence_status <- "SUCCESS"
          status_info$convergence_details <- sprintf("Model converged successfully. Theta sum: %.6f", theta_sum)
        } else if (!model_converged) {
          status_info$convergence_status <- "FAILED"
          status_info$convergence_details <- sprintf("Model failed to converge. Theta sum: %.6f", theta_sum)
        } else {
          status_info$convergence_status <- "FAILED"
          status_info$convergence_details <- sprintf("Variance components out of bounds. Theta sum: %.6f", theta_sum)
        }
        
        status_info$final_model_file <- rda_file
        
        # Check for variance ratio file
        var_ratio_file <- paste0(output_prefix, ".varianceRatio.txt")
        if (file.exists(var_ratio_file)) {
          status_info$variance_ratio_file <- var_ratio_file
        }
        
      } else {
        status_info$convergence_status <- "FAILED"
        status_info$convergence_details <- "Model object not found in .rda file"
      }
    }, error = function(e) {
      status_info$convergence_status <- "FAILED"
      status_info$convergence_details <- paste("Error loading model:", e$message)
    })
  } else {
    status_info$convergence_status <- "FAILED"
    status_info$convergence_details <- "Model .rda file not found"
  }
  
  # Override with final status if provided
  if (!is.null(final_status)) {
    status_info$convergence_status <- final_status$status
    status_info$convergence_details <- final_status$details
    # If final status indicates failure, ensure no output files remain
    if (final_status$status == "FAILED") {
      remove_failed_output_files(output_prefix, ifelse("outputPrefix_varRatio" %in% names(args_used), args_used$outputPrefix_varRatio, ""))
      status_info$final_model_file <- ""
      status_info$variance_ratio_file <- ""
    }
  }
  
  # Write status file
  cat("=== SAIGEQTL Step 1 Analysis Status ===\n", file = status_file)
  cat(sprintf("Timestamp: %s\n", status_info$timestamp), file = status_file, append = TRUE)
  cat(sprintf("SAIGEQTL Version: %s\n", status_info$package_version), file = status_file, append = TRUE)
  cat(sprintf("Convergence Status: %s\n", status_info$convergence_status), file = status_file, append = TRUE)
  cat(sprintf("Details: %s\n", status_info$convergence_details), file = status_file, append = TRUE)
  cat(sprintf("Final Model File: %s\n", ifelse(status_info$final_model_file != "", status_info$final_model_file, "None")), file = status_file, append = TRUE)
  cat(sprintf("Variance Ratio File: %s\n", ifelse(status_info$variance_ratio_file != "", status_info$variance_ratio_file, "None")), file = status_file, append = TRUE)
  cat("\n=== Arguments Used ===\n", file = status_file, append = TRUE)
  
  # Write all arguments
  for (arg_name in names(status_info$arguments_used)) {
    arg_value <- status_info$arguments_used[[arg_name]]
    if (is.vector(arg_value) && length(arg_value) > 1) {
      arg_value <- paste(arg_value, collapse = ",")
    }
    cat(sprintf("%s: %s\n", arg_name, arg_value), file = status_file, append = TRUE)
  }
  
  cat(sprintf("\nStatus file written to: %s\n", status_file))
  return(status_info)
}

# Get SAIGEQTL version
package_version <- tryCatch({
  as.character(packageVersion("SAIGEQTL"))
}, error = function(e) {
  "Unknown"
})

cat(sprintf("SAIGEQTL Version: %s\n", package_version))

# Call the main function
set.seed(1)

fit_success <- TRUE
initial_error <- NULL
tryCatch({
   do.call(fitNULLGLMM_multiV, args_list)
}, error = function(e) {
  message("Initial model failed with error: ", e$message)
  initial_error <<- e$message
  fit_success <<- FALSE  # Track failure
  # Clean up any partial output files created during failed attempt
  remove_failed_output_files(opt$outputPrefix, opt$outputPrefix_varRatio)
})

if(fit_success){
  # Check initial model results if we started without offset
  if (!opt$isCovariateOffset) {
    my_env <- new.env()
    load(paste0(opt$outputPrefix, ".rda"), envir = my_env)
    modglmm <- my_env$modglmm
    print(modglmm$theta)
    
    # Check if model failed or variance components are out of bounds
    theta_sum <- sum(modglmm$theta[2:length(modglmm$theta)])
    if (theta_sum <= 0 || theta_sum > 1 || !modglmm$converged) {
      cat("Initial model failed (convergence: ", modglmm$converged, ", theta sum: ", theta_sum, ")\n")
      cat("Retrying with all covariates as offset...\n")
      
      # Remove failed initial model files
      remove_failed_output_files(opt$outputPrefix, opt$outputPrefix_varRatio)
      
      # Retry with offset
      opt$isCovariateOffset <- TRUE
      args_list$isCovariateOffset <- opt$isCovariateOffset
      
      set.seed(1)
      tryCatch({
        do.call(fitNULLGLMM_multiV, args_list)
      }, error = function(e) {
        message("Offset model also failed with error: ", e$message)
        fit_success <<- FALSE
        # Clean up any partial output files from failed retry
        remove_failed_output_files(opt$outputPrefix, opt$outputPrefix_varRatio)
      })
    }
  } else {
    # If we started with offset, just check convergence and retry with more iterations if needed
    my_env <- new.env()
    load(paste0(opt$outputPrefix, ".rda"), envir = my_env)
    modglmm <- my_env$modglmm
    
    if(!modglmm$converged){
      cat("Model didn't converge successfully. Trying with increased maxiter (", opt$maxiter + 500, ")\n")
      opt$maxiter = opt$maxiter + 500
      args_list$maxiter <- opt$maxiter
      
      # Remove failed files before retry
      remove_failed_output_files(opt$outputPrefix, opt$outputPrefix_varRatio)
      
      tryCatch({
        do.call(fitNULLGLMM_multiV, args_list)
      }, error = function(e) {
        message("Retry with increased iterations also failed: ", e$message)
        fit_success <<- FALSE
        # Clean up any partial output files from failed retry
        remove_failed_output_files(opt$outputPrefix, opt$outputPrefix_varRatio)
      })
    }
  }
} else {
  # Initial fitting completely failed, try with offset
  cat("Initial model fitting failed completely\n")
  if (!opt$isCovariateOffset) {
    cat("Trying with all covariates as offset...\n")
    opt$isCovariateOffset <- TRUE
    args_list$isCovariateOffset <- opt$isCovariateOffset
    
    tryCatch({
      do.call(fitNULLGLMM_multiV, args_list)
      fit_success <<- TRUE  # Mark as successful if offset works
    }, error = function(e) {
      message("Offset model also failed with error: ", e$message)
      fit_success <<- FALSE
      # Clean up any partial output files from failed offset attempt
      remove_failed_output_files(opt$outputPrefix, opt$outputPrefix_varRatio)
    })
  }
}

# Final convergence check and cleanup
final_convergence_check <- function() {
  rda_file <- paste0(opt$outputPrefix, ".rda")
  if (file.exists(rda_file)) {
    tryCatch({
      my_env <- new.env()
      load(rda_file, envir = my_env)
      modglmm <- my_env$modglmm
      
      if (!is.null(modglmm)) {
        # Check both convergence flag and variance component bounds
        theta_sum <- sum(modglmm$theta[2:length(modglmm$theta)])
        theta_valid <- (theta_sum > 0 && theta_sum <= 1)
        
        if (!modglmm$converged || !theta_valid) {
          cat("Final check: Model convergence failed (converged:", modglmm$converged, ", theta_sum:", theta_sum, "). Removing output files.\n")
          remove_failed_output_files(opt$outputPrefix, opt$outputPrefix_varRatio)
          
          details_msg <- sprintf("Model failed (converged: %s, theta_sum: %.6f) and output files removed", 
                               modglmm$converged, theta_sum)
          return(list(status = "FAILED", details = details_msg))
        }
      } else {
        cat("Final check: Model object not found. Removing output files.\n")
        remove_failed_output_files(opt$outputPrefix, opt$outputPrefix_varRatio)
        return(list(status = "FAILED", details = "Model object not found and output files removed"))
      }
      return(NULL)
    }, error = function(e) {
      return(list(status = "FAILED", details = paste("Error in final check:", e$message)))
    })
  }
  return(NULL)
}

final_check_result <- final_convergence_check()

# Generate comprehensive status report
cat("\n=== Generating Final Status Report ===\n")

# Handle case where initial fitting completely failed
if (!fit_success && !is.null(initial_error)) {
  # Ensure all output files are cleaned up for complete failure
  remove_failed_output_files(opt$outputPrefix, opt$outputPrefix_varRatio)
  
  final_status <- list(
    status = "FAILED", 
    details = paste("Initial model fitting failed with error:", initial_error)
  )
  check_convergence_and_write_status(opt$outputPrefix, args_list, package_version, final_status)
} else {
  # Normal status checking
  status_result <- check_convergence_and_write_status(opt$outputPrefix, args_list, package_version, final_check_result)
  
  # Print final summary
  cat("\n=== FINAL SUMMARY ===\n")
  cat(sprintf("Analysis Status: %s\n", status_result$convergence_status))
  cat(sprintf("Details: %s\n", status_result$convergence_details))
  if (status_result$convergence_status == "SUCCESS") {
    cat("✓ Step 1 completed successfully!\n")
    cat(sprintf("✓ Model file: %s\n", status_result$final_model_file))
    if (status_result$variance_ratio_file != "") {
      cat(sprintf("✓ Variance ratio file: %s\n", status_result$variance_ratio_file))
    }
  } else {
    cat("✗ Step 1 failed. Check status file for details.\n")
  }
  cat(sprintf("📄 Status report: %s.status.txt\n", opt$outputPrefix))
}

if (BLASctl_installed) {
  # Restore originally configured BLAS thread count
  blas_set_num_threads(original_num_threads)
}
