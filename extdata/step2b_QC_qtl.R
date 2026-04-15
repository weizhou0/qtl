#!/usr/bin/env -S pixi run --manifest-path ../pixi.toml Rscript

options(stringsAsFactors = F)
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# =========================
# CLI options
# =========================
option_list <- list(
  make_option("--SAIGEstep1OutputPrefix",
              type = "character", default = "", help = "Step1 output prefix (no suffix)"),
  make_option("--SAIGEstep2OutputPrefix",
              type = "character", default = "", help = "Step2 output prefix (no suffix)"),
  make_option("--gene",
              type = "character", default = "NA", help = "Gene name (for reporting only)"),
  make_option("--outPrefix",
              type = "character", default = "step2b_qc", help = "Output prefix for QC summary"),
  make_option("--isPostStep2QC",
              type = "logical", default = FALSE, help = "Whether to remove step1 and step2 outputs if they fail QC standards"),
  make_option("--phiLower", type = "double", default = 0, help = "Lower bound for phi dispersion parameter"),
  make_option("--phiUpper", type = "double", default = 1.5, help = "Upper bound for phi dispersion parameter"),
  make_option("--lambdaLower", type = "double", default = 0.1, help = "Lower bound for GC lambda"),
  make_option("--lambdaUpper", type = "double", default = 1.5, help = "Upper bound for GC lambda")
)
opt <- parse_args(OptionParser(option_list = option_list))

# =========================
# Helper: compute phi
# =========================
get_phi_from_rda <- function(rda_file) {
  if (!file.exists(rda_file)) return(NA_real_)
  e <- new.env(parent = emptyenv())
  load(rda_file, envir = e)
  
  if (!exists("modglmm", envir = e)) return(NA_real_)
  modglmm <- e$modglmm
  Y  <- modglmm$y
  mu <- modglmm$fitted.values
  if (is.null(Y) || is.null(mu)) return(NA_real_)
  p <- length(modglmm$obj.noK$S_a)
  n <- length(Y)
  df <- n - p
  if (df <= 0) return(NA_real_)
  pearson <- (Y - mu)^2 / pmax(mu, 1e-8)
  phi <- sum(pearson, na.rm = TRUE) / df
  phi
}

# =========================
# Helper: compute GCLambda
# =========================
get_lambda_from_step2 <- function(file) {
  if (!file.exists(file)) return(NA_real_)
  dt <- tryCatch(fread(file), error = function(e) NULL)
  if (is.null(dt)) return(NA_real_)
  
  if (!"p.value" %in% names(dt)) return(NA_real_)
  pvals <- suppressWarnings(as.numeric(dt[["p.value"]]))
  pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]

  if (length(pvals) == 0) return(NA_real_)
  chisq_obs <- qchisq(1 - pvals, df = 1)
  obs_med <- median(chisq_obs, na.rm = TRUE)
  exp_med <- qchisq(0.9, df = 1)
  lambda <- obs_med / exp_med
  lambda
}

# =========================
# Output file paths
# =========================
step1_prefix <- opt$SAIGEstep1OutputPrefix
step2_prefix <- opt$SAIGEstep2OutputPrefix
rda_file <- paste0(step1_prefix, ".rda")
step2_file <- step2_prefix   # (main file, no suffix)

# =========================
# Compute metrics
# =========================
phi <- get_phi_from_rda(rda_file)
lambda <- get_lambda_from_step2(step2_file)

# =========================
# QC logic (flexible bounds)
# =========================
phi_outside <- !is.na(phi) && (phi < opt$phiLower || phi > opt$phiUpper)
lambda_outside <- !is.na(lambda) && (lambda < opt$lambdaLower || lambda > opt$lambdaUpper)
lambda_bad <- !is.na(lambda) && (lambda == 0 || is.infinite(lambda))
fail_both <- (phi_outside && lambda_outside) || lambda_bad
# for reporting
phi_fail <- phi_outside
lambda_fail <- lambda_outside || lambda_bad

# =========================
# Removal logic
# =========================
action <- "kept"

if (opt$isPostStep2QC && fail_both) {
  step1_files <- c(paste0(step1_prefix, ".rda"), paste0(step1_prefix, ".status.txt"), paste0(step1_prefix, ".varianceRatio.txt"))
  step2_files <- c(step2_prefix, paste0(step2_prefix, ".index"))
  files_to_remove <- c(step1_files, step2_files)
  files_to_remove <- files_to_remove[file.exists(files_to_remove)]
  if (length(files_to_remove) > 0) {
    file.remove(files_to_remove)
    action <- "removed"
  } else {
    action <- "no_files_found"
  }
}

# =========================
# Output summary
# =========================
res <- data.table(
  gene = opt$gene, phi = phi, lambda_0.9 = lambda,
  phi_fail = phi_fail, lambda_fail = lambda_fail, fail_both = fail_both,
  action = action,
  step1_prefix = step1_prefix, step2_prefix = step2_prefix,
  phiLower = opt$phiLower, phiUpper = opt$phiUpper,
  lambdaLower = opt$lambdaLower, lambdaUpper = opt$lambdaUpper
)
outfile <- paste0(opt$outPrefix, ".", opt$gene, ".qc.tsv")
fwrite(res, outfile, sep = "\t", quote = FALSE, na = "NA")

# =========================
# Logging
# =========================
cat("Gene:", opt$gene, "\n")
cat("Phi:", phi, "\n")
cat("Lambda_0.9:", lambda, "\n")
cat("Fail both:", fail_both, "\n")
cat("Action:", action, "\n")
cat("Output:", outfile, "\n")