#!/usr/bin/env Rscript

# SAIGEQTL Output Comparison Script

# Handle library path from environment
library_path <- Sys.getenv("SAIGEQTL_LIBRARY_PATH", "")
if (library_path != "") {
  .libPaths(c(library_path, .libPaths()))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: Rscript compare_outputs.R <plink_file> <vcf_file> <bgen_file>")
}

plink_file <- args[1]
vcf_file <- args[2]
bgen_file <- args[3]

cat("Reading output files...\n")

# Read the output files
plink_data <- tryCatch({
    read.table(plink_file, header = TRUE, stringsAsFactors = FALSE)
}, error = function(e) {
    cat("Error reading PLINK file:", e$message, "\n")
    return(NULL)
})

vcf_data <- tryCatch({
    read.table(vcf_file, header = TRUE, stringsAsFactors = FALSE)
}, error = function(e) {
    cat("Error reading VCF file:", e$message, "\n")
    return(NULL)
})

bgen_data <- tryCatch({
    read.table(bgen_file, header = TRUE, stringsAsFactors = FALSE)
}, error = function(e) {
    cat("Error reading BGEN file:", e$message, "\n")
    return(NULL)
})

if (is.null(plink_data) || is.null(vcf_data) || is.null(bgen_data)) {
    stop("Failed to read one or more input files")
}

cat("Comparing outputs...\n\n")

# Basic statistics
cat("=== BASIC STATISTICS ===\n")
cat(sprintf("PLINK: %d variants\n", nrow(plink_data)))
cat(sprintf("VCF:   %d variants\n", nrow(vcf_data)))
cat(sprintf("BGEN:  %d variants\n", nrow(bgen_data)))
cat("\n")

# Column names comparison
cat("=== COLUMN COMPARISON ===\n")
cat("PLINK columns:", paste(colnames(plink_data), collapse = ", "), "\n")
cat("VCF columns:  ", paste(colnames(vcf_data), collapse = ", "), "\n")
cat("BGEN columns: ", paste(colnames(bgen_data), collapse = ", "), "\n")
cat("\n")

# Find common variants (assuming CHR:POS or similar identifier)
# Try different potential ID columns
id_cols <- c("MarkerID", "SNPID", "CHR_POS", "rsID", "ID")
id_col <- NULL

for (col in id_cols) {
    if (col %in% colnames(plink_data) && col %in% colnames(vcf_data) && col %in% colnames(bgen_data)) {
        id_col <- col
        break
    }
}

if (is.null(id_col)) {
    # Try to create CHR:POS identifier if CHR and POS columns exist
    if ("CHR" %in% colnames(plink_data) && "POS" %in% colnames(plink_data)) {
        plink_data$ID <- paste(plink_data$CHR, plink_data$POS, sep = ":")
        vcf_data$ID <- paste(vcf_data$CHR, vcf_data$POS, sep = ":")
        bgen_data$ID <- paste(bgen_data$CHR, bgen_data$POS, sep = ":")
        id_col <- "ID"
    } else {
        cat("Warning: No common identifier column found. Using row order for comparison.\n")
        cat("This may not be accurate if variants are in different orders.\n\n")
        
        # Show first few rows for manual inspection
        cat("=== FIRST 5 ROWS COMPARISON ===\n")
        cat("PLINK (first 5 rows):\n")
        print(head(plink_data, 5))
        cat("\nVCF (first 5 rows):\n")
        print(head(vcf_data, 5))
        cat("\nBGEN (first 5 rows):\n")
        print(head(bgen_data, 5))
        
        return()
    }
}

# Compare common variants
cat(sprintf("Using '%s' as identifier column\n", id_col))

common_plink_vcf <- intersect(plink_data[[id_col]], vcf_data[[id_col]])
common_plink_bgen <- intersect(plink_data[[id_col]], bgen_data[[id_col]])
common_vcf_bgen <- intersect(vcf_data[[id_col]], bgen_data[[id_col]])
common_all <- intersect(common_plink_vcf, common_plink_bgen)

cat(sprintf("Common variants PLINK-VCF:  %d\n", length(common_plink_vcf)))
cat(sprintf("Common variants PLINK-BGEN: %d\n", length(common_plink_bgen)))
cat(sprintf("Common variants VCF-BGEN:   %d\n", length(common_vcf_bgen)))
cat(sprintf("Common variants (all):      %d\n", length(common_all)))
cat("\n")

if (length(common_all) == 0) {
    cat("No common variants found across all formats!\n")
    return()
}

# Subset to common variants and compare key statistics
plink_common <- plink_data[plink_data[[id_col]] %in% common_all, ]
vcf_common <- vcf_data[vcf_data[[id_col]] %in% common_all, ]
bgen_common <- bgen_data[bgen_data[[id_col]] %in% common_all, ]

# Sort by ID for proper comparison
plink_common <- plink_common[order(plink_common[[id_col]]), ]
vcf_common <- vcf_common[order(vcf_common[[id_col]]), ]
bgen_common <- bgen_common[order(bgen_common[[id_col]]), ]

# Compare p-values if available
p_cols <- c("p.value", "Pvalue", "P", "pval")
p_col <- NULL
for (col in p_cols) {
    if (col %in% colnames(plink_common)) {
        p_col <- col
        break
    }
}

if (!is.null(p_col)) {
    cat("=== P-VALUE COMPARISON ===\n")
    
    plink_p <- plink_common[[p_col]]
    vcf_p <- vcf_common[[p_col]]
    bgen_p <- bgen_common[[p_col]]
    
    # Correlation between p-values
    cor_plink_vcf <- cor(plink_p, vcf_p, use = "complete.obs")
    cor_plink_bgen <- cor(plink_p, bgen_p, use = "complete.obs")
    cor_vcf_bgen <- cor(vcf_p, bgen_p, use = "complete.obs")
    
    cat(sprintf("P-value correlation PLINK-VCF:  %.6f\n", cor_plink_vcf))
    cat(sprintf("P-value correlation PLINK-BGEN: %.6f\n", cor_plink_bgen))
    cat(sprintf("P-value correlation VCF-BGEN:   %.6f\n", cor_vcf_bgen))
    
    # Max differences
    max_diff_plink_vcf <- max(abs(plink_p - vcf_p), na.rm = TRUE)
    max_diff_plink_bgen <- max(abs(plink_p - bgen_p), na.rm = TRUE)
    max_diff_vcf_bgen <- max(abs(vcf_p - bgen_p), na.rm = TRUE)
    
    cat(sprintf("Max p-value diff PLINK-VCF:  %.2e\n", max_diff_plink_vcf))
    cat(sprintf("Max p-value diff PLINK-BGEN: %.2e\n", max_diff_plink_bgen))
    cat(sprintf("Max p-value diff VCF-BGEN:   %.2e\n", max_diff_vcf_bgen))
    cat("\n")
}

# Compare effect sizes if available
beta_cols <- c("BETA", "beta", "Beta", "Effect")
beta_col <- NULL
for (col in beta_cols) {
    if (col %in% colnames(plink_common)) {
        beta_col <- col
        break
    }
}

if (!is.null(beta_col)) {
    cat("=== EFFECT SIZE COMPARISON ===\n")
    
    plink_beta <- plink_common[[beta_col]]
    vcf_beta <- vcf_common[[beta_col]]
    bgen_beta <- bgen_common[[beta_col]]
    
    # Correlation between effect sizes
    cor_plink_vcf <- cor(plink_beta, vcf_beta, use = "complete.obs")
    cor_plink_bgen <- cor(plink_beta, bgen_beta, use = "complete.obs")
    cor_vcf_bgen <- cor(vcf_beta, bgen_beta, use = "complete.obs")
    
    cat(sprintf("Effect size correlation PLINK-VCF:  %.6f\n", cor_plink_vcf))
    cat(sprintf("Effect size correlation PLINK-BGEN: %.6f\n", cor_plink_bgen))
    cat(sprintf("Effect size correlation VCF-BGEN:   %.6f\n", cor_vcf_bgen))
    cat("\n")
}

cat("=== COMPARISON SUMMARY ===\n")
if (!is.null(p_col) && all(c(cor_plink_vcf, cor_plink_bgen, cor_vcf_bgen) > 0.999, na.rm = TRUE)) {
    cat("✓ EXCELLENT: All formats produce highly consistent results (r > 0.999)\n")
} else if (!is.null(p_col) && all(c(cor_plink_vcf, cor_plink_bgen, cor_vcf_bgen) > 0.99, na.rm = TRUE)) {
    cat("✓ GOOD: All formats produce consistent results (r > 0.99)\n")
} else if (!is.null(p_col)) {
    cat("⚠ WARNING: Some differences detected between formats\n")
} else {
    cat("? Unable to assess consistency - no p-value column found\n")
}

cat("Analysis complete.\n")
