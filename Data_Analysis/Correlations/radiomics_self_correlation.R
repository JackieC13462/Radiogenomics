# ===============================================================================
# Radiomic Features Self-Correlation Analysis
# ===============================================================================
# 
# Purpose: Computes Spearman correlation matrix among radiomic features to identify
#          highly correlated feature groups that may be redundant or represent
#          similar imaging characteristics.
#
# Description:
#   This script calculates pairwise Spearman correlations between all radiomic
#   features within a dataset. It filters features by allowed prefixes to focus
#   on clinically relevant original and transformed features while excluding
#   diagnostic metadata. Results help identify feature redundancy.
#
# Input Requirements:
#   1. Radiomic features file: CSV with samples as rows, radiomic features as columns
#   2. File must contain "SampleID" as first column header
#   3. Features should include original_, wavelet-, and transformation prefixes
#
# Output:
#   CSV file containing self-correlation matrix:
#   - Square matrix with radiomic features as both rows and columns
#   - Spearman correlation coefficients ranging from -1 to +1
#   - Diagonal values = 1 (perfect self-correlation)
#   - Used for identifying highly correlated feature clusters
#
# Analysis Method:
#   - Filters features by allowed prefixes (original_, wavelet-, square_, etc.)
#   - Uses Spearman rank correlation (robust to outliers)
#   - Computes pairwise complete observations
#   - Results used downstream for feature selection and redundancy removal
#
# Usage:
#   Rscript radiomics_self_correlation.R <radiomic_features_file> <output_correlation_file>
#
# Example:
#   Rscript radiomics_self_correlation.R HNSCC_radiomics.csv HNSCC_radiomics_correlation.csv
#
# Dependencies: Standard R libraries (stats, utils)
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

cat("Starting radiomics self-correlation script...\n")

# =========================
# Parse Command Line Arguments
# =========================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript radiomics_self_correlation.R <input_radiomics.csv> <output_corr.csv>")
}
input_radiomics <- args[1]
output_corr <- args[2]
cat("Input radiomics file:", input_radiomics, "\n")
cat("Output correlation file:", output_corr, "\n")
if (!dir.exists(dirname(output_corr))) {
  dir.create(dirname(output_corr), recursive = TRUE)
  cat("Created output directory:", dirname(output_corr), "\n")
}

# =========================
# Read Data
# =========================
cat("Reading radiomics data...\n")
radiomics_matrix <- fread(input_radiomics, data.table = FALSE)
cat("Radiomics data dimensions (including ID column):", dim(radiomics_matrix), "\n")
rownames(radiomics_matrix) <- radiomics_matrix[[1]]
radiomics_matrix <- radiomics_matrix[, -1, drop = FALSE]
cat("Radiomics matrix dimensions (samples x features):", dim(radiomics_matrix), "\n")

# =========================
# Filter features by allowed prefixes
# =========================
allowed_prefixes <- c("original_", "wavelet-", "square_", "squareroot_", "logarithm_", "exponential_", "gradient_")
prefix_pattern <- paste0("^(", paste(allowed_prefixes, collapse = "|"), ")")
filtered_feature_names <- grep(prefix_pattern, colnames(radiomics_matrix), value = TRUE)
radiomics_matrix <- radiomics_matrix[, filtered_feature_names, drop = FALSE]
cat("Filtered radiomics matrix dimensions (samples x allowed features):", dim(radiomics_matrix), "\n")

# =========================
# Ensure all values are numeric
# =========================
cat("Converting all values to numeric...\n")
radiomics_matrix[] <- lapply(radiomics_matrix, function(x) as.numeric(as.character(x)))

# =========================
# Compute Spearman Self-Correlation Matrix
# =========================
cat("Computing Spearman self-correlation matrix...\n")
if (ncol(radiomics_matrix) >= 2) {
  cor_mat <- cor(radiomics_matrix, method = "spearman", use = "pairwise.complete.obs")
  cat("Correlation matrix computed. Dimensions:", dim(cor_mat), "\n")
} else {
  stop("Not enough features (columns) for correlation analysis after filtering.")
}

# =========================
# Save Correlation Matrix
# =========================
cat("Saving correlation matrix to file...\n")
write.csv(
  cor_mat,
  file = output_corr,
  row.names = TRUE
)
cat("Radiomics self-correlation complete. Output saved at:", output_corr, "\n")
