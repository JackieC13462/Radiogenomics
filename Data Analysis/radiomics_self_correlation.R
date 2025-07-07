# =========================
# Radiomics Feature Self-Correlation (Spearman)
# =========================

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
