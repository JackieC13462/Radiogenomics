# =========================
# Genomic Feature Self-Correlation (Spearman) by Pathway Source
# =========================

suppressPackageStartupMessages({
  library(data.table)
})

cat("Starting genomics self-correlation script...\n")

# =========================
# Parse Command Line Arguments
# =========================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript genomics_self_correlation.R <input_genomics.csv> <prefix> <output_dir>")
}
input_genomics <- args[1]
prefix <- args[2]
output_dir <- args[3]
cat("Input genomics file:", input_genomics, "\n")
cat("Prefix:", prefix, "\n")
cat("Output directory:", output_dir, "\n")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# =========================
# Read Data
# =========================
cat("Reading genomics data...\n")
genomic_matrix <- fread(input_genomics, data.table = FALSE)
cat("Genomics data dimensions (including ID column):", dim(genomic_matrix), "\n")
rownames(genomic_matrix) <- genomic_matrix[[1]]
genomic_matrix <- genomic_matrix[, -1, drop = FALSE]
cat("Genomics matrix dimensions (pathways x samples):", dim(genomic_matrix), "\n")

# =========================
# Pathway source patterns
# =========================
sources <- c("KEGG", "HALLMARKS", "REACTOME", "BIOCARTA")

for (source in sources) {
  cat("\nProcessing pathway group:", source, "\n")
  # Identify pathway rows for this source
  source_rows <- grepl(source, rownames(genomic_matrix), ignore.case = TRUE)
  cat("Number of pathways found for", source, ":", sum(source_rows), "\n")
  if (!any(source_rows)) {
    cat("No pathways found for", source, "\n")
    next
  }
  # Subset matrix for this source
  sub_mat <- genomic_matrix[source_rows, , drop = FALSE]
  cat("Subset matrix dimensions (pathways x samples):", dim(sub_mat), "\n")
  # Transpose: pathways x samples -> samples x pathways
  sub_mat <- t(sub_mat)
  sub_mat <- as.data.frame(sub_mat)
  cat("Transposed matrix dimensions (samples x pathways):", dim(sub_mat), "\n")
  # Ensure all values are numeric
  cat("Converting all values to numeric...\n")
  sub_mat[] <- lapply(sub_mat, function(x) as.numeric(as.character(x)))
  # Compute correlation if enough features
  if (ncol(sub_mat) >= 2) {
    cat("Computing Spearman self-correlation matrix for", source, "...\n")
    cor_mat <- cor(sub_mat, method = "spearman", use = "pairwise.complete.obs")
    out_file <- file.path(output_dir, paste0(prefix, "_", source, "_self_correlation.csv"))
    write.csv(cor_mat, file = out_file, row.names = TRUE)
    cat("Self-correlation matrix saved for", source, "at", out_file, "\n")
  } else {
    cat("Not enough features (columns) for correlation analysis in", source, "\n")
  }
}
cat("Genomics self-correlation script complete.\n")
