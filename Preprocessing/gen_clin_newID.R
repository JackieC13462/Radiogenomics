# ===============================
# Clinical_filter_uniqueID.R
# -------------------------------
# Purpose:
#   Harmonizes clinical and multiple genomics pathway datasets by aligning sample IDs across files.
#   Ensures clinical and genomic datasets have matching sample IDs for downstream analysis.
#
# Inputs:
#   1. clinical_file: Path to the clinical data file (CSV).
#   2. hallmark_genomics_file: Path to Hallmark genomics matrix (CSV).
#   3. kegg_genomics_file: Path to Kegg genomics matrix (CSV).
#   4. reactome_genomics_file: Path to Reactome genomics matrix (CSV).
#   5. biocarta_genomics_file: Path to Biocarta genomics matrix (CSV).
#   6. output_prefix: Prefix for output files.
#   7. output_dir: Directory to write harmonized outputs.
#
# Outputs:
#   - Harmonized clinical and genomics pathway files (all with matching sample IDs).
#   - Console messages summarizing harmonization steps.
#
# Main Steps:
#   - Reads all input files.
#   - Finds common sample IDs between clinical and genomic datasets.
#   - Removes samples that don't exist in both clinical and genomic data.
#   - Outputs harmonized files for each data type.
# ===============================

suppressPackageStartupMessages(library(data.table))

# ---- SCRIPT DESCRIPTION ----
# This script harmonizes clinical and multiple genomics pathway datasets by ensuring consistent sample IDs across files.
# It performs the following steps:
#
# 1. Reads clinical and four genomics pathway files (Hallmark, Kegg, Reactome, Biocarta) provided as input.
# 2. Automatically identifies sample ID columns (looks for "cases.submitter_id" in clinical data and "SampleID" in genomic data).
# 3. Identifies common sample IDs between clinical and genomic datasets.
# 4. Removes samples from clinical data that don't exist in genomic datasets.
# 5. Removes samples from genomic datasets that don't exist in clinical data.
# 6. Writes harmonized clinical and all four genomics pathway datasets to separate CSV files.
#
# The output files ensure that clinical and genomic datasets have matching sample IDs for downstream analysis.

# ---- USER INPUTS ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript gen_clin_newID.R <clinical_file> <hallmark_genomics_file> <kegg_genomics_file> <reactome_genomics_file> <biocarta_genomics_file> <output_dir>")
}
clinical_file <- args[1]
hallmark_file <- args[2]
kegg_file <- args[3]
reactome_file <- args[4]
biocarta_file <- args[5]
output_dir <- args[6]
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- READ DATA ----
clinical <- fread(clinical_file, data.table = FALSE, quote = "\"", check.names = FALSE)
hallmark <- fread(hallmark_file, data.table = FALSE, quote = "\"", check.names = FALSE)
kegg <- fread(kegg_file, data.table = FALSE, quote = "\"", check.names = FALSE)
reactome <- fread(reactome_file, data.table = FALSE, quote = "\"", check.names = FALSE)
biocarta <- fread(biocarta_file, data.table = FALSE, quote = "\"", check.names = FALSE)

cat("Initial data dimensions:\n")
cat("Clinical:", nrow(clinical), "samples\n")
cat("Hallmark:", nrow(hallmark), "samples\n")
cat("KEGG:", nrow(kegg), "samples\n")
cat("Reactome:", nrow(reactome), "samples\n")
cat("BioCarta:", nrow(biocarta), "samples\n")

# ---- IDENTIFY SAMPLE ID COLUMNS ----
# For clinical data, look for cases.submitter_id column, fallback to first column
clinical_id_col <- if("cases.submitter_id" %in% colnames(clinical)) {
  "cases.submitter_id"
} else {
  colnames(clinical)[1]
}

# For genomic data, look for SampleID column, fallback to first column
genomic_id_col <- if("SampleID" %in% colnames(hallmark)) {
  "SampleID"
} else {
  colnames(hallmark)[1]
}

cat("Using sample ID columns:\n")
cat("Clinical:", clinical_id_col, "\n")
cat("Genomic datasets:", genomic_id_col, "\n")

# ---- FIND COMMON SAMPLE IDs BETWEEN CLINICAL AND GENOMIC DATASETS ----
clinical_ids <- clinical[[clinical_id_col]]
genomic_ids_list <- list(
  hallmark = hallmark[[genomic_id_col]],
  kegg = kegg[[genomic_id_col]],
  reactome = reactome[[genomic_id_col]],
  biocarta = biocarta[[genomic_id_col]]
)

# Find sample IDs that exist in ALL genomic datasets
common_genomic_ids <- Reduce(intersect, genomic_ids_list)
cat("Sample IDs common to all genomic datasets:", length(common_genomic_ids), "\n")

# Find sample IDs that exist in both clinical data and all genomic datasets
final_common_ids <- intersect(clinical_ids, common_genomic_ids)
cat("Sample IDs common to clinical and all genomic datasets:", length(final_common_ids), "\n")

# Report removed samples
clinical_only_ids <- setdiff(clinical_ids, final_common_ids)
genomic_only_ids <- setdiff(common_genomic_ids, final_common_ids)

cat("Samples in clinical but not in genomic data:", length(clinical_only_ids), "\n")
cat("Samples in genomic data but not in clinical:", length(genomic_only_ids), "\n")

# ---- FILTER ALL DATASETS TO COMMON SAMPLE IDs ----
# Filter clinical data
clinical <- clinical[clinical[[clinical_id_col]] %in% final_common_ids, , drop = FALSE]

# Filter genomic datasets
hallmark <- hallmark[hallmark[[genomic_id_col]] %in% final_common_ids, , drop = FALSE]
kegg <- kegg[kegg[[genomic_id_col]] %in% final_common_ids, , drop = FALSE]
reactome <- reactome[reactome[[genomic_id_col]] %in% final_common_ids, , drop = FALSE]
biocarta <- biocarta[biocarta[[genomic_id_col]] %in% final_common_ids, , drop = FALSE]

# ---- SUMMARY AFTER FILTERING ----
cat("\nFinal harmonized data dimensions:\n")
cat("Clinical:", nrow(clinical), "samples\n")
cat("Hallmark:", nrow(hallmark), "samples\n")
cat("KEGG:", nrow(kegg), "samples\n")
cat("Reactome:", nrow(reactome), "samples\n")
cat("BioCarta:", nrow(biocarta), "samples\n")

# ---- WRITE OUTPUT ----
# Extract cancer type from clinical file name for output file naming
clinical_basename <- basename(clinical_file)
cancer_type <- strsplit(clinical_basename, "_")[[1]][1]

# Write clinical dataset
fwrite(clinical, file.path(output_dir, paste0(cancer_type, "_harmonized_clinical.csv")), sep = ",", quote = TRUE, na = "NA")

# Write genomic datasets
fwrite(hallmark, file.path(output_dir, paste0(cancer_type, "_harmonized_hallmark_genomics.csv")), sep = ",", quote = TRUE, na = "NA")
fwrite(kegg, file.path(output_dir, paste0(cancer_type, "_harmonized_kegg_genomics.csv")), sep = ",", quote = TRUE, na = "NA")
fwrite(reactome, file.path(output_dir, paste0(cancer_type, "_harmonized_reactome_genomics.csv")), sep = ",", quote = TRUE, na = "NA")
fwrite(biocarta, file.path(output_dir, paste0(cancer_type, "_harmonized_biocarta_genomics.csv")), sep = ",", quote = TRUE, na = "NA")

cat("\nHarmonization completed successfully!\n")
cat("All datasets now have matching sample IDs for downstream analysis.\n")
cat("Harmonized files written to:", output_dir, "\n")
cat("Files created with prefix:", cancer_type, "\n")