# ===============================================================================
# Kaplan-Meier Survival Curve Generator for Radiogenomic Signatures
# ===============================================================================
# 
# Purpose: Generates Kaplan-Meier survival curves comparing high vs low expression
#          groups for genomic signatures and radiomic features using median split
#          stratification to visualize survival differences.
#
# Description:
#   This script creates survival curve plots for radiogenomic features by dividing
#   patients into high and low expression groups based on median values. It generates
#   publication-ready Kaplan-Meier curves with log-rank test statistics to assess
#   the prognostic value of different signatures and features.
#
# Input Requirements:
#   1. Feature expression data: CSV with samples as rows, features as columns
#   2. Clinical survival data: Must include overall survival time and event status
#   3. Matching sample IDs between feature and clinical datasets
#
# Output:
#   - Kaplan-Meier survival curve plots (PDF/PNG format)
#   - Log-rank test p-values for survival differences
#   - Risk tables showing number at risk over time
#   - Summary statistics for high vs low groups
#
# Analysis Method:
#   - Median split stratification for high/low group assignment
#   - Kaplan-Meier survival estimation
#   - Log-rank test for statistical significance testing
#   - Customizable plotting parameters and aesthetics
#
# Usage:
#   1. Configure input file paths and parameters in the script
#   2. Run: Rscript survival_curves.R
#   3. Review generated survival plots and statistics
#
# Dependencies: survival, survminer, ggplot2
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

library(survival)
library(survminer)
library(data.table)
library(ggplot2)

# ---- USER INPUTS ----
# File lists (update paths as needed)
cox_files <- list(
  kegg = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/CCRCC/CCRCC_KEGG_cox_results_binary.csv",
  hallmark = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/CCRCC/CCRCC_HALLMARK_cox_results_binary.csv",
  biocarta = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/CCRCC/CCRCC_BIOCARTA_cox_results_binary.csv",
  reactome = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/CCRCC/CCRCC_REACTOME_cox_results_binary.csv"
)
correlation_files <- list(
  kegg = "/Users/jackie-mac/Desktop/VSCode/outputs/correlations/CCRCC/CCRCC_KEGG_correlative_analysis.csv",
  hallmark = "/Users/jackie-mac/Desktop/VSCode/outputs/correlations/CCRCC/CCRCC_HALLMARK_correlative_analysis.csv",
  biocarta = "/Users/jackie-mac/Desktop/VSCode/outputs/correlations/CCRCC/CCRCC_BIOCARTA_correlative_analysis.csv",
  reactome = "/Users/jackie-mac/Desktop/VSCode/outputs/correlations/CCRCC/CCRCC_REACTOME_correlative_analysis.csv"
)
clinical_file <- "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/CCRCC/CCRCC_harmonized_clinical.csv"
genomic_files <- list(
  kegg = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/CCRCC/CCRCC_harmonized_kegg_genomics.csv",
  hallmark = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/CCRCC/CCRCC_harmonized_hallmark_genomics.csv",
  biocarta = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/CCRCC/CCRCC_harmonized_biocarta_genomics.csv",
  reactome = "/Users/jackie-mac/Desktop/VSCode/outputs/clinical_associations/CCRCC/CCRCC_harmonized_reactome_genomics.csv"
)
output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/plots/survival_curves/CCRCC"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- 1. STANDARDIZE SAMPLE ID COLUMN NAMES ----
standardize_sampleid <- function(file, old_col) {
  dt <- fread(file)
  if (colnames(dt)[1] != "SampleID") {
    setnames(dt, old_col, "SampleID")
  }
  return(dt)
}
clinical <- standardize_sampleid(clinical_file, "cases.submitter_id_CCRCC")

for (g in names(genomic_files)) {
  dt <- fread(genomic_files[[g]])
  if (colnames(dt)[1] != "SampleID") setnames(dt, 1, "SampleID")
  assign(paste0(g, "_genomic"), dt)
}

# ---- 2. FILTER COX MODEL RESULTS BY FDR ----
filtered_features <- list()
cox_results <- list()  # Store full cox results for C-index extraction
for (g in names(cox_files)) {
  cox <- fread(cox_files[[g]])
  cox <- cox[!is.na(FDR) & FDR <= 0.05]
  filtered_features[[g]] <- cox$Signature
  cox_results[[g]] <- cox  # Store for later C-index extraction
}

# --- 2.5 FILTER GENOMIC FEATURES BY FDR-PASSING SIGNATURES ---
for (g in names(cox_files)) {
  # Get the cox file for this group and select top 10 signatures by FDR
  cox <- cox_results[[g]]  # Use stored cox results
  cox <- cox[order(FDR, decreasing = FALSE)]
  top_signatures <- head(cox$Signature, 10)
  genomic_mat <- get(paste0(g, "_genomic"))
  filtered_features[[g]] <- top_signatures
}
# ---- 3. CORRELATION FILTERING ----
final_features <- list()
for (g in names(correlation_files)) {
  corr <- fread(correlation_files[[g]])
  # Robust matching: ensure both columns are character, trimmed, and case-matched
  corr$GenomicFeature <- trimws(as.character(corr$GenomicFeature))
  filtered_set <- trimws(as.character(filtered_features[[g]]))
  # Only keep pairs with |SpearmanRho| > 0.4 and FDR-passing genomic features
  keep <- corr[abs(SpearmanRho) > 0.4 & GenomicFeature %in% filtered_set]
  # For each GenomicFeature, keep only the row with the highest absolute SpearmanRho
  keep[, absRho := abs(SpearmanRho)]
  setorder(keep, GenomicFeature, -absRho)
  keep <- keep[!duplicated(GenomicFeature)]
  keep[, absRho := NULL]
  final_features[[g]] <- keep
}

# ---- 4. GENERATE SURVIVAL CURVES FOR EACH GENOMIC FEATURE ----
for (g in names(final_features)) {
  if (nrow(final_features[[g]]) == 0) next
  genomic_mat <- get(paste0(g, "_genomic"))
  cox_data <- cox_results[[g]]  # Get cox results for this pathway
  
  for (i in seq_len(nrow(final_features[[g]]))) {
    feature <- final_features[[g]]$GenomicFeature[i]
    radiomic <- final_features[[g]]$RadiomicFeature[i]
    
    # Extract C-index from Cox results (already calculated for binary high/low groups)
    cox_c_index <- cox_data[Signature == feature, C_index]
    if (length(cox_c_index) == 0 || is.na(cox_c_index)) {
      cox_c_index <- "N/A"
    } else {
      cox_c_index <- round(as.numeric(cox_c_index), 3)
    }
    
    # Merge clinical and feature
    merged <- merge(clinical, genomic_mat[, .(SampleID, value = get(feature))], by = "SampleID")
    
    # Remove rows with missing values
    merged <- merged[!is.na(value) & !is.na(OS_days_CCRCC) & !is.na(OS_event_CCRCC)]
    
    if (nrow(merged) < 10) {
      cat("Warning: Too few samples with complete data for", feature, ". Skipping.\n")
      next
    }
    
    # Split by median (same as original Cox model)
    median_val <- median(merged$value, na.rm = TRUE)
    merged$group <- ifelse(merged$value > median_val, "High", "Low")
    
    # Create annotation text with C-index from Cox model
    c_index_text <- paste0("C-index = ", cox_c_index)
    
    # Survival curve
    surv_obj <- Surv(time = merged$OS_days_CCRCC, event = merged$OS_event_CCRCC)
    plot_title <- paste0("Survival Curve: High vs Low ", feature, " (", g, ")\nCorrelated with radiomic: ", radiomic)
    p <- ggsurvplot(
      survfit(surv_obj ~ group, data = merged),
      data = merged, pval = TRUE, risk.table = TRUE, conf.int = TRUE,
      title = plot_title,
      legend.title = feature, legend.labs = c("Low", "High"),
      palette = c("#377EB8", "#E41A1C"),
      ggtheme = theme_minimal(base_size = 18) +
        theme(panel.background = element_rect(fill = "white", color = NA),
              plot.background = element_rect(fill = "white", color = NA),
              legend.background = element_rect(fill = "white", color = NA),
              legend.key = element_rect(fill = "white", color = NA),
              plot.title = element_text(size = 9),  # Much smaller title font
              plot.subtitle = element_text(size = 9)) # Smaller subtitle font if present
    )
    
    # Add C-index annotation to the plot
    p$plot <- p$plot + 
      annotate("text", 
               x = Inf, 
               y = Inf, 
               label = c_index_text, 
               hjust = 1.05, 
               vjust = 2.5, 
               size = 5, 
               color = "black",
               fontface = "bold")
    
    # Clean feature and radiomic names for filename
    clean_feature <- gsub("[^A-Za-z0-9_]", "_", feature)
    clean_radiomic <- gsub("[^A-Za-z0-9_]", "_", radiomic)
    
    ggsave(file.path(output_dir, paste0("survival_curve_", g, "_", clean_feature, "_corr_", clean_radiomic, ".png")), 
           p$plot, width = 14, height = 7, dpi = 300)
  }
}
cat("Survival curves generated for all selected feature pairs.\n")
