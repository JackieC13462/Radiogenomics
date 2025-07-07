suppressPackageStartupMessages({
  library(data.table)
  library(uwot)
  library(ggplot2)
  library(plotly)
})

# ---- USER INPUTS ----
rna_files <- list(
  BRCA = "/Users/jackie-mac/Desktop/VSCode/Outputs/Filtered_RNAseq/Protein_coding/MissingRemoved/TPM_Normalized/BRCA_filtered_TPM.csv",
  KIRC = "/Users/jackie-mac/Desktop/VSCode/Outputs/Filtered_RNAseq/Protein_coding/MissingRemoved/TPM_Normalized/KIRC_filtered_TPM.csv",
  LGG = "/Users/jackie-mac/Desktop/VSCode/Outputs/Filtered_RNAseq/Protein_coding/MissingRemoved/TPM_Normalized/LGG_filtered_TPM.csv",
  GBM = "/Users/jackie-mac/Desktop/VSCode/Outputs/Filtered_RNAseq/Protein_coding/MissingRemoved/TPM_Normalized/GBM_filtered_TPM.csv",
  CCRCC = "/Users/jackie-mac/Desktop/VSCode/Outputs/Filtered_RNAseq/Protein_coding/MissingRemoved/TPM_Normalized/CCRCC_filtered_TPM.csv",
  HNSCC = "/Users/jackie-mac/Desktop/VSCode/Outputs/Filtered_RNAseq/Protein_coding/MissingRemoved/TPM_Normalized/HNSCC_filtered_TPM.csv",
  PDA = "/Users/jackie-mac/Desktop/VSCode/Outputs/Filtered_RNAseq/Protein_coding/MissingRemoved/TPM_Normalized/PDA_filtered_TPM.csv"
  #NSCLC = "/Users/jackie-mac/Desktop/VSCode/Outputs/Filtered_RNAseq/Protein_coding/MissingRemoved/TPM_Normalized/NSCLC_filtered_TPM.csv"
)
output_umap_plot <- "/Users/jackie-mac/Desktop/VSCode/Outputs/umap_cancertype_plot.png"

# ---- CONSORTIUM LABELS ----
tcga_cancers <- c("BRCA", "GBM", "KIRC", "LGG")
cptac_cancers <- c("CCRCC", "HNSCC", "PDA")
#nsclc_cancer <- "NSCLC" # NSCLC will be its own shape

all_expr <- list()
all_labels <- c()
all_shapes <- c()
for (cancer in names(rna_files)) {
  expr <- fread(rna_files[[cancer]], data.table = FALSE)
  rownames(expr) <- expr[[1]]
  expr <- expr[,-1]
  all_expr[[cancer]] <- expr
  n_samples <- ncol(expr)
  all_labels <- c(all_labels, rep(cancer, n_samples))
  if (cancer == "NSCLC") {
    all_shapes <- c(all_shapes, rep("NSCLC", n_samples))
  } else if (cancer %in% cptac_cancers) {
    all_shapes <- c(all_shapes, rep("CPTAC", n_samples))
  } else {
    all_shapes <- c(all_shapes, rep("TCGA", n_samples))
  }
}
# Find common genes across all datasets
common_genes <- Reduce(intersect, lapply(all_expr, rownames))
# Subset and concatenate
all_expr <- lapply(all_expr, function(x) x[common_genes, , drop = FALSE])
expr_mat <- do.call(cbind, all_expr)
expr_mat <- t(expr_mat) # samples x genes

# ---- RUN UMAP ----
set.seed(123)
umap_res <- umap(expr_mat, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")

umap_df <- data.frame(
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2],
  CancerType = all_labels,
  Consortium = factor(all_shapes, levels = c("CPTAC", "TCGA", "NSCLC"))
)

# ---- PLOT ----
# Use more distinct shapes: CPTAC = X (shape 4), TCGA = triangle (17), NSCLC = circle (16)
shape_map <- c("CPTAC" = 4, "TCGA" = 17, "NSCLC" = 16) # X, triangle, circle

p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = CancerType, shape = Consortium, label = CancerType)) +
  geom_point(alpha = 0.85, size = 2.5, stroke = 1, fill = "white") + # stroke adds border, fill for open shapes
  scale_shape_manual(values = shape_map) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  ) +
  labs(title = "UMAP of RNAseq Data by Cancer Type and Consortium")

# Save static plot as PNG
ggsave(output_umap_plot, p, width = 7, height = 5, bg = "white")

# ---- PLOT (interactive) ----
p_interactive <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = CancerType, shape = Consortium, text = CancerType)) +
  geom_point(alpha = 0.85, size = 2.5, stroke = 1, fill = "white") +
  scale_shape_manual(values = shape_map) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  ) +
  labs(title = "UMAP of RNAseq Data by Cancer Type and Consortium")

# Convert ggplot to interactive plotly object
p_interactive <- ggplotly(p_interactive, tooltip = c("text", "shape"))

# Save as HTML (interactive)
htmlwidgets::saveWidget(p_interactive, file = sub("\\.png$", ".html", output_umap_plot))

cat("UMAP plot saved to:", output_umap_plot, "and interactive HTML.\n")