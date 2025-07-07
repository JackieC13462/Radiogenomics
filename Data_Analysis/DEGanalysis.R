suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

# ---- USER INPUTS ----
rna_files <- list(
  BRCA = "/path/to/BRCA_counts.tsv",
  LIHC = "/path/to/LIHC_counts.tsv"
  # Add more cancer types and their count files
)
metadata_files <- list(
  BRCA = "/path/to/BRCA_metadata.tsv",
  LIHC = "/path/to/LIHC_metadata.tsv"
  # Add more cancer types and their metadata files
)
group_column <- "condition"  # Column in metadata indicating group (e.g., "tumor" vs "normal")
top_k <- 20                  # Number of top DEGs to report per cancer type
output_dir <- "/path/to/DEG_results/"

dir.create(output_dir, showWarnings = FALSE)

# --- Combine all counts ---
all_counts <- list()
all_metadata <- list()
for (cancer in names(rna_files)) {
  counts <- fread(rna_files[[cancer]], data.table = FALSE)
  rownames(counts) <- counts[[1]]
  counts <- counts[,-1]
  all_counts[[cancer]] <- counts
  # Build metadata for this cancer
  meta <- data.frame(
    Sample = colnames(counts),
    TumourType = cancer,
    stringsAsFactors = FALSE
  )
  all_metadata[[cancer]] <- meta
}
# Find common genes
common_genes <- Reduce(intersect, lapply(all_counts, rownames))
all_counts <- lapply(all_counts, function(x) x[common_genes, , drop=FALSE])
counts_mat <- do.call(cbind, all_counts)
metadata <- do.call(rbind, all_metadata)
rownames(metadata) <- metadata$Sample

# --- Add BRCA_status column ---
metadata$BRCA_status <- ifelse(metadata$TumourType == "BRCA", 1, 0)

# --- Ensure order matches ---
counts_mat <- counts_mat[, metadata$Sample]

# --- DESeq2 analysis ---
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_mat)),
  colData = metadata,
  design = ~ BRCA_status + TumourType
)
dds <- DESeq(dds)
res <- results(dds, name = "BRCA_status")
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), "BRCA_vs_all_DEGs.csv")