# =========================
# Convert Entrez IDs to Gene Symbols using BioMart
# =========================

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
suppressPackageStartupMessages(library(biomaRt))

# ---- User Inputs ----
input_file <- "/Users/jackie-mac/Desktop/VSCode/outputs/Filtered_RNAseq/NSCLC_filtered_RNAseq.csv"   # <-- Set your input file path
output_file <- "/Users/jackie-mac/Desktop/VSCode/outputs/Filtered_RNAseq/gensym_NSCLC_filtered_RNAseq.csv" # <-- Set your output file path

# ---- Read RNAseq matrix (genes x samples) ----
rna_mat <- read.csv(input_file, row.names = 1, check.names = FALSE)

# ---- Get Entrez IDs ----
entrez_ids <- rownames(rna_mat)

# ---- Connect to BioMart ----
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# ---- Query BioMart for gene symbols ----
bm <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = entrez_ids,
  mart = mart
)

# ---- Create mapping ----
id_map <- setNames(bm$hgnc_symbol, as.character(bm$entrezgene_id))

# ---- Replace Entrez IDs with gene symbols ----
gene_symbols <- id_map[as.character(entrez_ids)]
# If mapping fails, keep original Entrez ID
gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- entrez_ids[is.na(gene_symbols) | gene_symbols == ""]

# ---- Make gene symbols unique ----
gene_symbols_unique <- make.unique(gene_symbols)
rownames(rna_mat) <- gene_symbols_unique

# ---- Write to new CSV ----
write.csv(rna_mat, file = output_file, row.names = TRUE)