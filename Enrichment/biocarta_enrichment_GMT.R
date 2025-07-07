suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(GSVA))

perform_biocarta_samplewise_gsva <- function(input_file, gmt_file, output_file) {
  cat("Loading dataset from:", input_file, "\n")
  expression_data <- fread(input_file, data.table = FALSE, check.names = FALSE)
  
  # Set rownames from the first column and remove it
  rownames(expression_data) <- expression_data[[1]]
  expression_data[[1]] <- NULL
  
  # Read BioCarta gene sets from GMT file
  gmt <- getGmt(gmt_file)
  gene_set_list <- geneIds(gmt)
  # pathway_names <- names(gene_set_list) # Not needed, will update after filtering
  
  # ---- Filter gene sets: keep only those with >=80% genes present in RNAseq data ----
  rnaseq_genes <- rownames(expression_data)
  filtered_gene_set_list <- lapply(gene_set_list, function(genes) {
    intersect(genes, rnaseq_genes)
  })
  filtered_gene_set_list <- filtered_gene_set_list[
    sapply(seq_along(filtered_gene_set_list), function(i) {
      length(filtered_gene_set_list[[i]]) / length(gene_set_list[[i]]) >= 0.8
    })
  ]

  cat("Number of gene sets after filtering for >=80% gene coverage:", length(filtered_gene_set_list), "\n")
  
  # Prepare expression matrix for GSVA (genes as rows, samples as columns)
  expr_mat <- as.matrix(expression_data)
  gsvaPar <- gsvaParam(expr_mat, filtered_gene_set_list)
  print(gsvaPar)
  cat("Running GSVA using gsvaParam object...\n")
  gsva_results <- gsva(gsvaPar, verbose = TRUE)
  results <- t(gsva_results)
  colnames(results) <- colnames(gsva_results)
  rownames(results) <- colnames(expression_data)
  write.csv(results, output_file, row.names = TRUE, quote = FALSE)
  cat("Sample-wise GSVA results saved to:", output_file, "\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript biocarta_enrichment_GMT.R <input_file> <gmt_file> [output_file]")
}
input_file <- args[1]
gmt_file <- args[2]
output_file <- ifelse(length(args) >= 3, args[3], "biocarta_gsva_GMT_results.csv")

cat("Starting sample-wise BioCarta GSVA analysis...\n")
perform_biocarta_samplewise_gsva(input_file, gmt_file, output_file)
cat("Process completed.\n")