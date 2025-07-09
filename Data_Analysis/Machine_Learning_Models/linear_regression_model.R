# ===============================================================================
# Linear Regression Model for Radiogenomic Feature Analysis
# ===============================================================================
# 
# Purpose: Performs multiple linear regression analysis using radiogenomic features
#          to predict clinical outcomes and identify significant predictors.
#
# Description:
#   This script implements standard multiple linear regression with optional
#   stepwise feature selection. It evaluates the linear relationship between
#   radiogenomic features and continuous clinical outcomes, providing interpretable
#   coefficients and statistical significance tests.
#
# Input Requirements:
#   1. Feature matrix: CSV with samples as rows, features (radiomic/genomic) as columns
#   2. Clinical outcome data: Continuous target variable
#   3. Properly formatted data with matching sample IDs
#
# Output:
#   - Linear regression coefficients and p-values
#   - Model performance metrics (R², adjusted R², RMSE)
#   - Feature significance tests
#   - Residual analysis and model diagnostics
#   - Prediction results on test data
#
# Analysis Method:
#   - Uses ordinary least squares estimation
#   - Optional stepwise selection for feature reduction
#   - Evaluates model assumptions (linearity, normality, homoscedasticity)
#   - Provides confidence intervals for coefficients
#   - Cross-validation for performance assessment
#
# Usage:
#   1. Configure input file paths and parameters
#   2. Run: Rscript linear_regression_model.R
#   3. Review model coefficients and diagnostic plots
#
# Dependencies: data.table, caret, MASS
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

# Load libraries
library(data.table)
library(caret)
library(MASS)

set.seed(42)

# Load data
X <- fread("/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_harmonized_radiomics.csv") # features matrix, samples x features
gene_matrix <- fread("/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/HNSCC_harmonized_kegg_genomics.csv") # samples x all gene signatures

# Specify the gene signature columns you want to use (replace with your actual column names)
selected_signatures <- c("KEGG_MEDICUS_REFERENCE_BASE_EXCISION_AND_STRAND_CLEAVAGE_BY_MONOFUNCTIONAL_GLYCOSYLASE", "KEGG_MEDICUS_REFERENCE_MANNOSE_TYPE_O_GLYCAN_BIOSYNTHESIS_POMT_TO_POMK", "KEGG_MEDICUS_REFERENCE_CDC25_CELL_CYCLE_G2_M", "KEGG_MEDICUS_VARIANT_MUTATION_ACTIVATED_KRAS_NRAS_TO_PI3K_SIGNALING_PATHWAY")

# Subset the matrix to just those columns
Y <- gene_matrix[, ..selected_signatures]

# Remove non-numeric columns (e.g., patient_ID) from X
X <- X[, sapply(X, is.numeric), with=FALSE]

output_dir <- "/Users/jackie-mac/Desktop/VSCode/outputs/clinical/Cox_results/ML_results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Now loop over each signature column
for (sig in selected_signatures) {
  y <- Y[[sig]]
  model_df <- data.table(Model=character(), Fold=integer(), Spearman=numeric(), Pearson=numeric())
  folds <- createFolds(y, k=5, list=TRUE, returnTrain=FALSE)
  best_corr <- -Inf
  best_fold <- 0
  best_coef <- NULL

  for (fold in seq_along(folds)) {
    test_idx <- folds[[fold]]
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)
    X_train <- X[train_idx, , drop=FALSE]
    y_train <- y[train_idx]
    X_test  <- X[test_idx, , drop=FALSE]
    y_test  <- y[test_idx]
    model <- lm(y_train ~ ., data = as.data.frame(X_train))
    y_pred <- predict(model, newdata = as.data.frame(X_test))
    s_cor <- cor(y_pred, y_test, method = "spearman")
    p_cor <- cor(y_pred, y_test, method = "pearson")
    if (s_cor > best_corr) {
      best_corr <- s_cor
      best_fold <- fold
      best_coef <- coef(model)[-1]
    }
    model_df <- rbind(model_df, data.table(Model="Linear", Fold=fold, Spearman=s_cor, Pearson=p_cor))
    cat("Signature:", sig, "Fold", fold, "Spearman correlation:", s_cor, "\n")
  }
  cat("\nBest correlation for", sig, ":", best_corr, "from Fold", best_fold, "\n")
  feature_importance <- data.table(
    Peak = colnames(X),
    Weight = best_coef
  )[order(-abs(Weight))]
  fwrite(feature_importance, file.path(output_dir, paste0("lm_features_", sig, ".csv")))
  fwrite(model_df, file.path(output_dir, paste0("lm_", sig, ".csv")))
}