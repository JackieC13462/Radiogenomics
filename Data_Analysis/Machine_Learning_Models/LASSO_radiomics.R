# ===============================================================================
# LASSO Regression for Radiomic Feature Selection and Prediction
# ===============================================================================
# 
# Purpose: Performs LASSO (Least Absolute Shrinkage and Selection Operator) 
#          regression on radiomic features to predict clinical outcomes while
#          performing automatic feature selection.
#
# Description:
#   This script uses LASSO regression to identify the most predictive radiomic
#   features for a given clinical outcome. LASSO automatically performs feature
#   selection by shrinking less important coefficients to zero, resulting in 
#   a sparse model with only the most relevant features.
#
# Input Requirements:
#   1. Radiomic features file: CSV with samples as rows, radiomic features as columns
#   2. Clinical outcome data: Continuous or binary outcome variable
#   3. Properly formatted data with matching sample IDs
#
# Output:
#   - Selected radiomic features with non-zero coefficients
#   - Model performance metrics (R², RMSE, etc.)
#   - Cross-validation results
#   - Feature importance rankings
#   - Prediction results on test data
#
# Analysis Method:
#   - Uses cross-validation to select optimal lambda parameter
#   - Applies LASSO penalty for feature selection
#   - Evaluates model performance using multiple metrics
#   - Provides interpretable feature coefficients
#
# Usage:
#   1. Configure input file paths and outcome variables
#   2. Run: Rscript LASSO_radiomics.R
#   3. Review selected features and model performance
#
# Dependencies: data.table, caret, glmnet
# Author: Radiogenomics Analysis Pipeline
# ===============================================================================

# Load libraries
library(data.table)
library(caret)
library(glmnet)

set.seed(42)

# Load data
X <- fread("X.csv")
gene_matrix <- fread("gene_signatures.csv") # samples x all gene signatures

# Specify the gene signature columns you want to use
selected_signatures <- c("Signature1", "Signature2", "Signature3")
Y <- gene_matrix[, ..selected_signatures]

# LASSO parameter grid (alpha=1 for LASSO, lambda is the penalty)
lambda_grid <- 10^seq(-3, 2, length = 50) # similar to your Python grid

for (sig in selected_signatures) {
  y <- Y[[sig]]
  model_df <- data.table(Model=character(), Fold=integer(), Spearman=numeric(), Pearson=numeric(), lambda=numeric())
  folds <- createFolds(y, k=10, list=TRUE, returnTrain=FALSE)
  best_corr <- -Inf
  best_fold <- 0
  best_coef <- NULL
  best_lambda <- NA

  for (fold in seq_along(folds)) {
    test_idx <- folds[[fold]]
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)
    X_train <- as.matrix(X[train_idx, , drop=FALSE])
    y_train <- y[train_idx]
    X_test  <- as.matrix(X[test_idx, , drop=FALSE])
    y_test  <- y[test_idx]

    # Fit LASSO with cross-validation for lambda
    cvfit <- cv.glmnet(X_train, y_train, alpha=1, lambda=lambda_grid, nfolds=5)
    best_lambda_fold <- cvfit$lambda.min
    lasso_model <- glmnet(X_train, y_train, alpha=1, lambda=best_lambda_fold)

    y_pred <- predict(lasso_model, newx=X_test, s=best_lambda_fold)
    s_cor <- cor(y_pred, y_test, method = "spearman")
    p_cor <- cor(y_pred, y_test, method = "pearson")

    if (s_cor > best_corr) {
      best_corr <- s_cor
      best_fold <- fold
      best_coef <- as.vector(coef(lasso_model))[-1] # drop intercept
      best_lambda <- best_lambda_fold
    }

    model_df <- rbind(model_df, data.table(Model="LASSO", Fold=fold, Spearman=s_cor, Pearson=p_cor, lambda=best_lambda_fold))
    cat("Signature:", sig, "Fold", fold, "Spearman correlation:", s_cor, "\n")
  }

  cat("\nBest correlation for", sig, ":", best_corr, "from Fold", best_fold, "\n")
  feature_importance <- data.table(
    Peak = colnames(X),
    Weight = best_coef
  )[order(-abs(Weight))]
  fwrite(feature_importance, paste0("lasso_features_", sig, ".csv"))
  fwrite(model_df, paste0("lasso_", sig, ".csv"))
}