# Ridge regression with cross-validation and hyperparameter tuning (R version)
# Matches the structure and input/output of your LASSO and linear regression R scripts

library(data.table)
library(caret)
library(glmnet)

set.seed(42)

# Load data
X <- fread("X.csv")
gene_matrix <- fread("gene_signatures.csv") # samples x all gene signatures

# Specify the gene signature columns you want to use (replace with your actual column names)
selected_signatures <- c("Signature1", "Signature2", "Signature3")
Y <- gene_matrix[, ..selected_signatures]

# Ridge parameter grid (alpha=0 for Ridge, lambda is the penalty)
lambda_grid <- c(0.1, 1, 10, 100) # match your Python grid

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

    # Fit Ridge with cross-validation for lambda
    cvfit <- cv.glmnet(X_train, y_train, alpha=0, lambda=lambda_grid, nfolds=5)
    best_lambda_fold <- cvfit$lambda.min
    ridge_model <- glmnet(X_train, y_train, alpha=0, lambda=best_lambda_fold)

    y_pred <- predict(ridge_model, newx=X_test, s=best_lambda_fold)
    s_cor <- cor(y_pred, y_test, method = "spearman")
    p_cor <- cor(y_pred, y_test, method = "pearson")

    if (s_cor > best_corr) {
      best_corr <- s_cor
      best_fold <- fold
      best_coef <- as.vector(coef(ridge_model))[-1] # drop intercept
      best_lambda <- best_lambda_fold
    }

    model_df <- rbind(model_df, data.table(Model="Ridge", Fold=fold, Spearman=s_cor, Pearson=p_cor, lambda=best_lambda_fold))
    cat("Signature:", sig, "Fold", fold, "Spearman correlation:", s_cor, "\n")
  }

  cat("\nBest correlation for", sig, ":", best_corr, "from Fold", best_fold, "\n")
  feature_importance <- data.table(
    Peak = colnames(X),
    Weight = best_coef
  )[order(-abs(Weight))]
  fwrite(feature_importance, paste0("ridge_features_", sig, ".csv"))
  fwrite(model_df, paste0("ridge_", sig, ".csv"))
}
