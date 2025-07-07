# Elastic Net regression with cross-validation and hyperparameter tuning (R version)
# Matches the structure and input/output of your other ML R scripts

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

# Elastic Net parameter grid
lambda_grid <- c(0.1, 1, 10, 100)
l1_ratio_grid <- c(0.2, 0.5, 0.8)

for (sig in selected_signatures) {
  y <- Y[[sig]]
  model_df <- data.table(Model=character(), Fold=integer(), Spearman=numeric(), Pearson=numeric(), lambda=numeric(), l1_ratio=numeric())
  folds <- createFolds(y, k=10, list=TRUE, returnTrain=FALSE)
  best_corr <- -Inf
  best_fold <- 0
  best_coef <- NULL
  best_lambda <- NA
  best_l1_ratio <- NA

  for (fold in seq_along(folds)) {
    test_idx <- folds[[fold]]
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)
    X_train <- as.matrix(X[train_idx, , drop=FALSE])
    y_train <- y[train_idx]
    X_test  <- as.matrix(X[test_idx, , drop=FALSE])
    y_test  <- y[test_idx]

    # Grid search for best alpha (l1_ratio) and lambda
    best_fold_corr <- -Inf
    best_fold_coef <- NULL
    best_fold_lambda <- NA
    best_fold_l1_ratio <- NA
    for (l1_ratio in l1_ratio_grid) {
      cvfit <- cv.glmnet(X_train, y_train, alpha=l1_ratio, lambda=lambda_grid, nfolds=5)
      lambda_min <- cvfit$lambda.min
      en_model <- glmnet(X_train, y_train, alpha=l1_ratio, lambda=lambda_min)
      y_pred <- predict(en_model, newx=X_test, s=lambda_min)
      s_cor <- cor(y_pred, y_test, method = "spearman")
      p_cor <- cor(y_pred, y_test, method = "pearson")
      if (s_cor > best_fold_corr) {
        best_fold_corr <- s_cor
        best_fold_coef <- as.vector(coef(en_model))[-1]
        best_fold_lambda <- lambda_min
        best_fold_l1_ratio <- l1_ratio
        best_fold_pcor <- p_cor
      }
    }
    if (best_fold_corr > best_corr) {
      best_corr <- best_fold_corr
      best_fold <- fold
      best_coef <- best_fold_coef
      best_lambda <- best_fold_lambda
      best_l1_ratio <- best_fold_l1_ratio
    }
    model_df <- rbind(model_df, data.table(Model="ElasticNet", Fold=fold, Spearman=best_fold_corr, Pearson=best_fold_pcor, lambda=best_fold_lambda, l1_ratio=best_fold_l1_ratio))
    cat("Signature:", sig, "Fold", fold, "Spearman correlation:", best_fold_corr, "\n")
  }

  cat("\nBest correlation for", sig, ":", best_corr, "from Fold", best_fold, "\n")
  feature_importance <- data.table(
    Peak = colnames(X),
    Weight = best_coef
  )[order(-abs(Weight))]
  fwrite(feature_importance, paste0("en_features_", sig, ".csv"))
  fwrite(model_df, paste0("en_", sig, ".csv"))
}
