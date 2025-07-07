library(survival)
library(survminer)

# ---- COMMAND LINE ARGUMENTS ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript Radiomics_cox_model.R <radiomics_features_file> <clinical_data_file> <output_file> <treatment_type>")
}

radiomics_features_file <- args[1]
clinical_data_file <- args[2]
output_file <- args[3]
treatment_type <- args[4]

# ---- LOAD DATA ----
radiomics_features <- read.csv(radiomics_features_file, row.names = 1)
clinical_data <- read.csv(clinical_data_file)

# ---- FILTER RADIOMICS FEATURES BASED ON PREFIXES ----
# Define allowed prefixes
allowed_prefixes <- c("patient_ID", "original-", "wavelet-", "square_", "squareroot_", "logarithm_", "exponential_", "gradient_")

# Filter columns based on prefixes
filtered_features <- radiomics_features[, grepl(paste0("^(", paste(allowed_prefixes, collapse = "|"), ")"), colnames(radiomics_features))]
cat(sprintf("Filtered radiomics features based on allowed prefixes. Remaining columns: %d\n", ncol(filtered_features)))

# ---- FILTER CLINICAL DATA BASED ON TREATMENT TYPE ----
if (!"treatments.treatment_type" %in% colnames(clinical_data)) {
  stop("Clinical data is missing the 'treatments.treatment_type' column.")
}
filtered_clinical_data <- clinical_data[clinical_data$treatments.treatment_type == treatment_type, ]
cat(sprintf("Filtered clinical data for treatment type '%s'. Remaining samples: %d\n", treatment_type, nrow(filtered_clinical_data)))

# ---- CONVERT AGE FROM DAYS TO YEARS ----
if ("diagnoses.age_at_diagnosis" %in% colnames(filtered_clinical_data)) {
  filtered_clinical_data$diagnoses.age_at_diagnosis <- filtered_clinical_data$diagnoses.age_at_diagnosis / 365.25
  cat("Converted age from days to years.\n")
} else {
  stop("Clinical data is missing the 'diagnoses.age_at_diagnosis' column.")
}

# ---- FILTER RADIOMICS FEATURES BASED ON SAMPLES ----
filtered_radiomics_features <- filtered_features[rownames(filtered_features) %in% filtered_clinical_data$cases.submitter_id, ]
cat(sprintf("Filtered radiomics features based on treatment type '%s'. Remaining samples: %d\n", treatment_type, nrow(filtered_radiomics_features)))

# ---- BINARIZE RADIOMICS FEATURES ----
binarized_features <- filtered_radiomics_features
for (feature in colnames(filtered_radiomics_features)) {
  median_value <- median(filtered_radiomics_features[[feature]], na.rm = TRUE)
  binarized_features[[feature]] <- ifelse(filtered_radiomics_features[[feature]] > median_value, 1, 0)
}
cat("Radiomics features have been binarized based on median scores.\n")

# ---- ENSURE CLINICAL DATA CONTAINS REQUIRED COLUMNS ----
required_columns <- c("OS_days", "diagnoses.age_at_diagnosis", "demographic.gender", "diagnoses.ajcc_pathologic_stage")
if (!all(required_columns %in% colnames(filtered_clinical_data))) {
  stop("Clinical data is missing required columns: OS_days, diagnoses.age_at_diagnosis, demographic.gender, diagnoses.ajcc_pathologic_stage")
}

# ---- MERGE DATA ----
merged_data <- merge(filtered_clinical_data, binarized_features, by.x = "cases.submitter_id", by.y = "row.names", all.x = TRUE)
rownames(merged_data) <- merged_data$cases.submitter_id

# Ensure no missing values in OS_days
if (any(is.na(merged_data$OS_days))) {
  stop("Missing values detected in OS_days column after merging.")
}

# ---- PREPARE SURVIVAL OBJECT ----
surv_object <- Surv(time = merged_data$OS_days)  # Continuous survival time

# ---- COX PROPORTIONAL HAZARDS MODEL ----
# Initialize results dataframe
results <- data.frame(
  Feature = character(),
  HR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through binarized radiomics features
for (feature in colnames(binarized_features)) {
  # Ensure no missing values in the feature column
  if (any(is.na(merged_data[[feature]]))) {
    cat("Skipping feature due to missing values:", feature, "\n")
    next
  }

  # Fit Cox model controlling for age, sex, and tumour stage
  cox_model <- coxph(surv_object ~ merged_data[[feature]] + diagnoses.age_at_diagnosis + demographic.gender + diagnoses.ajcc_pathologic_stage, data = merged_data)
  
  # Extract results
  summary_cox <- summary(cox_model)
  hr <- summary_cox$coefficients[1, "exp(coef)"]
  ci <- summary_cox$conf.int[1, c("lower .95", "upper .95")]
  p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
  
  # Append to results dataframe
  results <- rbind(results, data.frame(
    Feature = feature,
    HR = hr,
    CI_lower = ci[1],
    CI_upper = ci[2],
    p_value = p_value
  ))
}

# ---- SAVE RESULTS ----
write.csv(results, output_file, row.names = FALSE)
cat("Cox proportional hazards model results saved to:", output_file, "\n")