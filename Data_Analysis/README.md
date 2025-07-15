# Data Analysis

## Overview

This directory contains the core analytical components of the radiogenomics pipeline, encompassing differential gene expression analysis, pathway enrichment, correlation analyses, survival modeling, and machine learning approaches. These scripts form the analytical backbone for integrating genomic, radiomic, and clinical data to identify clinically relevant radiogenomic associations.

## Purpose

The Data_Analysis directory provides comprehensive analytical tools for:
- **Multi-omics Integration**: Correlating genomic signatures with radiomic features
- **Clinical Relevance Assessment**: Identifying biomarkers associated with survival outcomes
- **Pathway-Level Analysis**: Understanding biological mechanisms underlying radiogenomic associations
- **Predictive Modeling**: Building machine learning models for clinical prediction

## Directory Structure

### 📁 Core Analysis Scripts

#### `DEGanalysis.R`
**Differential Gene Expression Analysis**
- Performs DESeq2-based differential expression analysis between tumor and normal samples
- Processes RNA-seq count data for multiple cancer types
- Generates normalized count matrices and statistical summaries
- Applies variance stabilizing transformation and multiple testing correction
- **Output**: DESeq2 results tables with log2 fold changes, p-values, and adjusted p-values

### 📁 Subdirectories

#### [`Correlations/`](./Correlations/)
**Radiogenomic and Clinical Correlation Analysis**

Contains scripts for comprehensive correlation analyses between different data modalities:

- **`correlative_analysis.R`**: Core radiogenomic correlation analysis
  - Computes Spearman correlations between genomic signatures and radiomic features
  - Performs FDR correction for multiple testing
  - Focuses on original and transformed radiomic features

- **`clinical_correlations.R`**: Clinical outcome correlation analysis
  - Correlates genomic/radiomic features with clinical variables
  - Integrates survival data and clinical metadata
  - Identifies clinically relevant biomarkers

- **`clinical_correlation_filter.R`**: Clinical significance filtering
  - Filters correlation results based on clinical relevance thresholds
  - Applies statistical significance criteria
  - Prioritizes features with clinical associations

- **`genomics_self_correlation.R`**: Genomic feature inter-correlation analysis
  - Analyzes correlations within genomic signature space
  - Identifies co-regulated pathway signatures
  - Assesses pathway redundancy and independence

- **`radiomics_self_correlation.R`**: Radiomic feature inter-correlation analysis
  - Computes correlations within radiomic feature space
  - Identifies redundant imaging features
  - Supports feature selection for downstream analysis

#### [`CoxPH_models/`](./CoxPH_models/)
**Cox Proportional Hazards Survival Analysis**

Implements survival analysis using Cox regression models:

- **`Genomic_cox_model.R`**: Genomic survival modeling
  - Builds Cox proportional hazards models using genomic signatures
  - Assesses pathway associations with overall survival
  - Generates hazard ratios and survival statistics

- **`Radiomics_cox_model.R`**: Radiomic survival modeling
  - Implements Cox regression with radiomic features
  - Evaluates imaging biomarkers for survival prediction
  - Performs model validation and performance assessment

#### [`Enrichment/`](./Enrichment/)
**Gene Set Enrichment Analysis (GSEA)**

Comprehensive pathway enrichment analysis across multiple biological databases:

- **`kegg_enrichment_GMT.R`**: KEGG pathway enrichment
- **`hallmark_enrichment_GMT.R`**: MSigDB Hallmark pathway enrichment
- **`reactome_enrichment_GMT.R`**: Reactome pathway enrichment
- **`biocarta_enrichment_GMT.R`**: BioCarta pathway enrichment

Each script implements GSVA (Gene Set Variation Analysis) to:
- Calculate single-sample pathway enrichment scores
- Transform gene expression data to pathway-level signatures
- Enable pathway-focused downstream analysis

#### [`Machine_Learning_Models/`](./Machine_Learning_Models/)
**Predictive Modeling and Feature Selection**

Advanced machine learning approaches for radiogenomic prediction:

- **`linear_regression_model.R`**: Linear regression baseline models
- **`ridge_regression.R`**: Ridge regression with L2 regularization
- **`LASSO_radiomics.R`**: LASSO regression for feature selection
- **`elastic_net_regression_model.R`**: Elastic net combining Ridge and LASSO

These models provide:
- Feature selection capabilities for high-dimensional data
- Regularization to prevent overfitting
- Cross-validation for robust model evaluation
- Predictive performance assessment

## Analysis Workflow

### 1. **Preprocessing Phase**
- Differential expression analysis identifies cancer-specific gene signatures
- Quality control and normalization of multi-omics data

### 2. **Enrichment Phase**
- Gene expression data transformed to pathway-level signatures
- Multiple pathway databases provide comprehensive biological coverage

### 3. **Correlation Phase**
- Radiogenomic associations identified through correlation analysis
- Clinical relevance assessed through survival and metadata correlations

### 4. **Modeling Phase**
- Cox regression models assess survival associations
- Machine learning models provide predictive capabilities

### 5. **Validation Phase**
- Statistical significance testing with multiple comparison correction
- Clinical validation through outcome associations

## Key Features

### 🔬 **Multi-Modal Integration**
- Seamless integration of genomic, radiomic, and clinical data
- Standardized sample matching and harmonization procedures

### 📊 **Statistical Rigor**
- Multiple testing correction (FDR, Bonferroni)
- Cross-validation for model robustness
- Appropriate statistical tests for each data type

### 🎯 **Clinical Focus**
- Survival analysis integration throughout pipeline
- Clinical metadata incorporation for biomarker validation
- Outcome-driven feature prioritization

### ⚡ **Scalability**
- Designed for multiple cancer types and large datasets
- Efficient computational algorithms for high-dimensional data
- Parallelizable workflows for enhanced performance

## Dependencies

### R Packages
- **Statistical Analysis**: `DESeq2`, `survival`, `survminer`
- **Machine Learning**: `glmnet`, `caret`, `randomForest`
- **Pathway Analysis**: `GSVA`, `GSEABase`
- **Data Manipulation**: `data.table`, `dplyr`, `tidyr`
- **Visualization**: `ggplot2`, `pheatmap`, `corrplot`

### Data Requirements
- **Genomic Data**: RNA-seq count matrices, normalized expression data
- **Radiomic Data**: Quantitative imaging features from medical images
- **Clinical Data**: Survival outcomes, demographic and clinical variables
- **Pathway Data**: GMT files for gene set definitions

## Usage Guidelines

### Input Data Format
- All data files should use consistent sample identifiers
- CSV format preferred for data matrices
- Samples as rows, features as columns convention

### Output Organization
- Results organized by cancer type and analysis method
- Standardized file naming conventions
- Comprehensive logging and statistical summaries

### Quality Control
- Sample size requirements for statistical power
- Missing data handling strategies
- Outlier detection and management protocols

## Integration with Pipeline

The Data_Analysis directory integrates with:
- **[`Preprocessing/`](../Preprocessing/)**: Data harmonization and filtering
- **[`Snakemake/`](../Snakemake/)**: Automated workflow execution
- **[`Visualization_scripts/`](../Visualization_scripts/)**: Results visualization

## Contributing

When adding new analysis scripts:
1. Follow established naming conventions
2. Include comprehensive header documentation
3. Implement appropriate statistical methods
4. Add validation and quality control checks
5. Update this README with new functionality

---

**Note**: This directory represents the analytical core of the radiogenomics pipeline. Each script is designed for reproducibility and scalability across multiple cancer types and datasets.
