# Clinical Associations Pipeline (Snakemake)

## Overview

This Snakemake pipeline automates clinical outcome association analysis for radiogenomics studies. It performs Cox proportional hazards regression to identify radiomic and genomic features (from four pathway collections) that are significantly associated with patient survival outcomes, while controlling for clinical covariates such as age, gender, and tumor stage.

## Pipeline Workflow

The pipeline consists of the following sequential steps:

1. **Clinical Data Extraction** - Extracts and filters clinical data for specific cancer types and treatment modalities
2. **Data Harmonization** - Aligns sample IDs across clinical, radiomic, and genomic datasets
3. **Cox Proportional Hazards Modeling** - Performs survival analysis for:
   - **Genomic Features**: Pathway enrichment scores from 4 collections (Hallmark, KEGG, Reactome, BioCarta)
   - **Radiomic Features**: Quantitative imaging features extracted from medical images
4. **Statistical Analysis** - Calculates hazard ratios, confidence intervals, and FDR-corrected p-values

## File Structure

```
Snakemake/clinical_associations/
├── association_snakefile.snakefile    # Main pipeline definition
├── clinical_config.yaml               # Configuration file with paths and samples
├── clinical_outcome.sh                # SLURM submission script
└── README.md                          # This documentation
```

## Prerequisites

### Software Requirements
- **Snakemake** (workflow management)
- **R** with the following packages:
  - `survival` - Cox proportional hazards modeling
  - `survminer` - Survival analysis visualization
  - `data.table` - Fast data manipulation

### Data Requirements

#### Input Files (specified in config):

1. **Clinical Data**: Combined clinical outcome files
   - Format: CSV with patient demographics, treatment info, and survival data
   - Required columns:
     - `cases.submitter_id` - Patient identifiers
     - `treatments.treatment_type` - Treatment modality
     - `OS_days` - Overall survival time in days
     - `OS_event` - Survival event indicator (0/1)
     - `diagnoses.age_at_diagnosis` - Age at diagnosis
     - `demographic.gender` - Patient gender
     - `diagnoses.ajcc_pathologic_stage` - Tumor stage
   - Path example: `/cluster/projects/bhklab/rawdata/Radiogenomics/Clinical_data/combined_clinical/{cancer}_combined_clinical_data.csv`

2. **Sample ID Lists**: Patient IDs for analysis inclusion
   - Format: One patient ID per line
   - Path example: `/cluster/projects/bhklab/rawdata/Radiogenomics/sample_ids/{cancer}_patientids.txt`

3. **Radiomic Features**: Quantitative imaging features
   - Format: CSV with samples as rows, radiomic features as columns
   - Features filtered by prefix: `original_`, `wavelet-`, `square_`, `squareroot_`, `logarithm_`, `exponential_`, `gradient_`
   - Path example: `/cluster/projects/radiomics/PublicDatasets/procdata/{site}/{cancer}/features/pyradiomics/original_full_features.csv`

4. **Genomic Features**: Pathway enrichment scores (4 collections)
   - **Hallmark**: 50 well-characterized biological processes
   - **KEGG**: Metabolic and signaling pathways
   - **Reactome**: Comprehensive pathway database
   - **BioCarta**: Cancer-focused pathways
   - Format: CSV with samples as rows, pathway scores as columns
   - Path example: `/cluster/projects/bhklab/procdata/Radiogenomics/outputs/correlations/{cancer}/{cancer}_{pathway}_features_filtered.csv`

#### R Scripts (must be available in `/cluster/home/t138199uhn/Scripts/`):
- `clinical_data_extraction.R` - Clinical data filtering and extraction
- `Clinical_filter_uniqueID.R` - Multi-omics data harmonization
- `Genomic_cox_model.R` - Cox regression for genomic features
- `Radiomics_cox_model.R` - Cox regression for radiomic features

## Configuration

### `clinical_config.yaml`

This file contains all pipeline parameters:

```yaml
# Cancer types to analyze
cancer_types:
  - CCRCC    # Clear Cell Renal Cell Carcinoma
  - HNSCC    # Head and Neck Squamous Cell Carcinoma  
  - PDA      # Pancreatic Ductal Adenocarcinoma

# Clinical data file paths
clinical_data:
  CANCER_TYPE: "/path/to/clinical/data.csv"

# Sample ID file paths
sample_ids:
  CANCER_TYPE: "/path/to/sample/ids.txt"

# Radiomic feature file paths
radiomics_file:
  CANCER_TYPE: "/path/to/radiomic/features.csv"

# Genomic feature file paths (4 pathway collections)
genomics_file:
  CANCER_TYPE:
    HALLMARK: "/path/to/hallmark/features.csv"
    KEGG: "/path/to/kegg/features.csv"
    REACTOME: "/path/to/reactome/features.csv"
    BIOCARTA: "/path/to/biocarta/features.csv"

# Output directories
output_dir:
  CANCER_TYPE: "/cluster/projects/bhklab/procdata/Radiogenomics/clinical/{cancer}/"
```

**Supported Cancer Types:**
- **CCRCC** - Clear Cell Renal Cell Carcinoma (CPTAC)
- **HNSCC** - Head and Neck Squamous Cell Carcinoma (CPTAC)
- **PDA** - Pancreatic Ductal Adenocarcinoma (CPTAC)
- **BRCA** - Breast Invasive Carcinoma (TCGA) *[configurable]*
- **GBM** - Glioblastoma Multiforme (TCGA) *[configurable]*
- **KIRC** - Kidney Renal Clear Cell Carcinoma (TCGA) *[configurable]*
- **LGG** - Lower Grade Glioma (TCGA) *[configurable]*

## Usage

### 1. Update Configuration
Edit `clinical_config.yaml` to specify:
- Cancer types to analyze (`cancer_types`)
- Input file paths for your datasets
- Output directory locations

### 2. Submit Pipeline

#### Option A: Direct Snakemake Execution
```bash
cd /path/to/snakemake/clinical_associations
snakemake --snakefile association_snakefile.snakefile \
  --configfile clinical_config.yaml \
  --cores 8
```

#### Option B: SLURM Cluster Submission (Recommended)
```bash
sbatch clinical_outcome.sh
```

### 3. Monitor Progress
- Check SLURM log files in `/cluster/home/t138199uhn/slurm/slurm_logs/clinical_outcomes/`
- Monitor job status: `squeue -u $USER`

## Output Files

All outputs are organized by cancer type in: `/cluster/projects/bhklab/procdata/Radiogenomics/clinical/{cancer}/`

### Data Harmonization Outputs

#### Harmonized Datasets:
- **Clinical Data**: `{cancer}_harmonized_clinical.csv`
  - Filtered clinical data with standardized sample IDs
  - Includes survival endpoints and clinical covariates

- **Radiomic Data**: `{cancer}_harmonized_radiomics.csv`
  - Radiomic features aligned to harmonized sample IDs
  - Features filtered by allowed prefixes

- **Genomic Data** (4 files per cancer):
  - `{cancer}_harmonized_hallmark_genomics.csv`
  - `{cancer}_harmonized_kegg_genomics.csv`
  - `{cancer}_harmonized_reactome_genomics.csv`
  - `{cancer}_harmonized_biocarta_genomics.csv`

### Cox Regression Results

#### Genomic Cox Models (8 files per cancer):
- **Binary Models**: `{cancer}_{PATHWAY}_cox_results_binary.csv`
- **Continuous Models**: `{cancer}_{PATHWAY}_cox_results_continuous.csv`

Where `{PATHWAY}` includes: HALLMARK, KEGG, REACTOME, BIOCARTA

#### Radiomic Cox Models (2 files per cancer):
- **Binary Model**: `{cancer}_radiomics_cox_results_binary.csv`
- **Continuous Model**: `{cancer}_radiomics_cox_results_continuous.csv`

### Results File Structure

Each Cox regression results file contains:

| Column | Description |
|--------|-------------|
| `Feature/Signature` | Feature or pathway name |
| `HR` | Hazard Ratio |
| `CI_lower` | Lower bound of 95% confidence interval |
| `CI_upper` | Upper bound of 95% confidence interval |
| `p_value` | Statistical significance |
| `FDR_p_value` | FDR-corrected p-value for multiple testing |
| `C_index` | Concordance index (model performance) |

## Resource Requirements

### Computational Resources per Rule:

| Rule | Memory | Runtime | Purpose |
|------|--------|---------|---------|
| clinical_data_extraction | 8GB | 1 hour | Clinical data filtering |
| clinical_filter_unique_id | 4GB | 30 min | Data harmonization |
| genomic_cox_hallmark | 16GB | 14 hours | Hallmark Cox models |
| genomic_cox_kegg | 16GB | 14 hours | KEGG Cox models |
| genomic_cox_reactome | 16GB | 14 hours | Reactome Cox models |
| genomic_cox_biocarta | 16GB | 14 hours | BioCarta Cox models |
| radiomics_cox_model | 16GB | 14.5 hours | Radiomic Cox models |

**Note**: Cox modeling steps require substantial memory due to the large number of features and iterative model fitting.

## Methodology

### Clinical Data Processing
1. **Extraction**: Filters clinical data for specific treatment types:
   - Chemotherapy
   - Radiation Therapy, NOS
2. **Harmonization**: Aligns sample IDs across all data modalities
3. **Covariate Preparation**: Standardizes clinical variables for modeling

### Cox Proportional Hazards Modeling

#### Feature Processing:
- **Binary Analysis**: Features dichotomized using median cutoff
- **Continuous Analysis**: Features used as continuous variables
- **Covariate Adjustment**: All models control for:
  - Age at diagnosis (converted from days to years)
  - Gender
  - AJCC pathologic stage

#### Statistical Analysis:
- **Univariate Models**: Each feature modeled independently
- **Hazard Ratios**: Calculated with 95% confidence intervals
- **Multiple Testing Correction**: FDR (Benjamini-Hochberg) applied
- **Model Performance**: Concordance index (C-index) calculated

### Quality Control
- Sample ID consistency across all data types
- Missing value handling and validation
- Treatment type filtering for analysis relevance
- Minimum sample size requirements

## Troubleshooting

### Common Issues:

1. **Missing Required Columns**
   ```
   Error: Clinical data is missing required columns
   ```
   - Verify clinical data contains all required survival and covariate columns
   - Check column naming conventions match expected format

2. **Sample ID Mismatches**
   ```
   Error: No matching sample IDs found
   ```
   - Ensure consistent sample ID formats across all input files
   - Check for leading/trailing whitespace in sample ID files

3. **Insufficient Memory**
   ```
   Error: Memory allocation failed
   ```
   - Increase memory allocation in Snakefile resources
   - Consider running fewer cancer types in parallel

4. **Missing Input Files**
   ```
   Error: File not found
   ```
   - Verify all file paths in `clinical_config.yaml`
   - Check file permissions and accessibility on cluster

### Log Files:
- **Master job**: `clinical_pipeline_{jobid}.out/err`
- **Individual rules**: `{rule}.{cancer}.{jobid}.out/err`

### Debugging Tips:
1. **Check data quality**: Verify survival times and event indicators
2. **Validate harmonization**: Ensure sample counts match across data types
3. **Review convergence**: Check for Cox model convergence warnings
4. **Examine feature distributions**: Verify feature variance and outliers

## Integration with Radiogenomics Pipeline

The Cox regression results serve as input for:
- **Feature Prioritization**: Identifying survival-associated features for downstream analysis
- **Biomarker Discovery**: Finding radiomic/genomic signatures predictive of outcomes
- **Model Selection**: Choosing features for integrative machine learning models
- **Clinical Translation**: Developing clinically relevant prognostic signatures

### Downstream Analysis Connections:
- **Correlation Analysis**: Using survival-associated features for radiogenomic correlations
- **Machine Learning**: Incorporating significant features into predictive models
- **Visualization**: Creating survival curves and forest plots for significant associations
- **Meta-Analysis**: Combining results across multiple cancer types

## Clinical Interpretation

### Hazard Ratio Interpretation:
- **HR > 1**: Feature associated with increased risk (poor prognosis)
- **HR < 1**: Feature associated with decreased risk (good prognosis)
- **HR = 1**: No association with survival

### Significance Thresholds:
- **p-value < 0.05**: Nominally significant
- **FDR p-value < 0.05**: Significant after multiple testing correction
- **FDR p-value < 0.10**: Suggestive association

### Model Performance:
- **C-index > 0.7**: Good discriminative ability
- **C-index 0.6-0.7**: Moderate discriminative ability  
- **C-index < 0.6**: Poor discriminative ability

## Authors
- Jackie Chen
- Julia Nguyen

## Version
- Pipeline Version: 2024.1
- Last Updated: July 2025

## References

