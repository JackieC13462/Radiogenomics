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


---

# Pan-Cancer Survival Analysis Pipeline (Snakemake)

## Overview

This Snakemake pipeline performs pan-cancer survival analysis by combining genomic pathway enrichment data and clinical outcomes across multiple cancer types. The pipeline separates compiled genomic enrichment scores by pathway groups, harmonizes sample IDs across datasets, combines data across cancer types, and performs Cox proportional hazards regression to identify pathway signatures significantly associated with patient survival outcomes.

## Pipeline Workflow

The pipeline consists of the following sequential steps:

1. **Clinical Data Extraction** - Extracts and filters clinical data for specific cancer types and treatment modalities
2. **Genomic Data Separation & Harmonization** - Separates compiled enrichment scores by pathway groups (HALLMARK, KEGG, REACTOME, BIOCARTA) and harmonizes sample IDs
3. **Pan-Cancer Data Combination** - Combines harmonized datasets across all cancer types for pan-cancer analysis
4. **Pan-Cancer Cox Regression** - Performs survival analysis for each pathway group across all cancer types with covariate adjustment

## File Structure

```
Snakemake/clinical_associations/
├── pancancer_snakefile.snakefile     # Main pan-cancer pipeline definition
├── pancancer_config.yaml             # Configuration file with paths and samples
├── pancancer_outcome.sh              # SLURM submission script for pan-cancer analysis
└── README.md                         # This documentation
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

1. **Clinical Data**: Combined clinical outcome files per cancer type
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

2. **Sample ID Lists**: Patient IDs for analysis inclusion per cancer type
   - Format: One patient ID per line
   - Path example: `/cluster/projects/bhklab/rawdata/Radiogenomics/sample_ids/{dataset}_{cancer}_patientids.txt`

3. **Compiled Genomic Features**: Pathway enrichment scores (4 collections combined)
   - **Input Format**: CSV with pathways as rows, samples as columns
   - **Pathway Groups**: Identified by prefixes (HALLMARK_, KEGG_, REACTOME_, BIOCARTA_)
   - **Processing**: Script automatically separates pathway groups and transposes to samples × pathways format
   - **Collections**:
     - **Hallmark**: 50 well-characterized biological processes
     - **KEGG**: Metabolic and signaling pathways
     - **Reactome**: Comprehensive pathway database
     - **BioCarta**: Cancer-focused pathways
   - Path example: `/cluster/projects/bhklab/procdata/Radiogenomics/gene_signatures/{cancer}_compiled_enrichment.csv`

#### R Scripts (must be available in `/cluster/home/t138199uhn/Scripts/`):
- `clinical_data_extraction.R` - Clinical data filtering and extraction
- `gen_clin_newID.R` - Genomic data separation, transposition, and harmonization
- `combine_datasets.R` - Pan-cancer dataset combination
- `pancancer_gen_cox.R` - Pan-cancer Cox regression analysis

## Configuration

### `pancancer_config.yaml`

This file contains all pipeline parameters:

```yaml
cancer_types:
  - BRCA
  - LGG
  - KIRC
  - CCRCC
  - HNSCC
  - PDA

clinical_data:
  {cancer}: "/path/to/{cancer}_combined_clinical_data.csv"

sample_ids:
  {cancer}: "/path/to/{dataset}_{cancer}_patientids.txt"

genomics_file:
  {cancer}: "/path/to/{cancer}_compiled_enrichment.csv"

output_dir:
  {cancer}: "/cluster/projects/bhklab/procdata/Radiogenomics/clinical/PANCAN/{cancer}/"
```

#### Key Parameters:
- **cancer_types**: List of cancer types to include in pan-cancer analysis
- **clinical_data**: Paths to clinical data files per cancer type  
- **sample_ids**: Paths to sample ID lists per cancer type
- **genomics_file**: Paths to compiled genomic enrichment files per cancer type
- **output_dir**: Output directories for harmonized files per cancer type

## Usage

### Command Line Execution

#### Quick Start:
```bash
# Navigate to pipeline directory
cd /cluster/home/t138199uhn/Scripts

# Submit to SLURM cluster
sbatch pancancer_outcome.sh
```

#### Manual Execution:
```bash
# Load required modules
module load snakemake R

# Run pipeline
snakemake --snakefile pancancer_snakefile.snakefile \
          --configfile pancancer_config.yaml \
          --jobs 50 \
          --cluster "sbatch --mem={resources.mem_mb} --time={resources.time}"
```

### Input Data Preparation

1. **Prepare Clinical Data**: Ensure clinical files contain required survival and covariate columns
2. **Prepare Compiled Genomics**: Create single file per cancer type with all pathway enrichment scores
3. **Update Configuration**: Edit `pancancer_config.yaml` with correct file paths
4. **Verify Sample IDs**: Ensure sample ID lists match those in clinical and genomic data

## Output Files

### Individual Cancer Type Harmonized Datasets:
Located in cancer-specific subdirectories (`/PANCAN/{cancer}/`):

- **Clinical Data**: `{cancer}_harmonized_clinical.csv`
  - Filtered clinical data with standardized sample IDs
  - Includes survival endpoints and clinical covariates

- **Genomic Data** (4 files per cancer):
  - `{cancer}_harmonized_hallmark_genomics.csv`
  - `{cancer}_harmonized_kegg_genomics.csv`
  - `{cancer}_harmonized_reactome_genomics.csv`
  - `{cancer}_harmonized_biocarta_genomics.csv`
  - Format: Samples × pathways with harmonized sample IDs

### Combined Pan-Cancer Datasets:
Located in main PANCAN directory:

- **Combined Clinical**: `COMBINED_harmonized_clinical.csv`
  - Clinical data from all cancer types with cancer_type column added
  - Standardized column names (suffixes removed)

- **Combined Genomics** (4 files):
  - `COMBINED_harmonized_kegg_genomics.csv`
  - `COMBINED_harmonized_hallmark_genomics.csv`
  - `COMBINED_harmonized_reactome_genomics.csv`
  - `COMBINED_harmonized_biocarta_genomics.csv`
  - Format: All samples × pathways across cancer types

### Pan-Cancer Cox Regression Results:
Located in main PANCAN directory:

- **Survival Analysis Results** (4 files):
  - `COMBINED_KEGG_pancancer_cox_results_binary.csv`
  - `COMBINED_HALLMARK_pancancer_cox_results_binary.csv`
  - `COMBINED_REACTOME_pancancer_cox_results_binary.csv`
  - `COMBINED_BIOCARTA_pancancer_cox_results_binary.csv`

#### Cox Results Format:
- **Signature**: Pathway name
- **HR**: Hazard ratio
- **CI_lower, CI_upper**: 95% confidence interval bounds
- **p_value**: Statistical significance
- **C_index**: Concordance index (model performance)
- **FDR**: False discovery rate corrected p-value

## Methodology

### Genomic Data Processing
1. **Pathway Separation**: Compiled enrichment files separated by pathway group prefixes
2. **Data Transposition**: Converts from pathways × samples to samples × pathways format
3. **Harmonization**: Aligns sample IDs between clinical and genomic datasets
4. **Quality Control**: Removes samples missing in either clinical or genomic data

### Pan-Cancer Data Combination
1. **Clinical Integration**: Combines clinical data across cancer types with standardized column names
2. **Genomic Integration**: Merges genomic data by pathway group across cancer types
3. **Cancer Type Stratification**: Adds cancer_type column for downstream analysis

### Cox Proportional Hazards Modeling

#### Feature Processing:
- **Binary Analysis**: Pathways dichotomized using median cutoff across all samples
- **Survival Endpoints**: Uses both OS_days (time) and OS_event (event status)
- **Covariate Adjustment**: All models control for:
  - Age at diagnosis
  - Sex  
  - Tumour pathologic stage
  - Cancer type 

#### Statistical Analysis:
- **Univariate Models**: Each pathway modeled independently with covariate adjustment
- **Hazard Ratios**: Calculated with 95% confidence intervals
- **Multiple Testing Correction**: FDR (Benjamini-Hochberg) applied within each pathway group
- **Model Performance**: Concordance index (C-index) calculated for each model

## Troubleshooting

### Common Issues:

1. **Missing Sample IDs**: Ensure sample IDs are consistent across clinical and genomic files
2. **Pathway Separation Failures**: Verify pathway names contain correct prefixes (HALLMARK_, KEGG_, etc.)
3. **Memory Issues**: Increase memory allocation for large datasets in SLURM script
4. **Column Name Mismatches**: Check that clinical data contains all required columns

### Debug Mode:
```bash
# Run with verbose output
snakemake --snakefile pancancer_snakefile.snakefile \
          --configfile pancancer_config.yaml \
          --verbose \
          --dry-run
```

## Version
- Pan-Cancer Pipeline Version: 2025.1
- Last Updated: July 2025
- Compatible with: Snakemake ≥7.0, R ≥4.0


