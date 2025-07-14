# Radiogenomic Correlations Pipeline (Snakemake)

## Overview

This Snakemake pipeline automates comprehensive radiogenomic correlation analysis to identify associations between quantitative imaging features (radiomics) and genomic pathway signatures across multiple cancer datasets. The pipeline performs systematic correlation analysis between radiomic features and pathway enrichment scores from four major gene set collections, while also analyzing clinical outcome associations.

## Pipeline Workflow

The pipeline consists of the following sequential steps:

1. **Sample ID Harmonization** - Aligns sample identifiers between genomic and radiomic datasets
2. **Pathway Separation** - Separates compiled genomic enrichment scores by pathway source (KEGG, HALLMARK, REACTOME, BIOCARTA)
3. **Clinical Data Integration** - Harmonizes clinical outcome data with multi-omics datasets
4. **Self-Correlation Analysis** - Computes correlation matrices within radiomic and genomic feature sets
5. **Feature Filtering** - Removes highly correlated features (|r| > 0.9) to reduce redundancy
6. **Radiogenomic Correlation** - Performs Spearman correlation analysis between filtered radiomic and genomic features
7. **Clinical Correlation Analysis** - Analyzes correlations between genomic/radiomic signatures and survival outcomes
8. **Statistical Filtering** - Applies significance and correlation strength thresholds to identify meaningful associations

## File Structure

```
Snakemake/Correlations/
├── correlative_snakefile.snakefile    # Main pipeline definition
├── correlative_config.yaml            # Configuration file with dataset paths
├── correlative_snakefile.sh           # SLURM submission script
└── README.md                          # This documentation
```

## Prerequisites

### Software Requirements
- **Snakemake** (workflow management)
- **R** with the following packages:
  - `data.table` - Fast data manipulation
  - `corrplot` - Correlation visualization
  - `ggplot2` - Data visualization
  - `dplyr` - Data manipulation

### Data Requirements

#### Input Files (specified in config):

1. **Compiled Genomic Enrichment Scores**: Combined pathway scores from enrichment pipeline
   - Format: CSV with samples as rows, pathway signatures as columns
   - Content: GSVA enrichment scores from all four pathway collections (Hallmark, KEGG, Reactome, BioCarta)
   - Path example: `/cluster/projects/bhklab/procdata/Radiogenomics/gene_signatures/{cancer}_compiled_enrichment.csv`

2. **Radiomic Features**: Quantitative imaging features extracted using PyRadiomics
   - Format: CSV with samples as rows, radiomic features as columns
   - Feature categories: Original, wavelet-transformed, and mathematical transformations
   - Allowed prefixes: `original_`, `wavelet-`, `square_`, `squareroot_`, `logarithm_`, `exponential_`, `gradient_`
   - Path example: `/cluster/projects/radiomics/PublicDatasets/procdata/{site}/{cancer}/features/pyradiomics/original_full_features.csv`

3. **Clinical Data**: Patient clinical outcomes and demographics
   - Format: CSV with patient demographics and survival data
   - Required column: `OS_days` (overall survival in days)
   - Additional: Demographics, treatment information, staging
   - Path example: `/cluster/projects/bhklab/rawdata/Radiogenomics/Clinical_data/combined_clinical/{cancer}_combined_clinical_data.csv`

#### R Scripts (must be available in `/cluster/home/t138199uhn/Scripts/`):
- `unique_ID_generator.R` - Sample ID harmonization between datasets
- `enrichment_separator.R` - Pathway collection separation
- `Clinical_filter_uniqueID.R` - Multi-omics data harmonization
- `radiomics_self_correlation.R` - Radiomic feature correlation analysis
- `genomics_self_correlation.R` - Genomic signature correlation analysis
- `feature_filtering.R` - Redundant feature removal
- `correlative_analysis.R` - Radiogenomic correlation computation
- `clinical_correlations.R` - Clinical outcome correlation analysis
- `clinical_correlation_filter.R` - Statistical significance filtering

## Configuration

### `correlative_config.yaml`

This file contains dataset specifications:

```yaml
datasets:
  - name: CANCER_TYPE
    genomics_file: "/path/to/compiled/enrichment.csv"
    radiomics_file: "/path/to/radiomic/features.csv"
    clinical_file: "/path/to/clinical/data.csv"
```

**Supported Cancer Types:**
- **HNSCC** - Head and Neck Squamous Cell Carcinoma (CPTAC)
- **CCRCC** - Clear Cell Renal Cell Carcinoma (CPTAC)
- **PDA** - Pancreatic Ductal Adenocarcinoma (CPTAC)
- **BRCA** - Breast Invasive Carcinoma (TCGA) *[configurable]*
- **GBM** - Glioblastoma Multiforme (TCGA) *[configurable]*
- **KIRC** - Kidney Renal Clear Cell Carcinoma (TCGA) *[configurable]*
- **LGG** - Lower Grade Glioma (TCGA) *[configurable]*

## Usage

### 1. Update Configuration
Edit `correlative_config.yaml` to specify:
- Cancer datasets to analyze
- Input file paths for genomic, radiomic, and clinical data
- Ensure all file paths are accessible

### 2. Submit Pipeline

#### Option A: Direct Snakemake Execution
```bash
cd /path/to/snakemake/correlations
snakemake --snakefile correlative_snakefile.snakefile \
  --configfile correlative_config.yaml \
  --cores 8
```

#### Option B: SLURM Cluster Submission (Recommended)
```bash
sbatch correlative_snakefile.sh
```

### 3. Monitor Progress
- Check SLURM log files in `/cluster/home/t138199uhn/slurm/slurm_logs/`
- Monitor job status: `squeue -u $USER`

## Output Files

All outputs are organized by cancer type in: `/cluster/projects/bhklab/procdata/Radiogenomics/outputs/correlations/{cancer}/`

### Data Harmonization Outputs

#### Sample ID Alignment:
- **Genomic Data**: `{cancer}_radiogenomic_RNAseq.csv`
  - Genomic signatures with harmonized sample IDs
- **Radiomic Data**: `{cancer}_radiogenomic_features.csv`
  - Radiomic features with harmonized sample IDs

#### Pathway Separation:
- **Hallmark**: `{cancer}_HALLMARK_enrichment.csv`
- **KEGG**: `{cancer}_KEGG_enrichment.csv`
- **Reactome**: `{cancer}_REACTOME_enrichment.csv`
- **BioCarta**: `{cancer}_BIOCARTA_enrichment.csv`

#### Clinical Integration:
- **Harmonized Clinical**: `{cancer}_harmonized_clinical.csv`
- **Harmonized Genomics** (4 files): `{cancer}_harmonized_{pathway}_genomics.csv`
- **Harmonized Radiomics**: `{cancer}_harmonized_radiomics.csv`

### Self-Correlation Analysis

#### Correlation Matrices:
- **Radiomic Self-Correlation**: `{cancer}_radiomics_self_correlation.csv`
- **Genomic Self-Correlations** (4 files): `{cancer}_{PATHWAY}_self_correlation.csv`

### Feature Filtering Outputs

#### Filtered Feature Sets:
- **Filtered Genomics** (4 files): `{cancer}_{PATHWAY}_features_filtered.csv`
- **Filtered Radiomics**: `{cancer}_radiomics_features_filtered.csv`
- **Correlation Reports**: `{cancer}_{PATHWAY}_0.9correlation.csv` (highly correlated pairs)

### Radiogenomic Correlation Results

#### Primary Correlation Analysis (4 files per cancer):
- **File Format**: `{cancer}_{PATHWAY}_correlative_analysis.csv`
- **Content**: Spearman correlations between radiomic features and pathway signatures

#### Results Structure:
| Column | Description |
|--------|-------------|
| `Genomic_Signature` | Pathway or gene signature name |
| `Radiomic_Feature` | Radiomic feature name |
| `Correlation` | Spearman correlation coefficient |
| `P_value` | Statistical significance |
| `P_value_adjusted` | FDR-corrected p-value |
| `Sample_Size` | Number of samples used |

### Clinical Correlation Results

#### Full Clinical Correlations (5 files per cancer):
- **Genomic Pathways** (4 files): `{cancer}_{PATHWAY}_clinical_correlations.csv`
- **Radiomic Features**: `{cancer}_radiomics_clinical_correlations.csv`

#### Filtered Clinical Correlations (5 files per cancer):
- **Genomic Pathways** (4 files): `{cancer}_{PATHWAY}_clinical_correlations_filtered.csv`
- **Radiomic Features**: `{cancer}_radiomics_clinical_correlations_filtered.csv`

#### Clinical Results Structure:
| Column | Description |
|--------|-------------|
| `Signature` | Feature or pathway name |
| `Correlation` | Spearman correlation with OS_days |
| `P_value` | Statistical significance |
| `P_value_adjusted` | FDR-corrected p-value |

## Resource Requirements

### Computational Resources per Rule:

| Rule | Memory | Runtime | Purpose |
|------|--------|---------|---------|
| unique_ID_generator | 8GB | 40 min | Sample harmonization |
| enrichment_separation | 8GB | 30 min | Pathway separation |
| clinical_ID_filtering | 8GB | 30 min | Clinical integration |
| radiomics_self_correlation | 8GB | 2 hours | Radiomic correlation matrix |
| genomics_self_correlation | 8GB | 2 hours | Genomic correlation matrix |
| feature_filtering | 4GB | 50 min | Redundancy removal |
| correlative_analysis | 16GB | 2 days | Radiogenomic correlations |
| clinical_correlation_analysis | 8GB | 3 days | Clinical correlations |
| clinical_correlation_filtering | 4GB | 30 min | Statistical filtering |

**Note**: Correlative analysis steps are computationally intensive due to pairwise correlation calculations across thousands of features.

## Methodology

### Sample Harmonization
1. **ID Alignment**: Matches sample identifiers across genomic, radiomic, and clinical datasets
2. **Data Integration**: Ensures consistent sample ordering and removes samples missing from any modality
3. **Quality Control**: Validates data completeness and consistency

### Feature Processing
1. **Pathway Separation**: Divides compiled enrichment scores by gene set collection
2. **Radiomic Filtering**: Selects clinically relevant feature categories by prefix
3. **Self-Correlation**: Identifies redundant features within each data type
4. **Redundancy Removal**: Filters features with |correlation| > 0.9 to reduce multicollinearity

### Correlation Analysis
1. **Method**: Spearman rank correlation (robust to outliers and non-linear relationships)
2. **Scope**: Pairwise correlations between all filtered radiomic and genomic features
3. **Clinical Analysis**: Correlations with overall survival time (OS_days)
4. **Multiple Testing**: FDR (Benjamini-Hochberg) correction applied

### Statistical Filtering
- **Default Criteria**: |correlation| > 0.35 AND p-value < 1.0
- **Customizable**: Thresholds adjustable in pipeline configuration
- **Output**: Separate files for all correlations and filtered significant associations

## Troubleshooting

### Common Issues:

1. **Sample ID Mismatches**
   ```
   Error: No common samples found between datasets
   ```
   - Verify sample ID format consistency across input files
   - Check for extra characters, spaces, or case differences

2. **Missing Required Columns**
   ```
   Error: 'OS_days' column not found in clinical data
   ```
   - Ensure clinical data contains required survival endpoint
   - Verify column naming matches expected format

3. **Insufficient Memory for Correlation Analysis**
   ```
   Error: Memory allocation failed
   ```
   - Increase memory allocation for correlative_analysis rules
   - Consider running fewer datasets in parallel

4. **Feature Filtering Issues**
   ```
   Warning: No features remain after filtering
   ```
   - Check radiomic feature prefix patterns
   - Verify correlation threshold settings

### Log Files:
- **Master job**: `radiogenomics_pipeline_{jobid}.out/err`
- **Individual rules**: `{rule}.{cancer}.{jobid}.out/err`

### Debugging Tips:
1. **Check data dimensions**: Verify sample and feature counts after each step
2. **Validate harmonization**: Ensure consistent samples across all data types
3. **Review correlation distributions**: Check for reasonable correlation ranges
4. **Monitor convergence**: Watch for memory and time limit issues

## Integration with Radiogenomics Pipeline

The correlation results serve as input for:
- **Biomarker Discovery**: Identifying radiomic-genomic feature pairs for further validation
- **Pathway Analysis**: Understanding biological processes underlying imaging phenotypes
- **Machine Learning**: Selecting correlated features for integrative predictive models
- **Clinical Translation**: Prioritizing associations with survival outcomes

### Downstream Analysis Connections:
- **Survival Analysis**: Using correlated features in Cox proportional hazards models
- **Visualization**: Creating correlation heatmaps and network diagrams
- **Feature Selection**: Informing feature prioritization for machine learning models
- **Biological Interpretation**: Connecting imaging patterns to molecular mechanisms

## Clinical Interpretation

### Correlation Strength Guidelines:
- **|r| > 0.7**: Strong association
- **|r| 0.5-0.7**: Moderate association
- **|r| 0.3-0.5**: Weak to moderate association
- **|r| < 0.3**: Weak association

### Statistical Significance:
- **FDR p-value < 0.05**: Significant after multiple testing correction
- **FDR p-value < 0.10**: Suggestive association
- **Raw p-value < 0.05**: Nominally significant

### Biological Interpretation:
- **Positive Correlations**: Concordant radiomic and genomic patterns
- **Negative Correlations**: Inverse relationship between imaging and molecular features
- **Pathway-Specific Patterns**: Different pathway collections may show distinct correlation profiles

## Authors
- Jackie Chen
- Julia Nguyen

## Version
- Pipeline Version: 2025.1
- Created: June 2025
- Last Updated: July 2025

