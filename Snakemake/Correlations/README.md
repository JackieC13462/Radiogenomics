# Radiogenomics Snakemake Pipelines

This directory contains Snakemake workflows for automated radiogenomic correlation analysis.

## Files
- `correlative_snakefile.snakefile` - Main correlation analysis workflow
- `correlative_config.yaml` - Configuration file specifying datasets and file paths
- `correlative_snakefile.sh` - Bash script for pipeline execution
- Other specialized workflows for enrichment analysis

## Correlative Analysis Pipeline

The main correlative pipeline (`correlative_snakefile.snakefile`) performs integrated radiogenomic correlation analysis with the following steps:

1. **Sample ID Harmonization**: Aligns sample identifiers between genomic and radiomic datasets
2. **Pathway Separation**: Separates genomic data by pathway source (KEGG, HALLMARK, REACTOME, BIOCARTA)
3. **Feature Filtering**: Removes highly correlated features to reduce redundancy
4. **Radiogenomic Correlation**: Performs correlation analysis between imaging features and genomic pathway signatures
5. **Clinical Correlation**: Analyzes correlations between signatures and clinical outcomes
6. **Correlation Filtering**: Applies user-defined thresholds to filter significant correlations

### Key Features

- **Two-step correlation analysis**: Correlation calculation and filtering are separated for flexibility
- **Multi-pathway support**: Handles KEGG, HALLMARK, REACTOME, and BIOCARTA pathways
- **Configurable thresholds**: Filtering criteria can be adjusted in the workflow
- **Scalable**: Designed to process multiple cancer cohorts in parallel

### Usage

1. Update `correlative_config.yaml` with your dataset paths
2. Run the pipeline: `bash correlative_snakefile.sh`

### Output Files

- `{cohort}_{pathway}_clinical_correlations.csv` - Full correlation results
- `{cohort}_{pathway}_clinical_correlations_filtered.csv` - Filtered significant correlations
- `{cohort}_{pathway}_correlative_analysis.csv` - Radiogenomic correlation results
