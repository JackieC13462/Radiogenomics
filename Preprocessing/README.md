# Preprocessing

## Overview

The Preprocessing directory contains essential data preparation and harmonization scripts for the radiogenomics pipeline. These scripts handle the complex task of integrating multi-modal datasets (genomic, radiomic, and clinical data) from different sources, ensuring data quality, consistency, and compatibility for downstream analysis.

## Purpose

The Preprocessing directory provides comprehensive data preparation tools for:
- **Data Harmonization**: Standardizing sample identifiers across different data modalities
- **Quality Control**: Filtering and cleaning datasets to ensure analytical robustness
- **Feature Engineering**: Processing and transforming features for optimal analysis
- **Clinical Data Integration**: Extracting and formatting clinical metadata for survival analysis

## Key Features

### 🔄 **Multi-Modal Data Integration**
- Seamless integration of genomic, radiomic, and clinical datasets
- Robust sample matching algorithms with fallback mechanisms
- Consistent identifier standardization across data types

### 🧹 **Data Quality Assurance**
- Comprehensive filtering for data completeness and quality
- Duplicate detection and resolution
- Missing data handling and imputation strategies

### ⚙️ **Feature Processing**
- Gene identifier conversion and standardization
- Correlation-based feature filtering
- Pathway-specific data organization

### 📊 **Clinical Data Management**
- Treatment-specific filtering for therapy outcomes
- Survival data extraction and formatting
- Clinical metadata harmonization

## Scripts Overview

### 🔗 **Sample Harmonization**

#### `unique_ID_generator.R`
**Primary Sample ID Harmonization Tool**
- Harmonizes sample identifiers between genomic and radiomic datasets
- Implements fallback mechanisms for different UID column formats
- Ensures consistent sample matching across data modalities
- **Enhanced Features**: 
  - Prioritizes `seg_series_UID` with fallback to `SeriesInstanceUID_mask`
  - Robust error handling for missing columns
  - Comprehensive logging of harmonization process

#### `Clinical_filter_uniqueID.R`
**Clinical Data Sample Filtering**
- Filters clinical data to match available multi-omics samples
- Maintains consistency between clinical and molecular datasets
- Handles complex sample identifier mappings

### 🧬 **Genomic Data Processing**

#### `geneID_converter.R`
**Gene Identifier Standardization**
- Converts Entrez gene IDs to HGNC gene symbols using BioMart
- Essential for pathway enrichment analysis compatibility
- Provides mapping tables and conversion summaries
- **Key Features**:
  - Internet-based BioMart database queries
  - Human genome annotation support
  - Comprehensive conversion reporting

#### `protein_encoding_filtering.R`
**Gene Type Filtering**
- Filters gene expression data to include only protein-coding genes
- Removes non-coding RNAs and pseudogenes for focused analysis
- Improves pathway enrichment specificity

#### `combining_counts.R`
**RNA-seq Count Matrix Assembly**
- Combines individual sample count files into unified matrices
- Handles TCGA and other repository data structures
- Maps file names to sample identifiers using metadata
- **Processing Features**:
  - Duplicate sample detection and handling
  - Metadata-based file-to-sample mapping
  - Quality control logging

### 🏥 **Clinical Data Processing**

#### `clinical_data_extraction.R`
**Clinical Data Extraction and Filtering**
- Extracts clinical data for specific cancer types
- Filters for treatment types (Chemotherapy, Radiation Therapy)
- Matches clinical data to available molecular samples
- **Key Capabilities**:
  - Treatment-specific filtering
  - Cancer type prefix handling
  - Comprehensive extraction logging

#### `clinical_data_sampler.R`
**Clinical Data Sampling**
- Samples clinical data based on specific criteria
- Balances cohorts for statistical analysis
- Ensures representative sample selection

#### `clinical_duplicate_summary.R`
**Clinical Data Quality Assessment**
- Identifies and summarizes duplicate entries in clinical data
- Provides quality control metrics
- Guides data cleaning decisions

#### `clinical_outcome_intersect.R`
**Clinical Outcome Data Integration**
- Intersects clinical outcome data with molecular datasets
- Ensures survival data availability for analysis
- Handles missing outcome data appropriately

#### `TCGA_clinical_outcomes_sampling.R`
**TCGA-Specific Clinical Processing**
- Specialized processing for TCGA clinical outcome data
- Handles TCGA-specific data structures and identifiers
- Integrates with TCGA sample naming conventions

### 🔧 **Feature Engineering**

#### `feature_filtering.R`
**Advanced Feature Filtering System**
- Removes highly correlated features within pathway groups
- Separate filtering for genomic and radiomic features
- Pathway-specific correlation threshold application
- **Advanced Features**:
  - Multi-pathway correlation analysis (KEGG, HALLMARK, REACTOME, BIOCARTA)
  - Dynamic threshold setting based on data characteristics
  - Comprehensive filtering reports and statistics

#### `enrichment_separator.R`
**Pathway Enrichment Data Organization**
- Separates enrichment results by pathway database
- Organizes data for pathway-specific analysis
- Maintains traceability across enrichment methods

#### `Data_sampler.R`
**General Data Sampling Utility**
- Provides flexible data sampling capabilities
- Supports stratified and random sampling strategies
- Maintains data integrity during sampling

## Data Flow Architecture

### 1. **Raw Data Processing**
```
Raw Genomic Data → combining_counts.R → Unified Count Matrices
Raw Clinical Data → clinical_data_extraction.R → Filtered Clinical Data
Raw Radiomic Data → Direct Processing → Feature Matrices
```

### 2. **Identifier Harmonization**
```
Multi-Modal Data → unique_ID_generator.R → Harmonized Sample IDs
Clinical Data → Clinical_filter_uniqueID.R → Matched Clinical Data
```

### 3. **Feature/Signature Preparation**
```
Gene Expression → geneID_converter.R → Symbol-Based Expression
Protein Coding → protein_encoding_filtering.R → Filtered Genes to only protein coding genes (<50% missing and >80% of genes in the signature)
Features → feature_filtering.R → Correlation-Filtered Features
```

### 4. **Quality Control**
```
All Data → Duplicate Detection → Replace duplicate sample IDs with unique identifiers
```

## Input Data Requirements

### **Genomic Data**
- **Format**: TSV/CSV files with genes as rows, samples as columns
- **Identifiers**: Entrez IDs or Gene Symbols
- **Structure**: Count matrices or normalized expression data

### **Radiomic Data**
- **Format**: CSV files with samples as rows, features as columns
- **Identifiers**: Series instance UIDs or other imaging identifiers
- **Content**: Quantitative imaging features (shape, texture, intensity)

### **Clinical Data**
- **Format**: CSV/TSV files with comprehensive clinical variables
- **Required Fields**: Sample IDs, survival data, treatment information
- **Structure**: Patient-level clinical and demographic data

## Output Standards

### **Harmonized Datasets**
- Consistent sample identifier formatting across all data types
- Matched samples only (intersection of available data)
- Standardized file naming conventions

### **Quality Control Reports**
- Sample matching statistics and success rates
- Feature filtering summaries and retained features
- Data completeness and quality metrics

### **Processing Logs**
- Comprehensive logging of all preprocessing steps
- Error handling and resolution documentation
- Traceability of data transformations

## Integration Points

### **Upstream Integration**
- **Data Sources**: TCGA, CPTAC, institutional databases
- **File Formats**: Standard genomic and clinical data formats
- **Quality Requirements**: Minimum sample sizes and data completeness

### **Downstream Integration**
- **[`Data_Analysis/`](../Data_Analysis/)**: Provides clean, harmonized data for analysis
- **[`Snakemake/`](../Snakemake/)**: Automated workflow execution
- **[`Visualization_scripts/`](../Visualization_scripts/)**: Processed data for visualization

## Dependencies

### R Packages
- **Data Processing**: `data.table`, `dplyr`, `tidyr`
- **Bioinformatics**: `biomaRt`, `org.Hs.eg.db`
- **File I/O**: `readr`, `readxl`
- **String Processing**: `stringr`, `stringi`

### External Resources
- **BioMart Database**: For gene identifier conversion
- **Internet Connection**: Required for biomaRt queries
- **System Resources**: Sufficient memory for large dataset processing

## Usage Guidelines

### **Getting Started**
1. Ensure all dependencies are installed
2. Configure file paths in script headers
3. Run scripts in logical sequence (harmonization → filtering → feature processing)
4. Monitor logs for quality control feedback

### **Configuration**
- Modify file paths in script headers as needed
- Adjust filtering thresholds based on dataset characteristics
- Configure cancer type prefixes for multi-cohort analysis

### **Troubleshooting**
- Check sample identifier formats for matching issues
- Verify internet connectivity for BioMart queries
- Monitor memory usage for large dataset processing
- Review logs for detailed error information

## Best Practices

### **Data Management**
- Maintain original data backups before preprocessing
- Document all preprocessing parameters and decisions
- Use version control for preprocessing script modifications

### **Quality Control**
- Regularly review processing logs and summaries
- Validate sample matching across modalities
- Monitor feature retention after filtering steps

### **Performance Optimization**
- Use data.table for large dataset processing
- Implement parallel processing where appropriate
- Monitor memory usage and optimize accordingly

---

**Note**: The Preprocessing directory is a collection of foundational scripts that are essential to the entire radiogenomics pipeline. These scripts are not an automated process, and can be selectively executed, however, refer to each pipeline's input requirements to ensure proper filtering of the raw data is done before running any of the pipelines.
