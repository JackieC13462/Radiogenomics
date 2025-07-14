# Enrichment Analysis Pipeline (Snakemake)

## Overview

This Snakemake pipeline automates gene set enrichment analysis (GSEA) for radiogenomics studies across multiple cancer types. It performs pathway enrichment analysis using four major gene set collections from the Molecular Signatures Database (MSigDB) to calculate sample-wise enrichment scores that can be correlated with radiomic features.

## Pipeline Workflow

The pipeline consists of the following sequential steps:

1. **Patient Filtering** - Filters RNA-seq data to include only patients with both genomic and radiomic data
2. **Pathway Enrichment Analysis** - Calculates enrichment scores using four gene set collections:
   - **Hallmark Pathways** - 50 well-characterized biological processes
   - **KEGG Pathways** - Metabolic and signaling pathways from KEGG database  
   - **Reactome Pathways** - Comprehensive pathway database with detailed biological processes
   - **BioCarta Pathways** - Cancer-focused pathway collection
3. **Score Compilation** - Combines all enrichment scores into a single file per cancer type

## File Structure

```
Snakemake/Enrichment/
├── Snakefile_Enrichment.snakefile    # Main pipeline definition
├── enrichment_config.yaml            # Configuration file with paths and samples
├── Snakemake_enrichment.sh          # SLURM submission script
└── README.md                         # This documentation
```

## Prerequisites

### Software Requirements
- **Snakemake** (workflow management)
- **R** with the following packages:
  - `data.table` - Fast data manipulation
  - `GSEABase` - Gene set operations
  - `GSVA` - Gene Set Variation Analysis

### Data Requirements

#### Input Files (specified in config):
1. **RNA-seq Expression Data**: TPM-normalized gene expression matrices
   - Format: CSV with genes as rows, samples as columns
   - First column: Gene symbols
   - Path example: `/cluster/projects/bhklab/rawdata/Radiogenomics/TPM_RNAseq_data/{cancer}_filtered_TPM.csv`

2. **Patient ID Lists**: Text files containing patient IDs to include in analysis
   - Format: One patient ID per line
   - Path example: `/cluster/projects/bhklab/rawdata/Radiogenomics/sample_ids/{cancer}_patientids.txt`

3. **GMT Files**: Gene set definitions in GMT format
   - Hallmark: `h.all.v2024.1.Hs.symbols.gmt`
   - KEGG: `c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt`
   - Reactome: `c2.cp.reactome.v2024.1.Hs.symbols.gmt`
   - BioCarta: `c2.cp.biocarta.v2024.1.Hs.symbols.gmt`

#### R Scripts (must be available in `/cluster/home/t138199uhn/Scripts/`):
- `Data_sampler.R` - Patient filtering
- `hallmark_enrichment_GMT.R` - Hallmark pathway analysis
- `kegg_enrichment_GMT.R` - KEGG pathway analysis  
- `reactome_enrichment_GMT.R` - Reactome pathway analysis
- `biocarta_enrichment_GMT.R` - BioCarta pathway analysis
- `Enrichment_compiler.R` - Score compilation

## Configuration

### `enrichment_config.yaml`

This file contains all pipeline parameters:

```yaml
# Cancer types to analyze
SAMPLES: GBM,BRCA,KIRC,LGG,HNSCC,PDA,CCRCC

# RNA-seq expression data paths
datasets:
  CANCER_TYPE: "/path/to/expression/data.csv"

# Patient ID file paths  
patient_ids:
  CANCER_TYPE: "/path/to/patient/ids.txt"

# GMT file paths
hallmark_gmt: "/path/to/hallmark.gmt"
kegg_gmt: "/path/to/kegg.gmt"
reactome_gmt: "/path/to/reactome.gmt"
biocarta_gmt: "/path/to/biocarta.gmt"
```

**Supported Cancer Types:**
- **GBM** - Glioblastoma Multiforme (TCGA)
- **BRCA** - Breast Invasive Carcinoma (TCGA)
- **KIRC** - Kidney Renal Clear Cell Carcinoma (TCGA)
- **LGG** - Lower Grade Glioma (TCGA)
- **HNSCC** - Head and Neck Squamous Cell Carcinoma (CPTAC)
- **PDA** - Pancreatic Ductal Adenocarcinoma (CPTAC)
- **CCRCC** - Clear Cell Renal Cell Carcinoma (CPTAC)

## Usage

### 1. Update Configuration
Edit `enrichment_config.yaml` to specify:
- Cancer types to analyze (SAMPLES)
- Input file paths for your data
- GMT file locations

### 2. Submit Pipeline

#### Option A: Direct Snakemake Execution
```bash
cd /path/to/snakemake/enrichment
snakemake --snakefile Snakefile_Enrichment.snakefile \
  --configfile enrichment_config.yaml \
  --cores 8
```

#### Option B: SLURM Cluster Submission (Recommended)
```bash
sbatch Snakemake_enrichment.sh
```

### 3. Monitor Progress
- Check SLURM log files in `/cluster/home/t138199uhn/slurm/slurm_logs/GSVA_run/`
- Monitor job status: `squeue -u $USER`

## Output Files

### Intermediate Outputs

#### Filtered Data:
- **Location**: `/cluster/projects/bhklab/procdata/Radiogenomics/outputs/filtered_RNAseq/`
- **Format**: `{cancer}_filtered.csv`
- **Content**: Expression data filtered to patients with radiomic data

#### Individual Enrichment Scores:
- **Hallmark**: `/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Hallmark_scores/{cancer}_hallmark_signatures.csv`
- **KEGG**: `/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Kegg_scores/{cancer}_kegg_enrichment.csv`
- **Reactome**: `/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Reactome_scores/{cancer}_reactome_enrichment.csv`
- **BioCarta**: `/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Biocarta_scores/{cancer}_biocarta_enrichment.csv`

### Final Output

#### Compiled Enrichment Scores:
- **Location**: `/cluster/projects/bhklab/procdata/Radiogenomics/gene_signatures/`
- **Format**: `{cancer}_compiled_enrichment.csv`
- **Content**: Combined enrichment scores from all four pathway collections
- **Structure**: 
  - Rows: All pathways from all collections
  - Columns: Patient samples
  - Values: GSVA enrichment scores (-1 to +1 range)

## Resource Requirements

### Computational Resources per Rule:

| Rule | Memory | Runtime | Purpose |
|------|--------|---------|---------|
| filter_patients | 2GB | 30 min | Patient filtering |
| hallmark_signature_extraction | 8GB | 2 hours | Hallmark analysis |
| kegg_gene_enrichment | 8GB | 2 hours | KEGG analysis |
| reactome_gene_enrichment | 12GB | 14 hours | Reactome analysis (largest) |
| biocarta_gene_enrichment | 10GB | 2 hours | BioCarta analysis |
| compile_enrichment_scores | 4GB | 30 min | Score compilation |

**Note**: Reactome analysis requires the most resources due to the large number of pathways (~2000).

## Methodology

### Gene Set Variation Analysis (GSVA)
- **Algorithm**: Kernel density estimation-based enrichment scoring
- **Output**: Single-sample enrichment scores (no phenotype comparison required)
- **Range**: Scores normalized to approximately -1 to +1
- **Interpretation**: 
  - Positive scores = pathway activation
  - Negative scores = pathway suppression
  - Magnitude indicates strength of enrichment

### Quality Control
- Gene symbol matching between expression data and GMT files
- Patient ID consistency across genomic and radiomic datasets
- Pathway coverage validation (minimum gene overlap requirements)

## Troubleshooting

### Common Issues:

1. **Missing Input Files**
   - Verify all paths in `enrichment_config.yaml`
   - Check file permissions and accessibility

2. **Memory Errors**
   - Increase memory allocation in Snakefile resources
   - Consider running fewer samples in parallel

3. **Gene Symbol Mismatches**
   - Ensure expression data uses HUGO gene symbols
   - Update GMT files to match gene symbol versions

4. **SLURM Job Failures**
   - Check individual job logs in slurm_logs directory
   - Verify module loading (R, snakemake)
   - Confirm script paths exist

### Log Files:
- **Master job**: `snakemakemaster_{jobid}.out/err`
- **Individual rules**: `{rule}.{cancer}.out/err`

## Integration with Radiogenomics Pipeline

The compiled enrichment scores serve as input for:
- **Correlation Analysis**: Correlating pathway scores with radiomic features
- **Machine Learning**: Using pathway scores as genomic features for predictive modeling
- **Survival Analysis**: Analyzing pathway associations with clinical outcomes
- **Visualization**: Creating pathway-radiomic correlation heatmaps

## Authors
- Jackie Chen
- Julia Nguyen

## Version
- Pipeline Version: 2025.1
- MSigDB Version:  May 2024 v2024.1.Hs (https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)
- Created: May 2025
- Last Updated: July 2025
