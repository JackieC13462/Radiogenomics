# Load configuration
configfile: "pancancer_config.yaml"

# Define script directory for R scripts
SCRIPT_DIR = "/cluster/home/t138199uhn/Scripts"

# Harmonize pathway group names to uppercase for consistency
cancer_types = config["cancer_types"]
pathway_groups = ["HALLMARK", "KEGG", "REACTOME", "BIOCARTA"]
# Define output directory structure
OUTDIR = "/cluster/projects/bhklab/procdata/Radiogenomics/clinical/PANCAN/{cancer}/"
# Extract the parent PANCAN directory for combined outputs
PANCAN_OUTDIR = "/cluster/projects/bhklab/procdata/Radiogenomics/clinical/PANCAN/"

# Add rule all to expand all expected outputs for all cancer types and pathway groups
rule all:
    input:
        # Individual cancer type harmonized files (in cancer-specific output directories)
        expand(OUTDIR + "{cancer}_harmonized_clinical.csv", cancer=cancer_types) +
        expand(OUTDIR + "{cancer}_harmonized_hallmark_genomics.csv", cancer=cancer_types) +
        expand(OUTDIR + "{cancer}_harmonized_kegg_genomics.csv", cancer=cancer_types) +
        expand(OUTDIR + "{cancer}_harmonized_reactome_genomics.csv", cancer=cancer_types) +
        expand(OUTDIR + "{cancer}_harmonized_biocarta_genomics.csv", cancer=cancer_types) +
        # Combined pan-cancer datasets (in main PANCAN directory)
        [PANCAN_OUTDIR + "COMBINED_harmonized_clinical.csv",
         PANCAN_OUTDIR + "COMBINED_harmonized_kegg_genomics.csv",
         PANCAN_OUTDIR + "COMBINED_harmonized_hallmark_genomics.csv",
         PANCAN_OUTDIR + "COMBINED_harmonized_reactome_genomics.csv",
         PANCAN_OUTDIR + "COMBINED_harmonized_biocarta_genomics.csv"] +
        # Pan-cancer Cox model outputs (in main PANCAN directory)
        [PANCAN_OUTDIR + "COMBINED_KEGG_pancancer_cox_results_binary.csv",
         PANCAN_OUTDIR + "COMBINED_HALLMARK_pancancer_cox_results_binary.csv",
         PANCAN_OUTDIR + "COMBINED_REACTOME_pancancer_cox_results_binary.csv",
         PANCAN_OUTDIR + "COMBINED_BIOCARTA_pancancer_cox_results_binary.csv"]

# Rule: Clinical Data Extraction
rule clinical_data_extraction:
    input:
        dataset_files=lambda wildcards: config["clinical_data"][wildcards.cancer],
        sample_id_files=lambda wildcards: config["sample_ids"][wildcards.cancer]
    output:
        output_file=OUTDIR + "{cancer}_extracted_clinical.csv"
    params:
        cancer="{cancer}",
        output_dir=OUTDIR
    resources:
        mem_mb=8000,
        time="01:00:00"
    shell:
        """
        module load R
        Rscript {SCRIPT_DIR}/clinical_data_extraction.R \
            {input.dataset_files} {input.sample_id_files} {output.output_file} {params.cancer}
        """
# Rule: Clinical Filter Unique ID
rule clinical_filter_unique_id:
    input:
        clinical_file=OUTDIR + "{cancer}_extracted_clinical.csv",
        compiled_genomics_file=lambda wildcards: config["genomics_file"][wildcards.cancer]
    output:
        harmonized_clinical=OUTDIR + "{cancer}_harmonized_clinical.csv",
        harmonized_hallmark=OUTDIR + "{cancer}_harmonized_hallmark_genomics.csv",
        harmonized_kegg=OUTDIR + "{cancer}_harmonized_kegg_genomics.csv",
        harmonized_reactome=OUTDIR + "{cancer}_harmonized_reactome_genomics.csv",
        harmonized_biocarta=OUTDIR + "{cancer}_harmonized_biocarta_genomics.csv"
    params:
        output_dir=OUTDIR
    resources:
        mem_mb=4000,
        time="00:30:00"
    shell:
        """
        module load R
        Rscript {SCRIPT_DIR}/gen_clin_newID.R \
            {input.clinical_file} \
            {input.compiled_genomics_file} \
            $(echo {params.output_dir} | sed "s/{{cancer}}/{wildcards.cancer}/g")
        """
rule combine_datasets:
    input:
        # All harmonized files that need to be combined (from cancer-specific output directories)
        clinical_files=expand(OUTDIR + "{cancer}_harmonized_clinical.csv", cancer=cancer_types),
        kegg_files=expand(OUTDIR + "{cancer}_harmonized_kegg_genomics.csv", cancer=cancer_types),
        hallmark_files=expand(OUTDIR + "{cancer}_harmonized_hallmark_genomics.csv", cancer=cancer_types),
        reactome_files=expand(OUTDIR + "{cancer}_harmonized_reactome_genomics.csv", cancer=cancer_types),
        biocarta_files=expand(OUTDIR + "{cancer}_harmonized_biocarta_genomics.csv", cancer=cancer_types)
    output:
        combined_clinical=PANCAN_OUTDIR + "COMBINED_harmonized_clinical.csv",
        combined_kegg=PANCAN_OUTDIR + "COMBINED_harmonized_kegg_genomics.csv",
        combined_hallmark=PANCAN_OUTDIR + "COMBINED_harmonized_hallmark_genomics.csv",
        combined_reactome=PANCAN_OUTDIR + "COMBINED_harmonized_reactome_genomics.csv",
        combined_biocarta=PANCAN_OUTDIR + "COMBINED_harmonized_biocarta_genomics.csv"
    params:
        input_dir="/cluster/projects/bhklab/procdata/Radiogenomics/clinical/PANCAN",
        output_dir=PANCAN_OUTDIR
    resources:
        mem_mb=4000,
        time="00:30:00"
    shell:
        """
        # Create output directory
        mkdir -p {params.output_dir}
        module load R
        Rscript {SCRIPT_DIR}/combine_datasets.R "{params.input_dir}" "{params.output_dir}"
        """

# Rule: Pan-cancer Cox modeling for KEGG pathways
rule pancancer_kegg_cox:
    input:
        genomic_file=PANCAN_OUTDIR + "COMBINED_harmonized_kegg_genomics.csv",
        clinical_file=PANCAN_OUTDIR + "COMBINED_harmonized_clinical.csv"
    output:
        cox_results=PANCAN_OUTDIR + "COMBINED_KEGG_pancancer_cox_results_binary.csv"
    params:
        output_prefix="COMBINED_KEGG_pancancer_cox_results",
        output_dir=PANCAN_OUTDIR
    resources:
        mem_mb=8000,
        time="02:00:00"
    shell:
        """
        module load R
        Rscript {SCRIPT_DIR}/pancancer_gen_cox.R \
            {input.genomic_file} {input.clinical_file} {params.output_dir} {params.output_prefix}
        """

# Rule: Pan-cancer Cox modeling for Hallmark pathways
rule pancancer_hallmark_cox:
    input:
        genomic_file=PANCAN_OUTDIR + "COMBINED_harmonized_hallmark_genomics.csv",
        clinical_file=PANCAN_OUTDIR + "COMBINED_harmonized_clinical.csv"
    output:
        cox_results=PANCAN_OUTDIR + "COMBINED_HALLMARK_pancancer_cox_results_binary.csv"
    params:
        output_prefix="COMBINED_HALLMARK_pancancer_cox_results",
        output_dir=PANCAN_OUTDIR
    resources:
        mem_mb=8000,
        time="02:00:00"
    shell:
        """
        module load R
        Rscript {SCRIPT_DIR}/pancancer_gen_cox.R \
            {input.genomic_file} {input.clinical_file} {params.output_dir} {params.output_prefix}
        """

# Rule: Pan-cancer Cox modeling for Reactome pathways
rule pancancer_reactome_cox:
    input:
        genomic_file=PANCAN_OUTDIR + "COMBINED_harmonized_reactome_genomics.csv",
        clinical_file=PANCAN_OUTDIR + "COMBINED_harmonized_clinical.csv"
    output:
        cox_results=PANCAN_OUTDIR + "COMBINED_REACTOME_pancancer_cox_results_binary.csv"
    params:
        output_prefix="COMBINED_REACTOME_pancancer_cox_results",
        output_dir=PANCAN_OUTDIR
    resources:
        mem_mb=8000,
        time="02:00:00"
    shell:
        """
        module load R
        Rscript {SCRIPT_DIR}/pancancer_gen_cox.R \
            {input.genomic_file} {input.clinical_file} {params.output_dir} {params.output_prefix}
        """

# Rule: Pan-cancer Cox modeling for BioCarta pathways
rule pancancer_biocarta_cox:
    input:
        genomic_file=PANCAN_OUTDIR + "COMBINED_harmonized_biocarta_genomics.csv",
        clinical_file=PANCAN_OUTDIR + "COMBINED_harmonized_clinical.csv"
    output:
        cox_results=PANCAN_OUTDIR + "COMBINED_BIOCARTA_pancancer_cox_results_binary.csv"
    params:
        output_prefix="COMBINED_BIOCARTA_pancancer_cox_results",
        output_dir=PANCAN_OUTDIR
    resources:
        mem_mb=8000,
        time="02:00:00"
    shell:
        """
        module load R
        Rscript {SCRIPT_DIR}/pancancer_gen_cox.R \
            {input.genomic_file} {input.clinical_file} {params.output_dir} {params.output_prefix}
        """
