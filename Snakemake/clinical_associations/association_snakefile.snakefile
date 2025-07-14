# Load configuration
configfile: "clinical_config.yaml"

# Harmonize pathway group names to uppercase for consistency
cancer_types = config["cancer_types"]
pathway_groups = ["HALLMARK", "KEGG", "REACTOME", "BIOCARTA"]
OUTDIR = "/cluster/projects/bhklab/procdata/Radiogenomics/clinical/{cancer}/"

# Add rule all to expand all expected outputs for all cancer types and pathway groups, including both binary and continuous Cox model outputs and harmonized clinical data
rule all:
    input:
        # Harmonized clinical data
        expand(OUTDIR + "{cancer}_harmonized_clinical.csv", cancer=cancer_types) +
        # Genomic Cox model outputs (binary and continuous for each pathway)
        expand(OUTDIR + "{cancer}_HALLMARK_cox_results_binary.csv", cancer=cancer_types) +
        expand(OUTDIR + "{cancer}_HALLMARK_cox_results_continuous.csv", cancer=cancer_types) +
        expand(OUTDIR + "{cancer}_KEGG_cox_results_binary.csv", cancer=cancer_types) +
        expand(OUTDIR + "{cancer}_KEGG_cox_results_continuous.csv", cancer=cancer_types) +
        expand(OUTDIR + "{cancer}_REACTOME_cox_results_binary.csv", cancer=cancer_types) +
        expand(OUTDIR + "{cancer}_REACTOME_cox_results_continuous.csv", cancer=cancer_types) +
        expand(OUTDIR + "{cancer}_BIOCARTA_cox_results_binary.csv", cancer=cancer_types) +
        expand(OUTDIR + "{cancer}_BIOCARTA_cox_results_continuous.csv", cancer=cancer_types) +
        # Radiomics Cox model outputs (binary and continuous)
        expand(OUTDIR + "{cancer}_radiomics_cox_results_binary.csv", cancer=cancer_types) +
        expand(OUTDIR + "{cancer}_radiomics_cox_results_continuous.csv", cancer=cancer_types)

# Rule: Clinical Data Extraction
rule clinical_data_extraction:
    input:
        dataset_files=lambda wildcards: config["clinical_data"][wildcards.cancer],
        sample_id_files=lambda wildcards: config["sample_ids"][wildcards.cancer]
    output:
        output_file=OUTDIR + "{cancer}_extracted_clinical.csv"
    params:
        cancer="{cancer}"
    resources:
        mem_mb=8000,
        time="01:00:00"
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/clinical_data_extraction.R \
            {input.dataset_files} {input.sample_id_files} {output.output_file} {params.cancer}
        """
# Rule: Clinical Filter Unique ID
rule clinical_filter_unique_id:
    input:
        clinical_file=OUTDIR + "{cancer}_extracted_clinical.csv",
        radiomics_file=lambda wildcards: config["radiomics_file"][wildcards.cancer],
        hallmark_file=lambda wildcards: config["genomics_file"][wildcards.cancer]["HALLMARK"],
        kegg_file=lambda wildcards: config["genomics_file"][wildcards.cancer]["KEGG"],
        reactome_file=lambda wildcards: config["genomics_file"][wildcards.cancer]["REACTOME"],
        biocarta_file=lambda wildcards: config["genomics_file"][wildcards.cancer]["BIOCARTA"]
    output:
        harmonized_clinical=OUTDIR + "{cancer}_harmonized_clinical.csv",
        harmonized_radiomics=OUTDIR + "{cancer}_harmonized_radiomics.csv",
        harmonized_hallmark=OUTDIR + "{cancer}_harmonized_hallmark_genomics.csv",
        harmonized_kegg=OUTDIR + "{cancer}_harmonized_kegg_genomics.csv",
        harmonized_reactome=OUTDIR + "{cancer}_harmonized_reactome_genomics.csv",
        harmonized_biocarta=OUTDIR + "{cancer}_harmonized_biocarta_genomics.csv"
    params:
        output_prefix="{cancer}",
        output_dir="/cluster/projects/bhklab/procdata/Radiogenomics/clinical/{cancer}"
    resources:
        mem_mb=4000,
        time="00:30:00"
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/Clinical_filter_uniqueID.R \
            {input.clinical_file} {input.radiomics_file} \
            {input.hallmark_file} {input.kegg_file} {input.reactome_file} {input.biocarta_file} \
            {params.output_prefix} {params.output_dir}
        """

# Rule: Genomic Cox Model for HALLMARK
rule genomic_cox_hallmark:
    input:
        harmonized_clinical=OUTDIR + "{cancer}_harmonized_clinical.csv",
        harmonized_genomics=OUTDIR + "{cancer}_harmonized_hallmark_genomics.csv"
    output:
        binary=OUTDIR + "{cancer}_HALLMARK_cox_results_binary.csv",
        continuous=OUTDIR + "{cancer}_HALLMARK_cox_results_continuous.csv"
    params:
        column_prefix="{cancer}",
        output_prefix="{cancer}_HALLMARK_cox_results",
        output_dir="/cluster/projects/bhklab/procdata/Radiogenomics/clinical/{cancer}"
    resources:
        mem_mb=16000,
        time="14:00:00"
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/Genomic_cox_model.R \
            {input.harmonized_genomics} {input.harmonized_clinical} {params.output_dir} {params.column_prefix} {params.output_prefix}
        """

# Rule: Genomic Cox Model for KEGG
rule genomic_cox_kegg:
    input:
        harmonized_clinical=OUTDIR + "{cancer}_harmonized_clinical.csv",
        harmonized_genomics=OUTDIR + "{cancer}_harmonized_kegg_genomics.csv"
    output:
        binary=OUTDIR + "{cancer}_KEGG_cox_results_binary.csv",
        continuous=OUTDIR + "{cancer}_KEGG_cox_results_continuous.csv"
    params:
        column_prefix="{cancer}",
        output_prefix="{cancer}_KEGG_cox_results",
        output_dir="/cluster/projects/bhklab/procdata/Radiogenomics/clinical/{cancer}"
    resources:
        mem_mb=16000,
        time="14:00:00"
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/Genomic_cox_model.R \
            {input.harmonized_genomics} {input.harmonized_clinical} {params.output_dir} {params.column_prefix} {params.output_prefix}
        """

# Rule: Genomic Cox Model for REACTOME
rule genomic_cox_reactome:
    input:
        harmonized_clinical=OUTDIR + "{cancer}_harmonized_clinical.csv",
        harmonized_genomics=OUTDIR + "{cancer}_harmonized_reactome_genomics.csv"
    output:
        binary=OUTDIR + "{cancer}_REACTOME_cox_results_binary.csv",
        continuous=OUTDIR + "{cancer}_REACTOME_cox_results_continuous.csv"
    params:
        column_prefix="{cancer}",
        output_prefix="{cancer}_REACTOME_cox_results",
        output_dir="/cluster/projects/bhklab/procdata/Radiogenomics/clinical/{cancer}"
    resources:
        mem_mb=16000,
        time="14:00:00"
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/Genomic_cox_model.R \
            {input.harmonized_genomics} {input.harmonized_clinical} {params.output_dir} {params.column_prefix} {params.output_prefix}
        """

# Rule: Genomic Cox Model for BIOCARTA
rule genomic_cox_biocarta:
    input:
        harmonized_clinical=OUTDIR + "{cancer}_harmonized_clinical.csv",
        harmonized_genomics=OUTDIR + "{cancer}_harmonized_biocarta_genomics.csv"
    output:
        binary=OUTDIR + "{cancer}_BIOCARTA_cox_results_binary.csv",
        continuous=OUTDIR + "{cancer}_BIOCARTA_cox_results_continuous.csv"
    params:
        column_prefix="{cancer}",
        output_prefix="{cancer}_BIOCARTA_cox_results",
        output_dir="/cluster/projects/bhklab/procdata/Radiogenomics/clinical/{cancer}"
    resources:
        mem_mb=16000,
        time="14:00:00"
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/Genomic_cox_model.R \
            {input.harmonized_genomics} {input.harmonized_clinical} {params.output_dir} {params.column_prefix} {params.output_prefix}
        """

# Rule: Radiomics Cox Model
rule radiomics_cox_model:
    input:
        harmonized_clinical=OUTDIR + "{cancer}_harmonized_clinical.csv",
        harmonized_radiomics=OUTDIR + "{cancer}_harmonized_radiomics.csv"
    output:
        binary=OUTDIR + "{cancer}_radiomics_cox_results_binary.csv",
        continuous=OUTDIR + "{cancer}_radiomics_cox_results_continuous.csv"
    params:
        column_prefix="{cancer}",
        output_prefix="{cancer}_radiomics_cox_results",
        output_dir="/cluster/projects/bhklab/procdata/Radiogenomics/clinical/{cancer}"
    resources:
        mem_mb=16000,
        time="14:30:00"
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/Radiomics_cox_model.R \
            {input.harmonized_radiomics} {input.harmonized_clinical} {params.output_dir} {params.column_prefix} {params.output_prefix}
        """
