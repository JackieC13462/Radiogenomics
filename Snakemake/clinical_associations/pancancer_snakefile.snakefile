# Load configuration
configfile: "pancancer_config.yaml"

# Harmonize pathway group names to uppercase for consistency
cancer_types = config["cancer_types"]
pathway_groups = ["HALLMARK", "KEGG", "REACTOME", "BIOCARTA"]
# Use the first cancer type's output directory since all point to the same PANCAN directory
PANCAN_OUTDIR = config["output_dir"][cancer_types[0]]

# Add rule all to expand all expected outputs for all cancer types and pathway groups, including both binary and continuous Cox model outputs and harmonized clinical data
rule all:
    input:
        # Harmonized clinical data
        [config["output_dir"][cancer] + f"{cancer}_harmonized_clinical.csv" for cancer in cancer_types] +
        # Combined datasets
        [PANCAN_OUTDIR + "COMBINED_harmonized_clinical.csv",
         PANCAN_OUTDIR + "COMBINED_harmonized_kegg_genomics.csv",
         PANCAN_OUTDIR + "COMBINED_harmonized_hallmark_genomics.csv",
         PANCAN_OUTDIR + "COMBINED_harmonized_reactome_genomics.csv",
         PANCAN_OUTDIR + "COMBINED_harmonized_biocarta_genomics.csv"] +
        # Pan-cancer Cox model outputs
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
        output_file=lambda wildcards: config["output_dir"][wildcards.cancer] + f"{wildcards.cancer}_extracted_clinical.csv"
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
        clinical_file=lambda wildcards: config["output_dir"][wildcards.cancer] + f"{wildcards.cancer}_extracted_clinical.csv",
        hallmark_file=lambda wildcards: config["genomics_file"][wildcards.cancer]["HALLMARK"],
        kegg_file=lambda wildcards: config["genomics_file"][wildcards.cancer]["KEGG"],
        reactome_file=lambda wildcards: config["genomics_file"][wildcards.cancer]["REACTOME"],
        biocarta_file=lambda wildcards: config["genomics_file"][wildcards.cancer]["BIOCARTA"]
    output:
        harmonized_clinical=lambda wildcards: config["output_dir"][wildcards.cancer] + f"{wildcards.cancer}_harmonized_clinical.csv",
        harmonized_hallmark=lambda wildcards: config["output_dir"][wildcards.cancer] + f"{wildcards.cancer}_harmonized_hallmark_genomics.csv",
        harmonized_kegg=lambda wildcards: config["output_dir"][wildcards.cancer] + f"{wildcards.cancer}_harmonized_kegg_genomics.csv",
        harmonized_reactome=lambda wildcards: config["output_dir"][wildcards.cancer] + f"{wildcards.cancer}_harmonized_reactome_genomics.csv",
        harmonized_biocarta=lambda wildcards: config["output_dir"][wildcards.cancer] + f"{wildcards.cancer}_harmonized_biocarta_genomics.csv"
    params:
        output_dir=lambda wildcards: config["output_dir"][wildcards.cancer]
    resources:
        mem_mb=4000,
        time="00:30:00"
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/gen_clin_newID.R \
            {input.clinical_file} \
            {input.hallmark_file} {input.kegg_file} {input.reactome_file} {input.biocarta_file} \
            {params.output_dir}
        """
rule combine_datasets:
    input:
        # All harmonized files that need to be combined
        clinical_files=[config["output_dir"][cancer] + f"{cancer}_harmonized_clinical.csv" for cancer in cancer_types],
        genomic_files=[config["output_dir"][cancer] + f"{cancer}_harmonized_{pathway}_genomics.csv" 
                      for cancer in cancer_types 
                      for pathway in ["kegg", "hallmark", "reactome", "biocarta"]]
    output:
        combined_clinical=PANCAN_OUTDIR + "COMBINED_harmonized_clinical.csv",
        combined_kegg=PANCAN_OUTDIR + "COMBINED_harmonized_kegg_genomics.csv",
        combined_hallmark=PANCAN_OUTDIR + "COMBINED_harmonized_hallmark_genomics.csv",
        combined_reactome=PANCAN_OUTDIR + "COMBINED_harmonized_reactome_genomics.csv",
        combined_biocarta=PANCAN_OUTDIR + "COMBINED_harmonized_biocarta_genomics.csv"
    params:
        input_dir="/cluster/projects/bhklab/procdata/Radiogenomics/clinical",
        output_dir=PANCAN_OUTDIR
    resources:
        mem_mb=4000,
        time="00:30:00"
    shell:
        """
        # Create output directory
        mkdir -p {params.output_dir}
        
        # Run the combine script
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/combine_datasets.R "{params.input_dir}" "{params.output_dir}"
        """

# Rule: Pan-cancer Cox modeling for KEGG pathways
rule pancancer_kegg_cox:
    input:
        genomic_file=PANCAN_OUTDIR + "COMBINED_harmonized_kegg_genomics.csv",
        clinical_file=PANCAN_OUTDIR + "COMBINED_harmonized_clinical.csv"
    output:
        cox_results=PANCAN_OUTDIR + "COMBINED_KEGG_pancancer_cox_results_binary.csv"
    params:
        output_dir=PANCAN_OUTDIR,
        output_prefix="COMBINED_KEGG_pancancer_cox_results"
    resources:
        mem_mb=8000,
        time="02:00:00"
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/pancancer_gen_cox.R \
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
        output_dir=PANCAN_OUTDIR,
        output_prefix="COMBINED_HALLMARK_pancancer_cox_results"
    resources:
        mem_mb=8000,
        time="02:00:00"
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/pancancer_gen_cox.R \
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
        output_dir=PANCAN_OUTDIR,
        output_prefix="COMBINED_REACTOME_pancancer_cox_results"
    resources:
        mem_mb=8000,
        time="02:00:00"
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/pancancer_gen_cox.R \
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
        output_dir=PANCAN_OUTDIR,
        output_prefix="COMBINED_BIOCARTA_pancancer_cox_results"
    resources:
        mem_mb=8000,
        time="02:00:00"
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/pancancer_gen_cox.R \
            {input.genomic_file} {input.clinical_file} {params.output_dir} {params.output_prefix}
        """
