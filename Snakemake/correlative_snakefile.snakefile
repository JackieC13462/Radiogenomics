# Define config variables for input/output locations
configfile: "correlative_config.yaml"

COHORTS = [d["name"] for d in config["datasets"]]
PATHWAY_GROUPS = ["KEGG", "HALLMARK", "REACTOME", "BIOCARTA"]

def get_genomics_file(wildcards):
    return next(d["genomics_file"] for d in config["datasets"] if d["name"] == wildcards.prefix)

def get_radiomics_file(wildcards):
    return next(d["radiomics_file"] for d in config["datasets"] if d["name"] == wildcards.prefix)
def get_clinical_file(wildcards):
    return next(d["clinical_file"] for d in config["datasets"] if d["name"] == wildcards.prefix)

OUTDIR = "/cluster/projects/bhklab/procdata/Radiogenomics/outputs/correlations/{prefix}/"
RSCRIPTDIR = "/cluster/home/t138199uhn/Scripts/"

rule all:
    input:
        expand(OUTDIR + "{prefix}_{group}_correlative_analysis.csv", prefix=COHORTS, group=PATHWAY_GROUPS),
        expand(OUTDIR + "{prefix}_KEGG_clinical_correlations.csv", prefix=COHORTS),
        expand(OUTDIR + "{prefix}_KEGG_clinical_correlations_filtered.csv", prefix=COHORTS),
        expand(OUTDIR + "{prefix}_HALLMARK_clinical_correlations.csv", prefix=COHORTS),
        expand(OUTDIR + "{prefix}_HALLMARK_clinical_correlations_filtered.csv", prefix=COHORTS),
        expand(OUTDIR + "{prefix}_REACTOME_clinical_correlations.csv", prefix=COHORTS),
        expand(OUTDIR + "{prefix}_REACTOME_clinical_correlations_filtered.csv", prefix=COHORTS),
        expand(OUTDIR + "{prefix}_BIOCARTA_clinical_correlations.csv", prefix=COHORTS),
        expand(OUTDIR + "{prefix}_BIOCARTA_clinical_correlations_filtered.csv", prefix=COHORTS)

rule unique_ID_generator:
    input:
        genomics = get_genomics_file,
        radiomics = get_radiomics_file
    output:
        genomics_out = OUTDIR + "{prefix}_radiogenomic_RNAseq.csv",
        radiomics_out = OUTDIR + "{prefix}_radiogenomic_features.csv"
    params:
        prefix = "{prefix}",
        outdir = OUTDIR
    resources:
        mem_mb=8000,
        runtime="00:40:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p {params.outdir}
        Rscript {RSCRIPTDIR}unique_ID_generator.R {input.genomics} {input.radiomics} {params.prefix} {params.outdir}
        """
rule enrichment_separation:
    input:
        genomics = OUTDIR + "{prefix}_radiogenomic_RNAseq.csv",
    output:
        kegg = OUTDIR + "{prefix}_KEGG_enrichment.csv",
        hallmark = OUTDIR + "{prefix}_HALLMARK_enrichment.csv",
        reactome = OUTDIR + "{prefix}_REACTOME_enrichment.csv",
        biocarta = OUTDIR + "{prefix}_BIOCARTA_enrichment.csv"
    params:
        prefix = "{prefix}",
        outdir = OUTDIR
    resources:
        mem_mb=8000,
        runtime="00:30:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p {params.outdir}
        Rscript {RSCRIPTDIR}enrichment_separator.R {input.genomics} {params.outdir} {params.prefix}
        """
rule clinical_ID_filtering:
    input:
        kegg = OUTDIR + "{prefix}_KEGG_enrichment.csv",
        hallmark = OUTDIR + "{prefix}_HALLMARK_enrichment.csv",
        reactome = OUTDIR + "{prefix}_REACTOME_enrichment.csv",
        biocarta = OUTDIR + "{prefix}_BIOCARTA_enrichment.csv",
        radiomics = OUTDIR + "{prefix}_radiogenomic_features.csv",
        clinical = get_clinical_file
    output:
        kegg_ID = OUTDIR + "{prefix}_harmonized_kegg_genomics.csv",
        hallmark_ID = OUTDIR + "{prefix}_harmonized_hallmark_genomics.csv",
        reactome_ID = OUTDIR + "{prefix}_harmonized_reactome_genomics.csv",
        biocarta_ID = OUTDIR + "{prefix}_harmonized_biocarta_genomics.csv",
        radiomics_ID = OUTDIR + "{prefix}_harmonized_radiomics.csv",
        clinical_ID = OUTDIR + "{prefix}_harmonized_clinical.csv"
    params:
        prefix = "{prefix}",
        outdir = OUTDIR
    resources:
        mem_mb=8000,
        runtime="00:30:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p {params.outdir}
        Rscript {RSCRIPTDIR}Clinical_filter_uniqueID.R {input.clinical} {input.radiomics} {input.hallmark} {input.kegg} {input.reactome} {input.biocarta} {params.prefix} {params.outdir}
        """
rule radiomics_self_correlation:
    input:
        radiomics = OUTDIR + "{prefix}_harmonized_radiomics.csv"
    output:
        corr = OUTDIR + "{prefix}_radiomics_self_correlation.csv"
    resources:
        mem_mb=8000,
        runtime="2:00:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p $(dirname {output.corr})
        Rscript {RSCRIPTDIR}radiomics_self_correlation.R {input.radiomics} {output.corr}
        """

rule genomics_self_correlation_HALLMARK:
    input:
        genomics = OUTDIR + "{prefix}_harmonized_hallmark_genomics.csv"
    output:
        hallmark_cor = OUTDIR + "{prefix}_HALLMARK_self_correlation.csv"
    params:
        outdir = lambda wildcards: OUTDIR.format(prefix=wildcards.prefix),
        prefix = "{prefix}"
    resources:
        mem_mb=8000,
        runtime="2:00:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p {params.outdir}
        Rscript {RSCRIPTDIR}genomics_self_correlation.R {input.genomics} {params.outdir} {params.prefix}
        """
rule genomics_self_correlation_KEGG:
    input:
        genomics = OUTDIR + "{prefix}_harmonized_kegg_genomics.csv"
    output:
        kegg_cor = OUTDIR + "{prefix}_KEGG_self_correlation.csv"
    params:
        outdir = lambda wildcards: OUTDIR.format(prefix=wildcards.prefix),
        prefix = "{prefix}"
    resources:
        mem_mb=8000,
        runtime="2:00:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p {params.outdir}
        Rscript {RSCRIPTDIR}genomics_self_correlation.R {input.genomics} {params.outdir} {params.prefix}
        """
rule genomics_self_correlation_REACTOME:
    input:
        genomics = OUTDIR + "{prefix}_harmonized_reactome_genomics.csv"
    output:
        reactome_cor = OUTDIR + "{prefix}_REACTOME_self_correlation.csv"
    params:
        outdir = lambda wildcards: OUTDIR.format(prefix=wildcards.prefix),
        prefix = "{prefix}"
    resources:
        mem_mb=8000,
        runtime="2:00:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p {params.outdir}
        Rscript {RSCRIPTDIR}genomics_self_correlation.R {input.genomics} {params.outdir} {params.prefix}
        """
rule genomics_self_correlation_BIOCARTA:
    input:
        genomics = OUTDIR + "{prefix}_harmonized_biocarta_genomics.csv"
    output:
        biocarta_cor = OUTDIR + "{prefix}_BIOCARTA_self_correlation.csv"
    params:
        outdir = lambda wildcards: OUTDIR.format(prefix=wildcards.prefix),
        prefix = "{prefix}"
    resources:
        mem_mb=8000,
        runtime="2:00:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p {params.outdir}
        Rscript {RSCRIPTDIR}genomics_self_correlation.R {input.genomics} {params.outdir} {params.prefix}
        """
rule feature_filtering:
    input:
        kegg_cor = OUTDIR + "{prefix}_KEGG_self_correlation.csv",
        hallmark_cor = OUTDIR + "{prefix}_HALLMARK_self_correlation.csv",
        reactome_cor = OUTDIR + "{prefix}_REACTOME_self_correlation.csv",
        biocarta_cor = OUTDIR + "{prefix}_BIOCARTA_self_correlation.csv",
        radiomics_cor = OUTDIR + "{prefix}_radiomics_self_correlation.csv",
        kegg_data = OUTDIR + "{prefix}_harmonized_kegg_genomics.csv",
        hallmark_data = OUTDIR + "{prefix}_harmonized_hallmark_genomics.csv",
        reactome_data = OUTDIR + "{prefix}_harmonized_reactome_genomics.csv",
        biocarta_data = OUTDIR + "{prefix}_harmonized_biocarta_genomics.csv",
        radiomics_data = OUTDIR + "{prefix}_harmonized_radiomics.csv"
    output:
        kegg_filtered = OUTDIR + "{prefix}_KEGG_features_filtered.csv",
        hallmark_filtered = OUTDIR + "{prefix}_HALLMARK_features_filtered.csv",
        reactome_filtered = OUTDIR + "{prefix}_REACTOME_features_filtered.csv",
        biocarta_filtered = OUTDIR + "{prefix}_BIOCARTA_features_filtered.csv",
        radiomics_filtered = OUTDIR + "{prefix}_radiomics_features_filtered.csv"
    params:
        outdir = OUTDIR,
        prefix = "{prefix}"
    resources:
        mem_mb=4000,
        runtime="00:50:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p {params.outdir}
        Rscript {RSCRIPTDIR}feature_filtering.R \
            {input.kegg_cor} {input.hallmark_cor} {input.reactome_cor} {input.biocarta_cor} \
            {input.kegg_data} {input.hallmark_data} {input.reactome_data} {input.biocarta_data} \
            {input.radiomics_cor} {input.radiomics_data} {params.outdir} {params.prefix}
        """

rule correlative_analysis_KEGG:
    input:
        genomics_filtered = OUTDIR + "{prefix}_KEGG_features_filtered.csv",
        radiomics_filtered = OUTDIR + "{prefix}_radiomics_features_filtered.csv"
    output:
        corr = OUTDIR + "{prefix}_KEGG_correlative_analysis.csv"
    params:
        group = "KEGG"
    resources:
        mem_mb=16000,
        runtime="2-00:00:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p $(dirname {output.corr})
        Rscript {RSCRIPTDIR}correlative_analysis.R {input.genomics_filtered} {input.radiomics_filtered} {wildcards.prefix} {output.corr}
        """

rule correlative_analysis_HALLMARK:
    input:
        genomics_filtered = OUTDIR + "{prefix}_HALLMARK_features_filtered.csv",
        radiomics_filtered = OUTDIR + "{prefix}_radiomics_features_filtered.csv"
    output:
        corr = OUTDIR + "{prefix}_HALLMARK_correlative_analysis.csv"
    params:
        group = "HALLMARK"
    resources:
        mem_mb=16000,
        runtime="2-00:00:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p $(dirname {output.corr})
        Rscript {RSCRIPTDIR}correlative_analysis.R {input.genomics_filtered} {input.radiomics_filtered} {wildcards.prefix} {output.corr}
        """

rule correlative_analysis_REACTOME:
    input:
        genomics_filtered = OUTDIR + "{prefix}_REACTOME_features_filtered.csv",
        radiomics_filtered = OUTDIR + "{prefix}_radiomics_features_filtered.csv"
    output:
        corr = OUTDIR + "{prefix}_REACTOME_correlative_analysis.csv"
    params:
        group = "REACTOME"
    resources:
        mem_mb=16000,
        runtime="2-00:00:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p $(dirname {output.corr})
        Rscript {RSCRIPTDIR}correlative_analysis.R {input.genomics_filtered} {input.radiomics_filtered} {wildcards.prefix} {output.corr}
        """

rule correlative_analysis_BIOCARTA:
    input:
        genomics_filtered = OUTDIR + "{prefix}_BIOCARTA_features_filtered.csv",
        radiomics_filtered = OUTDIR + "{prefix}_radiomics_features_filtered.csv"
    output:
        corr = OUTDIR + "{prefix}_BIOCARTA_correlative_analysis.csv"
    params:
        group = "BIOCARTA"
    resources:
        mem_mb=16000,
        runtime="2-00:00:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p $(dirname {output.corr})
        Rscript {RSCRIPTDIR}correlative_analysis.R {input.genomics_filtered} {input.radiomics_filtered} {wildcards.prefix} {output.corr}
        """
rule clinical_correlation_analysis:
    input:
        kegg = OUTDIR + "{prefix}_KEGG_features_filtered.csv",
        hallmark = OUTDIR + "{prefix}_HALLMARK_features_filtered.csv",
        reactome = OUTDIR + "{prefix}_REACTOME_features_filtered.csv",
        biocarta = OUTDIR + "{prefix}_BIOCARTA_features_filtered.csv",
        clinical = OUTDIR + "{prefix}_harmonized_clinical.csv"
    output:
        kegg_full = OUTDIR + "{prefix}_KEGG_clinical_correlations.csv",
        kegg_filtered = OUTDIR + "{prefix}_KEGG_clinical_correlations_filtered.csv",
        hallmark_full = OUTDIR + "{prefix}_HALLMARK_clinical_correlations.csv",
        hallmark_filtered = OUTDIR + "{prefix}_HALLMARK_clinical_correlations_filtered.csv",
        reactome_full = OUTDIR + "{prefix}_REACTOME_clinical_correlations.csv",
        reactome_filtered = OUTDIR + "{prefix}_REACTOME_clinical_correlations_filtered.csv",
        biocarta_full = OUTDIR + "{prefix}_BIOCARTA_clinical_correlations.csv",
        biocarta_filtered = OUTDIR + "{prefix}_BIOCARTA_clinical_correlations_filtered.csv"
    params:
        prefix = "{prefix}",
        outdir = lambda wildcards: OUTDIR.format(prefix=wildcards.prefix),
        output_prefix = lambda wildcards: OUTDIR.format(prefix=wildcards.prefix) + wildcards.prefix
    resources:
        mem_mb=8000,
        runtime="3-00:00:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p {params.outdir}
        Rscript {RSCRIPTDIR}clinical_correlations.R {input.clinical} {params.output_prefix} {input.kegg} {input.hallmark} {input.reactome} {input.biocarta}
        """