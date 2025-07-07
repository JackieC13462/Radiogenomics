# Define config variables for input/output locations
configfile: "correlative_config.yaml"

COHORTS = [d["name"] for d in config["datasets"]]
PATHWAY_GROUPS = ["KEGG", "HALLMARK", "REACTOME", "BIOCARTA"]

def get_genomics_file(wildcards):
    return next(d["genomics_file"] for d in config["datasets"] if d["name"] == wildcards.prefix)

def get_radiomics_file(wildcards):
    return next(d["radiomics_file"] for d in config["datasets"] if d["name"] == wildcards.prefix)

OUTDIR = "/cluster/projects/bhklab/procdata/Radiogenomics/outputs/correlations/{prefix}/"
RSCRIPTDIR = "/cluster/home/t138199uhn/Scripts/"

rule all:
    input:
        expand(OUTDIR + "{prefix}_{group}_correlative_analysis.csv", prefix=COHORTS, group=PATHWAY_GROUPS)

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

rule radiomics_self_correlation:
    input:
        radiomics = OUTDIR + "{prefix}_radiogenomic_features.csv"
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

rule genomics_self_correlation:
    input:
        genomics = OUTDIR + "{prefix}_radiogenomic_RNAseq.csv"
    output:
        kegg = OUTDIR + "{prefix}_KEGG_self_correlation.csv",
        hallmark = OUTDIR + "{prefix}_HALLMARK_self_correlation.csv",
        reactome = OUTDIR + "{prefix}_REACTOME_self_correlation.csv",
        biocarta = OUTDIR + "{prefix}_BIOCARTA_self_correlation.csv"
    params:
        outdir = OUTDIR
    resources:
        mem_mb=8000,
        runtime="2:00:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p $(dirname {output.kegg})
        Rscript {RSCRIPTDIR}genomics_self_correlation.R {input.genomics} {wildcards.prefix} {params.outdir}
        """

rule feature_filtering:
    input:
        kegg = OUTDIR + "{prefix}_KEGG_self_correlation.csv",
        hallmark = OUTDIR + "{prefix}_HALLMARK_self_correlation.csv",
        reactome = OUTDIR + "{prefix}_REACTOME_self_correlation.csv",
        biocarta = OUTDIR + "{prefix}_BIOCARTA_self_correlation.csv",
        radiomics = OUTDIR + "{prefix}_radiomics_self_correlation.csv",
        genomics_data = OUTDIR + "{prefix}_radiogenomic_RNAseq.csv",
        radiomics_data = OUTDIR + "{prefix}_radiogenomic_features.csv"
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
            {input.kegg} {input.hallmark} {input.reactome} {input.biocarta} {input.radiomics} \
            {input.genomics_data} {input.radiomics_data} {params.outdir} {params.prefix}
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
        runtime="4-00:00:00",
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
        runtime="4-00:00:00",
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
        runtime="4-00:00:00",
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
        runtime="4-00:00:00",
        cpus=1
    shell:
        """
        module load R
        mkdir -p $(dirname {output.corr})
        Rscript {RSCRIPTDIR}correlative_analysis.R {input.genomics_filtered} {input.radiomics_filtered} {wildcards.prefix} {output.corr}
        """