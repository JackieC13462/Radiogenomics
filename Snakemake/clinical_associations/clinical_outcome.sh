#!/bin/bash
#SBATCH --output=/cluster/home/t138199uhn/slurm/slurm_logs/clinical_outcomes/clinical_pipeline_%j.out
#SBATCH --error=/cluster/home/t138199uhn/slurm/slurm_logs/clinical_outcomes/clinical_pipeline_%j.err
#SBATCH --time=2-12:00:00
#SBATCH --mem=8000
#SBATCH --cpus-per-task=1

# Load Snakemake and R modules
module load snakemake
module load R

# Move to the directory containing your Snakefile
cd /cluster/home/t138199uhn/Scripts

# Make sure this directory exists for logs!
mkdir -p /cluster/home/t138199uhn/slurm/slurm_logs/clinical_outcomes

# Run Snakemake with SLURM job submission for each rule
snakemake --snakefile association_snakefile.snakefile \
          --configfile clinical_config.yaml \
          --jobs 100 \
          --latency-wait 60 \
          --rerun-incomplete \
          --keep-going \
          --use-conda \
          --cluster "sbatch --job-name=smk_{rule} --mem={resources.mem_mb} --cpus-per-task=1 --time={resources.time} --output=/cluster/home/t138199uhn/slurm/slurm_logs/clinical_outcomes/{rule}.{wildcards}.%j.out --error=/cluster/home/t138199uhn/slurm/slurm_logs/clinical_outcomes/{rule}.{wildcards}.%j.err"

# Notes:
# --jobs 100 allows up to 100 jobs in parallel (adjust as needed for your cluster).
# --cluster submits each job as its own SLURM job, using the resources specified in your Snakefile.
# Log files will be named by rule and wildcards for easier tracking.
# Make sure the slurm_logs directory exists or is created automatically.
# Adjust --cpus-per-task if your rules use more than 1 CPU.
