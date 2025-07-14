#!/bin/bash
#SBATCH --output=radiogenomics_pipeline_%j.out
#SBATCH --error=radiogenomics_pipeline_%j.err
#SBATCH --time=4-12:00:00
#SBATCH --mem=8000
#SBATCH --cpus-per-task=1

# Load Snakemake 
module load snakemake
module load R

# Move to the directory containing your Snakefile
cd /cluster/home/t138199uhn/Scripts 

# Make sure this directory exists!
mkdir -p /cluster/home/t138199uhn/slurm/slurm_logs

# Run Snakemake with SLURM job submission for each rule
snakemake --snakefile correlative_snakefile.snakefile \
          --configfile correlative_config.yaml \
          --jobs 50 \
          --latency-wait 120 \
          --cluster "sbatch --job-name=smk_{rule} --mem={resources.mem_mb} --cpus-per-task={resources.cpus} --time={resources.runtime} --output=/cluster/home/t138199uhn/slurm/slurm_logs/{rule}.{wildcards}.%j.out --error=/cluster/home/t138199uhn/slurm/slurm_logs/{rule}.{wildcards}.%j.err"

# Notes:
# --jobs 50 allows up to 50 jobs in parallel (adjust as needed for your cluster).
# --cluster submits each job as its own SLURM job, using the resources specified in your Snakefile.
# Log files will be named by rule and wildcards for easier tracking.
# Make sure the slurm_logs directory exists or is created automatically.
