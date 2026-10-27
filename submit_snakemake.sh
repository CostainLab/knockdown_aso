#!/bin/bash
#SBATCH --job-name=snakemake_master
#SBATCH --output=logs/master_%j.log
#SBATCH --error=logs/master_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --partition=general

# Load profile and activate the master environment
source ~/.bashrc
conda activate kd_snakemake_env

# Unload conflicting system modules
module unload python 2>/dev/null

# Submit rules as independent Slurm jobs
snakemake \
    --executor slurm \
    --jobs 10 \
    --use-conda \
    --use-envmodules \
    --latency-wait 60 \
    --rerun-incomplete \
    --printshellcmds