#!/bin/bash
#SBATCH --job-name=clinvar_filtering
#SBATCH --output=../logs/clinvar_filtering_%j.out
#SBATCH --error=../logs/clinvar_filtering_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

module purge

python filter_clinvar_variants.py