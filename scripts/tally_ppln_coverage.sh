#!/bin/bash
#SBATCH --job-name=tally_cov
#SBATCH --output=logs/tally_cov_%j.out
#SBATCH --error=logs/tally_cov_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

for cov in greedy_out/*.sample_coverage.tsv; do
    base=$(basename "$cov" .sample_coverage.tsv)
    gene="${base%%__*}"

    echo "Tallying coverage: $gene"

    python scripts/tally_ppln_coverage.py \
        --coverage "$cov" \
        --popdata "metadata/20130606_g1k_3202_samples_ped_population.txt" \
        --gene $gene \
        --out "coverage/${base}.sample_coverage_tally.tsv"
done