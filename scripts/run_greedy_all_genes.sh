#!/bin/bash
#SBATCH --job-name=greedy_all
#SBATCH --output=logs/greedy_all_%j.out
#SBATCH --error=logs/greedy_all_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

mkdir -p greedy_out logs

for vcf in 1000G_snps_filtered/*.common_snps_filtered.vcf.gz; do
    base=$(basename "$vcf" .common_snps.vcf.gz)

    echo "Running greedy algorithm for: $base"

    python scripts/greedy_snp_coverage_v2.py \
        --vcf "$vcf" \
        --out-prefix "greedy_out/${base}"
done
