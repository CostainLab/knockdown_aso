#!/bin/bash
#SBATCH --job-name=greedy_in_repeats
#SBATCH --output=../logs/greedy_repeats_%j.out
#SBATCH --error=../logs/greedy_repeats_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
set -euo pipefail

module load bedtools

GREEDY="/hpf/projects/gregorycostain/kd_aso/greedy_out/g2p_logofunc_analysis/master_greedy_selected_snps.tsv"
OUTDIR="/hpf/projects/gregorycostain/kd_aso/repeats"
REPEATS=(/hpf/projects/gregorycostain/kd_aso/repeats/*.bed.gz)
SUMMARY="${OUTDIR}/greedy_snps_in_repeats.tsv"

# make header for summary file
echo -e "gene\tchrom\tpos\tvariant_id" > "$SUMMARY"

# find overlap with repeats
awk 'NR>1{OFS="\t"; print $4, $5-1, $5, $1, $9}' "$GREEDY" | \
bedtools intersect \
    -a - \
    -b ${REPEATS} \
    -wa -u | \
awk 'BEGIN{OFS="\t"} {print $4, $1, $3, $5}' >> "$SUMMARY"

echo "Finished finding overlaps using: " $REPEATS