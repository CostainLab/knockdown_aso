#!/bin/bash
#SBATCH --job-name=proximal_snps
#SBATCH --output=../logs/proximal_snps_%j.out
#SBATCH --error=../logs/proximal_snps_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2

set -euo pipefail

# --- input files ---
GREEDY_TSV='/hpf/projects/gregorycostain/kd_aso/greedy_out/master_greedy_selected_snps_repeat_filtered_30Jun26.tsv'
COMMON_TSV='/hpf/projects/gregorycostain/kd_aso/1000G_snps/common_snp_coordinates.tsv'

# --- outputs ---
GREEDY_BED="/hpf/projects/gregorycostain/kd_aso/proximal_snps/greedy_snps.bed"
COMMON_BED="/hpf/projects/gregorycostain/kd_aso/proximal_snps/common_snps.bed"
OUTPUT="/hpf/projects/gregorycostain/kd_aso/proximal_snps/greedy_with_proximal_snps.tsv"

# --- load tools ---
module load bedtools

echo "Converting TSV to BED..."

# convert greedy TSV → BED
# needs to convert to 0-based w start/end
awk 'BEGIN{OFS="\t"} NR>1 {
    gene=$1
    chr=$4
    pos=$5
    snp_id=$6
    print chr, pos-1, pos, snp_id, gene
}' "$GREEDY_TSV" > "$GREEDY_BED"

# Convert common TSV → BED
awk 'BEGIN{OFS="\t"} NR>1 {
    chr=$2
    pos=$3
    snp_id=$4
    print chr, pos-1, pos, snp_id
}' "$COMMON_TSV" | sort -u > "$COMMON_BED"

echo "Running bedtools window (±19 bp)..."

# count common SNPs within ±19 bp of each greedy SNP, exclude the greedy SNP itself
bedtools window -c -a "$GREEDY_BED" -b "$COMMON_BED" -w 19 |
awk 'BEGIN{OFS="\t"} {
    corrected_count = $6 - 1

    if (corrected_count < 0) corrected_count = 0

    print $1, $2, $3, $4, $5, corrected_count
}' > /hpf/projects/gregorycostain/kd_aso/proximal_snps/counts.bed

echo "Formatting final output..."

# add header + convert back to readable format
echo -e "chr\tstart\tend\tsnp_id\tgene\tn_common_within_20bp" > "$OUTPUT"
sort -k1,1V -k2,2n -k5,5 /hpf/projects/gregorycostain/kd_aso/proximal_snps/counts.bed >> "$OUTPUT"

echo "Done!"
echo "Output: $OUTPUT"