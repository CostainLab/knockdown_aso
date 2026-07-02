#!/bin/bash
#SBATCH --job-name=proximal_clinvar
#SBATCH --output=../logs/proximal_clinvar_%j.out
#SBATCH --error=../logs/proximal_clinvar_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2

set -euo pipefail

# --- input files ---
GREEDY_TSV='/hpf/projects/gregorycostain/kd_aso/greedy_out/master_greedy_selected_snps_repeat_filtered_30Jun26.tsv'
CLINVAR_TSV='/hpf/projects/gregorycostain/kd_aso/clinvar/variant_summary/clinvar_gof_12Jun26.tsv'

# --- outputs ---
GREEDY_BED="/hpf/projects/gregorycostain/kd_aso/clinvar/greedy_snps.bed"
CLINVAR_BED="/hpf/projects/gregorycostain/kd_aso/clinvar/clinvar.bed"
OUTPUT="/hpf/projects/gregorycostain/kd_aso/clinvar/greedy_with_proximal_clinvar.tsv"

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

# Convert clinvar TSV → BED
# Add chr prefix
awk 'BEGIN{OFS="\t"} NR>1 {
    gene=$6
    chr="chr"$1
    pos=$2
    ref=$3
    alt=$4
    print chr, pos-1, pos, gene, ref, alt
}' "$CLINVAR_TSV" | sort -u > "$CLINVAR_BED"

# Remove duplicate clinvar variants if present
sort -k1,1 -k2,2n -k3,3n "$CLINVAR_BED" | uniq > "${CLINVAR_BED}.tmp" && mv "${CLINVAR_BED}.tmp" "$CLINVAR_BED"

echo "Running bedtools window (±19 bp + variant = 20 bp)..."

# count common SNPs within ±19 bp of each greedy SNP, exclude the greedy SNP itself
bedtools window -c -a "$GREEDY_BED" -b "$CLINVAR_BED" -w 19 \
> /hpf/projects/gregorycostain/kd_aso/clinvar/counts.bed

echo "Formatting final output..."

# add header + convert back to readable format
echo -e "chr\tstart\tend\tsnp_id\tgene\tn_common_within_20bp" > "$OUTPUT"
sort -k1,1V -k2,2n -k5,5 /hpf/projects/gregorycostain/kd_aso/clinvar/counts.bed >> "$OUTPUT"

echo "Done!"
echo "Output: $OUTPUT"