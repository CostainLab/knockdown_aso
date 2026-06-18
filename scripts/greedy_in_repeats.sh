#!/bin/bash
#SBATCH --job-name=greedy_in_repeats
#SBATCH --output=../logs/greedy_repeats_%j.out
#SBATCH --error=../logs/greedy_repeats_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2

module load bedtools

INDIR="/hpf/projects/gregorycostain/kd_aso/greedy_out/g2p_logofunc_analysis"
OUTDIR="/hpf/projects/gregorycostain/kd_aso/repeats"
REPEATS="/hpf/projects/gregorycostain/kd_aso/repeats/master_repeats.bed.gz"
SUMMARY="${OUTDIR}/greedy_snps_in_repeats.tsv"

# make header for summary file
echo -e "gene\tchr\tstart\tend\tsnp_id" > "$SUMMARY"

for gene_file in "$INDIR"/*.greedy_selected_snps.tsv; do

    #get basename and gene name
    base=$(basename "$gene_file" .tsv)
    gene=${base%%__*}

    # define bed name
    bed="${OUTDIR}/${base}_in_repeats.bed"

    # convert file to bed
    awk 'BEGIN{OFS="\t"} {print $2, $3-1, $3, $4}' "$gene_file" > "$bed"

    # find overlap with repeats
    bedtools intersect \
        -a "$bed" \
        -b "$REPEATS" \
        -wa | \
    awk -v g="$gene" 'BEGIN{OFS="\t"} {print g, $1, $2, $3, $4}' >> "$SUMMARY"

done