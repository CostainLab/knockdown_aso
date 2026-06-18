#!/bin/bash
#SBATCH --job-name=1000G_filter_vcfs
#SBATCH --output=../logs/1000G_filter_%j.out
#SBATCH --error=../logs/1000G_filter_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

module load bcftools

# define directories
INDIR='/hpf/projects/gregorycostain/kd_aso/1000G_highcov'
OUTDIR='/hpf/projects/gregorycostain/kd_aso/1000G_highcov_unrelated'
SAMPLES='/hpf/projects/gregorycostain/kd_aso/1000G_unrelated_ids.txt'

# loop through files
for vcf in "$INDIR"/*vcf.gz; do
    # extract file base name
    base=$(basename "$vcf" .vcf.gz)
    # generate output name
    out_file="${OUTDIR}/${base}_unrelated.vcf.gz"

    echo "Processing $base -> $out_file"

    # use bcf tools to filter samples
    bcftools view -S "$SAMPLES" -O z "$vcf" > "$out_file"

    # index vcf
    bcftools index -t "$out_file"

done

echo "Filtering complete"