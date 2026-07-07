#!/bin/bash
#SBATCH --job-name=common_snp_coords
#SBATCH --output=../logs/common_snp_coords_%j.out
#SBATCH --error=../logs/common_snp_coords_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2

# Integrated into extract_common_snps_per_gene

module load bcftools

INDIR='/hpf/projects/gregorycostain/kd_aso/1000G_snps'

# initialize TSV file
OUT='/hpf/projects/gregorycostain/kd_aso/1000G_snps/common_snp_coordinates.tsv'
echo -e "gene\tchr\tpos\tsnp_id" > "$OUT"

for gene_snps in "$INDIR"/*.common_snps.vcf.gz; do

    # extract gene name
    file=${gene_snps##*/}   # strip path
    gene=${file%%__*}

    echo "Compiling "$gene" common SNPs"

    # extract info
    bcftools query -f '%CHROM\t%POS\t%ID\n' "$gene_snps" |
    awk -v gene="$gene" '{print gene "\t" $1 "\t" $2 "\t" $3}' >> "$OUT"

done
