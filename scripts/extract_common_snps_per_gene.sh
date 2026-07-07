#!/bin/bash
#SBATCH --job-name=gof_common_snps
#SBATCH --output=logs/gof_common_snps_%j.out
#SBATCH --error=logs/gof_common_snps_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

module purge
module load htslib
module load bcftools


BED="bed_files/gof_gene_regions_g2p_logo.bed"
VCF_DIR="1000G_highcov_unrelated"
OUT_DIR="1000G_snps"

mkdir -p "$OUT_DIR" logs

SUMMARY="${OUT_DIR}/common_snp_counts.tsv"
echo -e "gene\ttranscript\tchrom\tstart\tend\toutput_vcf\tn_snps" > "$SUMMARY"

# Generate tsv with location of each SNP for later ------------------------
TSV_OUT="${OUT_DIR}/common_snp_coordinates.tsv"
echo -e "gene\tchr\tstart\tsnp_id" > "$TSV_OUT"

while read -r chrom start end name score strand; do
    # Skip empty lines
    [[ -z "${chrom:-}" ]] && continue

    chrom_no="${chrom#chr}"

    if [[ "$chrom_no" == "X" ]]; then
        vcf="${VCF_DIR}/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2_unrelated.vcf.gz"
    else
        vcf="${VCF_DIR}/1kGP_high_coverage_Illumina.chr${chrom_no}.filtered.SNV_INDEL_SV_phased_panel_unrelated.vcf.gz"
    fi

    safe_name=$(echo "$name" | sed 's/|/__/g' | sed 's/[^A-Za-z0-9_.-]/_/g')
    out_vcf="${OUT_DIR}/${safe_name}.common_snps.vcf.gz"

    tmp_bed=$(mktemp)
    echo -e "${chrom}\t${start}\t${end}\t${name}\t${score}\t${strand}" > "$tmp_bed"

    bcftools view \
        -R "$tmp_bed" \
        -v snps \
        -i 'AF>0.01' \
        "$vcf" \
        -Oz -o "$out_vcf"

    tabix -p vcf "$out_vcf"

    n_snps=$(bcftools view -H "$out_vcf" | wc -l)

    gene=$(echo "$name" | cut -d'|' -f1)
    transcript=$(echo "$name" | cut -d'|' -f2)

    # append location of each common snp to generated file ------------------------
    bcftools query -f '%CHROM\t%POS\t%ID\n' "$out_vcf" |
    awk -v gene="$gene" '{print gene "\t" $1 "\t" $2 "\t" $3}' >> "$TSV_OUT"
    # ------------------------------------------------------------------------

    echo -e "${gene}\t${transcript}\t${chrom}\t${start}\t${end}\t${out_vcf}\t${n_snps}" >> "$SUMMARY"

    rm "$tmp_bed"

done < "$BED"
