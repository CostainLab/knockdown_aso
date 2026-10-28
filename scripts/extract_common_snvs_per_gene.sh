#!/bin/bash

set -euo pipefail

# Initialize variables
BED=$1
VCF_DIR=$2
OUT_DIR=$3

VCF_OUTDIR="${OUT_DIR}/vcfs"
mkdir -p "$VCF_OUTDIR"
mkdir -p "$OUT_DIR"

# Initialize summary file per gene
SUMMARY="${OUT_DIR}/common_snv_counts.tsv"
echo -e "gene\ttranscript\tchrom\tstart\tend\toutput_vcf\tn_snvs" > "$SUMMARY"

# Initialize file with coordinates for each SNV
TSV_OUT="${OUT_DIR}/common_snv_coordinates.tsv"
echo -e "gene\tchr\tstart\tref\talt\tsnv_id" > "$TSV_OUT"

while read -r chrom start end name score strand; do
    # Skip empty lines
    [[ -z "${chrom:-}" ]] && continue

    chrom_no="${chrom#chr}"

    # File for X chrom is named slightly differently
    if [[ "$chrom_no" == "X" ]]; then
        vcf="${VCF_DIR}/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    else
        vcf="${VCF_DIR}/1kGP_high_coverage_Illumina.chr${chrom_no}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    fi

    safe_name=$(echo "$name" | sed 's/|/__/g' | sed 's/[^A-Za-z0-9_.-]/_/g')
    out_vcf="${VCF_OUTDIR}/${safe_name}.common_snvs.vcf.gz"

    tmp_bed=$(mktemp)
    echo -e "${chrom}\t${start}\t${end}\t${name}\t${score}\t${strand}" > "$tmp_bed"

    bcftools view \
        -R "$tmp_bed" \
        -v snps \
        -i 'AF>0.01' \
        "$vcf" \
        -Oz -o "$out_vcf"

    tabix -p vcf "$out_vcf"

    n_snvs=$(bcftools view -H "$out_vcf" | wc -l)

    gene=$(echo "$name" | cut -d'|' -f1)
    transcript=$(echo "$name" | cut -d'|' -f2)

    # Append location of each common snv to generated file
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' "$out_vcf" |
    awk -v gene="$gene" '{print gene "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' >> "$TSV_OUT"

    # Generate tsv file with per-gene summary
    echo -e "${gene}\t${transcript}\t${chrom}\t${start}\t${end}\t${out_vcf}\t${n_snvs}" >> "$SUMMARY"

    rm "$tmp_bed"

done < "$BED"
