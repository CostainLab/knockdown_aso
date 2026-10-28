#!/bin/bash

# Initialize variables
REPEATS=$1
INDIR=$2
OUTDIR=$3

mkdir -p "$OUTDIR"

VCF_OUTDIR="${OUTDIR}/filtered_vcfs"
mkdir -p "$VCF_OUTDIR"

# Initialize table
SUMMARY="${OUTDIR}/filtered_snv_summary.tsv"
echo -e "gene\tbefore\tafter\tremoved\tprop_removed" > "$SUMMARY"

# Loop through vcf files
shopt -s nullglob

vcfs=("$INDIR"/*.vcf.gz)

if [[ ${#vcfs[@]} -eq 0 ]]; then
    echo "No VCF files found in $INDIR"
    exit 1
fi

for vcf in "${vcfs[@]}"; do
    # Extract base name
    base=$(basename "$vcf" .common_snvs.vcf.gz)

    # Generate output names
    filtered_vcf="${VCF_OUTDIR}/${base}.common_snvs_filtered.vcf.gz"

    echo "Processing $base..."

    # Count variants before filtering
    before_count=$(bcftools index -n "$vcf")

    # Filter out variants in repeats
    bcftools view \
        -T ^"$REPEATS" \
        "$vcf" \
        -Oz -o "$filtered_vcf"

    # Count variants after filtering and index
    after_count=0
    if [[ -s "$filtered_vcf" ]]; then
        # index
        bcftools index -t "$filtered_vcf"
        # count number of variants left
        after_count=$(bcftools index -n "$filtered_vcf")
    fi

    # Add filtering summary to table
    removed=$((before_count - after_count))

    if [[ "$before_count" -gt 0 ]]; then
        prop_removed=$(awk -v r="$removed" -v b="$before_count" 'BEGIN {printf "%.1f", r/b}')
    else
        prop_removed="NA"
    fi

    echo -e "${base}\t${before_count}\t${after_count}\t${removed}\t${prop_removed}" >> "$SUMMARY"


done

echo "All genes filtered"