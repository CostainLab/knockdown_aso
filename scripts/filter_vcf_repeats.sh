#!/bin/bash
#SBATCH --job-name=filter_common_snps
#SBATCH --output=../logs/repeat_filter_%j.out
#SBATCH --error=../logs/repeat_filter_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2

module load bcftools

INDIR='/hpf/projects/gregorycostain/kd_aso/1000G_snps/g2p_logofunc_analysis'
OUTDIR='/hpf/projects/gregorycostain/kd_aso/1000G_snps_filtered'
REPEATS='/hpf/projects/gregorycostain/kd_aso/repeats/master_repeats.bed.gz'
SUMMARY="${OUTDIR}/filtered_snp_summary.tsv"

mkdir -p "$OUTDIR"

# initialize table
echo -e "gene\tbefore\tafter\tremoved\tprop_removed" > "$SUMMARY"

# can use sbatch --export=TEST_VCF=/path/to/your/test.vcf.gz script.sh to test script on a single vcf
if [[ -n "$TEST_VCF" ]]; then
    vcfs="$TEST_VCF"
else
    vcfs="$INDIR"/*.vcf.gz
fi

# loop through vcf files
for vcf in $vcfs; do

    # extract base name
    base=$(basename "$vcf" .vcf.gz)

    # generate output names
    filtered_vcf="${OUTDIR}/${base}_filtered.vcf.gz"

    echo "Processing $base..."

    # count variants before filtering
    before_count=$(bcftools index -n "$vcf")

    # filter out variants in repeats
    bcftools view \
        -T ^"$REPEATS" \
        "$vcf" \
        -Oz -o "$filtered_vcf"

    # count variants before and after filtering and index
    after_count=0
    if [[ -s "$filtered_vcf" ]]; then
        # index
        bcftools index -t "$filtered_vcf"
        # count number of variants left
        after_count=$(bcftools index -n "$filtered_vcf")
    fi

    # add filtering summary to table
    removed=$((before_count - after_count))

    if [[ "$before_count" -gt 0 ]]; then
        prop_removed=$(awk -v r="$removed" -v b="$before_count" 'BEGIN {printf "%.1f", r/b}')
    else
        prop_removed="NA"
    fi

    echo -e "${base}\t${before_count}\t${after_count}\t${removed}\t${prop_removed}" >> "$SUMMARY"


done

echo "All genes filtered"