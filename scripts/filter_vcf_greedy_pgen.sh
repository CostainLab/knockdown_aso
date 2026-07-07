#!/bin/bash
#SBATCH --job-name=1000G_filter_vcfs_snps
#SBATCH --output=../logs/1000G_greedy_snps_%j.out
#SBATCH --error=../logs/1000G_greedy_snps_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --constraint=avx512

module purge
module load bcftools
module load plink2/a6LM

set -euo pipefail
shopt -s nullglob

# define directories
INDIR='/hpf/projects/gregorycostain/kd_aso/1000G_highcov_unrelated'
OUTDIR='/hpf/projects/gregorycostain/kd_aso/1000G_LD'
SNPS='/hpf/projects/gregorycostain/kd_aso/greedy_out/greedy_snp_ids.txt'
POP_DIR='/hpf/projects/gregorycostain/kd_aso/metadata/population_lists'
SEX='/hpf/projects/gregorycostain/kd_aso/metadata/1000G_sex.txt'

echo "VCFs found:"
printf '%s\n' "$INDIR"/*unrelated.vcf.gz

#--------- loop through populations -------#
for pop_file in "$POP_DIR"/*.txt; do

    #extract pop
    pop_name=$(basename "$pop_file" .txt)

    echo "Initiating pipeline for population: $pop_name"

    # output folder for pop
    POP_OUTDIR="${OUTDIR}/${pop_name}"
    mkdir -p "$POP_OUTDIR"

#--------- inner loop loops through vcfs -------#    
# loop through vcfs
    for vcf in "$INDIR"/*unrelated.vcf.gz; do
        # extract base name
        base=$(basename "$vcf" .vcf.gz)
        # generate output name
        filtered_vcf="${POP_OUTDIR}/${base}_${pop_name}_greedy_snps.vcf.gz"
        pgen_prefix="${POP_OUTDIR}/${base}_${pop_name}"

        echo "Processing $base..."

       
        # filter
        bcftools view \
            -i "ID=@$SNPS" \
            -S "$pop_file" \
            "$vcf" \
            --force-samples \
            -Oz -o "$filtered_vcf"

        if [[ -s "$filtered_vcf" ]] && [[ $(bcftools view -H "$filtered_vcf" | wc -l) -gt 0 ]]; then
            bcftools index -t "$filtered_vcf"

            # convert to pgen
            plink2 --vcf "$filtered_vcf" \
                   --make-pgen \
                   --split-par hg38 \
                   --update-sex "$SEX" \
                   --out "$pgen_prefix"
        else
            echo "Skipping empty VCF"
        fi
    done
done
