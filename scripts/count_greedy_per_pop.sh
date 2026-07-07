#!/bin/bash
#SBATCH --job-name=count_greedy
#SBATCH --output=../logs/count_greedy_%j.out
#SBATCH --error=../logs/count_greedy_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2

# counts hom and het individuals per super pop using file produced for LD analysis

shopt -s nullglob

module purge
module load bcftools

# define directories
INDIR='/hpf/projects/gregorycostain/kd_aso/1000G_LD'
OUTDIR='/hpf/projects/gregorycostain/kd_aso/pop_counts'

for pop_dir in "$INDIR"/*/; do

    # extract pop
    pop=$(basename "$pop_dir")
    echo "Processing $pop"

    # initialize output
    OUTPUT="${OUTDIR}/pop_counts_greedy_${pop}.tsv"
    echo -e "chr\tpos\tsnp_id\tsuper_pop\tnum_het\tnum_hom_alt\tpop_total" > "$OUTPUT"

    for vcf in "$pop_dir"/*_greedy_snps.vcf.gz; do

        bcftools query -f '%CHROM\t%POS\t%ID[\t%GT]\n' "$vcf" |
        awk  -v pop="$pop" 'BEGIN{OFS="\t"}
        {
            het=0
            homalt=0
            total=NF-3
        
            for(i=4;i<=NF;i++) {
                if($i~/^(0[\/|]1|1[\/|]0)$/) het++
                else if($i~/^1[\/|]1$/) homalt++
            }

            print $1,$2,$3,pop,het,homalt,total
        }' >> "$OUTPUT"
    done
done