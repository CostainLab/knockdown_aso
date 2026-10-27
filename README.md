# SNV target selection for knockdown ASO pipeline

This repository contains the scripts used for the manuscript "Population-scale genome-wide prioritization of targets for allele-specific antisense oligonucleotide therapy design". The main scripts to run the greedy algorithm with 1kGP data have been integrated into a Snakemake workflow (Snakefile) that provides a reproducible run order that mirrors the manual steps used for the published analysis.

**Top-level overview (as implemented in the `Snakefile`):**

NOTE: Directories #3-6 will be created during the pipeline.

- `1-BED_File`: Generates BED file containing coordinates of genes of interest (MANE-Select transcript) with `generate_gene_bed_files.py`.
- `2-1kGP_VCF`: Directory containing 1kGP VCF files used for greedy algorithm.
- `3-Common_SNVs`: Extracts common SNVs (MAF >1%) from VCF files with `extract_common_snvs_per_gene.sh`.
- `4-Repeat_Filtered_SNVs`: Filters common SNVs to exclude those in repeat regions with `filter_vcf_repeats.sh`.
- `5-Greedy_Output`: Applies greedy algorithm to iteratively select SNVs that maximize population coverage with `greedy_snv_coverage.py`. Additionally, outputs per-superpopulation coverage per gene.
- `6-Overall_Results`: Contains summary files for all greedy-selected SNVs and superpopulation coverage with `summarize_greedy_snvs.py` and `tally_ppln_coverage.py`.

**Prerequisites:**

Snakemake

- Install Snakemake (recommended via conda or mamba):

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake mamba -y
conda activate snakemake
```

- Ensure `Rscript` and `python3` are available. To use per-rule conda environments, Snakemake will create the envs from the `envs/` YAML files.

Gene list & Gene annotation file

- To run this pipeline on a custom gene set, modify `1-BED_File/gene_list.csv` to contain your list of genes (one per row) under the column title "gene_symbol".
- `gencode.v49.annotation.gtf` must be downloaded and added to `1-BED_File`. This annotation is publicly available here: https://www.gencodegenes.org/human/release_49.html

1kGP VCF files

- Per-chromosome 1kGP VCF and index files (.vcf.gz and .vcf.gz.tbi) must be added to "2-1kGP_VCF" before running the pipeline.
- VCF files are publicly available and can be downloaded from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/.
- After downloading, you can filter to exclude the offspring of trios (to reduce relatedness in the dataset) with `metadata/1kGP_unrelated_ids.txt` and, for example, `bcftools` using 1kGP VCF files.

dbSNP file

- To pull rsIDs for SNVs of interest, the file containing rsIDs for each SNV must be downloaded form dbSNP and added to `metadata`
- Download GCF_000001405.40.gz and GCF_000001405.40.gz.tbi from https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/

**Running the full workflow:**

NOTE: This pipeline has been optimized for use on a compute cluster and will need to be modified for local use or different clusters.

From the repository root in an interactive session, run:

```bash
sbatch submit_snakemake.sh
```
