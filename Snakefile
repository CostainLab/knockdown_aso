configfile: "config.yaml"

P = config["paths"]
R = config["references"]
SCRIPTS = P["scripts_dir"]

rule all:
    input:
        f"{P['results_dir']}/greedy_selected_snvs.tsv",
        f"{P['results_dir']}/gene_coverage_by_pop.tsv",
        f"{P['results_dir']}/gene_coverage_by_pop_long.tsv"

rule generate_bed:
    input:
        gene_list=R["gene_list"],
        annotation=R["annotation_gtf"]
    output:
        f"{P['bed_dir']}/gene_regions.bed"
    log:
        f"{P['logs_dir']}/generate_bed_file.log"
    conda: "envs/py.yaml"
    threads: 1
    resources:
        mem_mb=4000,
        time_min=30
    shell:
        "python3 {SCRIPTS}/generate_gene_bed_files.py "
        "--gene_list '{input.gene_list}' "
        "--annotation '{input.annotation}' "
        "--output '{output}' > {log} 2>&1"

rule extract_snvs:
    input:
        bed=f"{P['bed_dir']}/gene_regions.bed"
    output:
        summary=f"{P['common_snv_dir']}/common_snv_counts.tsv",
        coords=f"{P['common_snv_dir']}/common_snv_coordinates.tsv",
        vcf_dir=directory(f"{P['common_snv_dir']}/vcfs")
    params:
        vcf_dir=P["vcf_dir"],
        out_dir=P["common_snv_dir"]
    log:
        f"{P['logs_dir']}/extract_common_snvs.log"
    envmodules:
        "bcftools",
        "htslib"
    threads: 2
    resources:
        mem_mb=16000,
        time_min=480
    shell:
        "bash {SCRIPTS}/extract_common_snvs_per_gene.sh "
        "'{input.bed}' '{params.vcf_dir}' '{params.out_dir}' > {log} 2>&1"

rule filter_snvs:
    input:
        vcf_dir=f"{P['common_snv_dir']}/vcfs",
        repeats=R["repeats_bed"]
    output:
        summary=f"{P['filtered_snv_dir']}/filtered_snv_summary.tsv",
        filtered_vcf_dir=directory(f"{P['filtered_snv_dir']}/filtered_vcfs")
    params:
        outdir=P["filtered_snv_dir"]
    log:
        f"{P['logs_dir']}/filter_common_snvs.log"
    envmodules:
        "bcftools"
    threads: 2
    resources:
        mem_mb=16000,
        time_min=480
    shell:
        "bash {SCRIPTS}/filter_vcf_repeats.sh "
        "'{input.repeats}' '{input.vcf_dir}' '{params.outdir}' > {log} 2>&1"

rule run_greedy:
    input:
        vcf_dir=f"{P['filtered_snv_dir']}/filtered_vcfs"
    output:
        greedy_dir=directory(P["greedy_dir"])
    log:
        f"{P['logs_dir']}/greedy_algorithm.log"
    conda: "envs/py.yaml"
    threads: 2
    resources:
        mem_mb=16000,
        time_min=480
    shell:
        """
        mkdir -p "{output.greedy_dir}"

        for vcf in "{input.vcf_dir}"/*.common_snvs_filtered.vcf.gz; do
            base=$(basename "$vcf" .common_snvs_filtered.vcf.gz)

            echo "Running greedy algorithm for: $base"

            python3 {SCRIPTS}/greedy_snv_coverage.py \
                --vcf "$vcf" \
                --out-prefix "{output.greedy_dir}/${{base}}"
        done > {log} 2>&1
        """

rule merge_greedy:
    input:
        greedy_dir=P["greedy_dir"],
        dbsnp=R["dbsnp"]
    output:
        f"{P['results_dir']}/greedy_selected_snvs.tsv"
    log:
        f"{P['logs_dir']}/summarize_greedy.log"
    conda: "envs/py.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        time_min=480
    shell:
        "python3 {SCRIPTS}/summarize_greedy_snvs.py "
        "--input '{input.greedy_dir}' "
        "--output '{output}' "
        "--dbsnp '{input.dbsnp}' > {log} 2>&1"

rule tally_coverage:
    input:
        greedy_dir=P["greedy_dir"],
        popdata=R["popdata"]
    output:
        wide=f"{P['results_dir']}/gene_coverage_by_pop.tsv",
        long=f"{P['results_dir']}/gene_coverage_by_pop_long.tsv"
    params:
        out_prefix=f"{P['results_dir']}/gene_coverage"
    log:
        f"{P['logs_dir']}/tally_ppl_coverage.log"
    conda: "envs/py.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        time_min=480
    shell:
        """
        python3 {SCRIPTS}/tally_ppln_coverage.py \
            --popdata '{input.popdata}' \
            --coverage_files "{input.greedy_dir}"/*.sample_coverage.tsv \
            --out_prefix '{params.out_prefix}' > {log} 2>&1
        """