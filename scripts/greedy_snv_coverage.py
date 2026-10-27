#!/usr/bin/env python3

import argparse
import gzip
from typing import Dict, Set, List


def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def is_het(gt: str) -> bool:
    """
    Return True if genotype is heterozygous.

    Handles phased and unphased GTs:
    0|1, 1|0, 0/1, 1/0

    For now, we treat any two different non-missing alleles as heterozygous.
    """
    gt = gt.split(":")[0]

    if gt in {".", "./.", ".|."}:
        return False

    sep = "|" if "|" in gt else "/"

    alleles = gt.split(sep)

    if len(alleles) != 2:
        return False

    if "." in alleles:
        return False

    return alleles[0] != alleles[1]


def parse_vcf(vcf_path):
    samples: List[str] = []
    snp_to_hets: Dict[str, Set[str]] = {}
    snp_info: Dict[str, Dict[str, str]] = {}

    with open_vcf(vcf_path) as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.split("\t")
                samples = header[9:]
                continue

            fields = line.split("\t")

            chrom = fields[0]
            pos = fields[1]
            snp_id = fields[2]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            genotypes = fields[9:]

            # Make a stable SNP name
            variant_key = f"{chrom}:{pos}:{ref}:{alt}"

            het_samples = {
                sample
                for sample, gt in zip(samples, genotypes)
                if is_het(gt)
            }

            snp_to_hets[variant_key] = het_samples
            snp_info[variant_key] = {
                "chrom": chrom,
                "pos": pos,
                "id": snp_id,
                "ref": ref,
                "alt": alt,
                "info": info,
                "total_het": str(len(het_samples)),
            }

    return samples, snp_to_hets, snp_info


def run_greedy(samples, snp_to_hets, snp_info, coverage_threshold):
    total_n = len(samples)
    remaining = set(samples)
    selected = []
    rank = 1

    while len(remaining) / total_n > (1 - coverage_threshold):
        best_snp = None
        best_newly_covered = set()

        for snp, hets in snp_to_hets.items():
            newly_covered = hets & remaining

            if len(newly_covered) > len(best_newly_covered):
                best_snp = snp
                best_newly_covered = newly_covered

        if best_snp is None or len(best_newly_covered) == 0:
            break

        remaining -= best_newly_covered

        cumulative_covered = total_n - len(remaining)
        cumulative_fraction = cumulative_covered / total_n
        remaining_fraction = len(remaining) / total_n

        selected.append({
            "rank": rank,
            "variant_key": best_snp,
            **snp_info[best_snp],
            "newly_covered": len(best_newly_covered),
            "cumulative_covered": cumulative_covered,
            "cumulative_fraction": cumulative_fraction,
            "remaining_fraction": remaining_fraction,
        })

        # Remove selected SNP from future consideration
        del snp_to_hets[best_snp]

        rank += 1

    # MODIFICATION: Return the 'remaining' set alongside 'selected'
    return selected, remaining


def write_outputs(selected, remaining, samples, snp_to_hets_original, snp_info, out_prefix):
    selected_out = f"{out_prefix}.greedy_selected_snvs.tsv"
    all_counts_out = f"{out_prefix}.all_snv_het_counts.tsv"
    # MODIFICATION: Define new output file for sample coverage
    sample_coverage_out = f"{out_prefix}.sample_coverage.tsv"

    with open(selected_out, "w") as out:
        out.write(
            "rank\tchrom\tpos\tsnp_id\tref\talt\tvariant_key\t"
            "total_het\tnewly_covered\tcumulative_covered\t"
            "cumulative_fraction\tremaining_fraction\n"
        )

        for row in selected:
            out.write(
                f"{row['rank']}\t{row['chrom']}\t{row['pos']}\t{row['id']}\t"
                f"{row['ref']}\t{row['alt']}\t{row['variant_key']}\t"
                f"{row['total_het']}\t{row['newly_covered']}\t"
                f"{row['cumulative_covered']}\t"
                f"{row['cumulative_fraction']:.6f}\t"
                f"{row['remaining_fraction']:.6f}\n"
            )

    with open(all_counts_out, "w") as out:
        out.write("chrom\tpos\tsnp_id\tref\talt\tvariant_key\ttotal_het\n")

        for snp, hets in sorted(
            snp_to_hets_original.items(),
            key=lambda x: len(x[1]),
            reverse=True
        ):
            info = snp_info[snp]
            out.write(
                f"{info['chrom']}\t{info['pos']}\t{info['id']}\t"
                f"{info['ref']}\t{info['alt']}\t{snp}\t{len(hets)}\n"
            )

    # Write the sample coverage status to the new file
    with open(sample_coverage_out, "w") as out:
        out.write("sample_id\tcoverage_status\n")
        for sample in samples:
            status = "Not_Covered" if sample in remaining else "Covered"
            out.write(f"{sample}\t{status}\n")

    print(f"Wrote: {selected_out}")
    print(f"Wrote: {all_counts_out}")
    print(f"Wrote: {sample_coverage_out}")


def main():
    parser = argparse.ArgumentParser(
        description="Run greedy SNP heterozygote coverage algorithm on a per-gene VCF."
    )

    parser.add_argument(
        "--vcf",
        required=True,
        help="Input per-gene common SNP VCF.gz"
    )

    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix"
    )

    parser.add_argument(
        "--coverage-threshold",
        type=float,
        default=0.95,
        help="Stop when cumulative coverage reaches this fraction. Default: 0.95"
    )

    args = parser.parse_args()

    samples, snp_to_hets, snp_info = parse_vcf(args.vcf)

    if len(samples) == 0:
        raise ValueError(f"No samples found in {args.vcf}")

    # Keep original copy for QC output
    snp_to_hets_original = {k: set(v) for k, v in snp_to_hets.items()}

    print(f"Input VCF: {args.vcf}")
    print(f"Number of samples: {len(samples)}")
    print(f"Number of SNVs: {len(snp_to_hets)}")

    # MODIFICATION: Unpack both variables returned by run_greedy
    selected, remaining = run_greedy(
        samples=samples,
        snp_to_hets=snp_to_hets,
        snp_info=snp_info,
        coverage_threshold=args.coverage_threshold
    )

    if selected:
        final = selected[-1]
        print(f"Selected SNVs: {len(selected)}")
        print(f"Final cumulative coverage: {final['cumulative_fraction']:.4f}")
        print(f"Final remaining fraction: {final['remaining_fraction']:.4f}")
    else:
        print("No SNVs selected.")

    # Pass remaining and samples down to write_outputs
    write_outputs(
        selected=selected,
        remaining=remaining,
        samples=samples,
        snp_to_hets_original=snp_to_hets_original,
        snp_info=snp_info,
        out_prefix=args.out_prefix
    )


if __name__ == "__main__":
    main()