#!/usr/bin/env python3
"""
Generate a single regions.bed file for all unique genes in specified gene list csv file.
Genomic coordinates are pulled from the GENCODE GTF using the MANE Select transcript.
Each row covers the full transcript span (UTRs + introns + exons) for one gene.
Output: autosomal_dominant_regions.bed
"""

import os
import logging
import pandas as pd
from pathlib import Path
import argparse

# ---- Logging ----

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s"
)

# ---- User inputs ----

parser = argparse.ArgumentParser(description="Generate bed file for coordinates of genes using MANE Select transcript")
parser.add_argument("--gene_list", "-g", default="gene_list.csv", help="List of gene symbols")
parser.add_argument("--annotation", "-a", default="gencode.v49.annotation.gtf", help="GTF annotation file")
parser.add_argument("--output", "-o", default="gene_regions.bed", help="Output BED file")

args = parser.parse_args()

GENE_LIST_CSV = Path(args.gene_list)
GTF_FILE = Path(args.annotation)
OUTPUT_BED = Path(args.output)

# ── Step 1: collect unique gene symbols ──────────────────────────────────────
def load_gene_symbols(csv_path: str) -> list[str]:
    df = pd.read_csv(csv_path, usecols=["gene_symbol"], dtype=str)
    df["gene_symbol"] = df["gene_symbol"].str.strip()
    return df["gene_symbol"].dropna().unique().tolist()


# ── Step 2: parse GTF for MANE Select transcript coordinates ─────────────────
GTF_COLS = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]

def extract_attr(series: pd.Series, key: str) -> pd.Series:
    """Extract a single-value attribute (e.g. gene_name) from the attributes column."""
    pattern = rf'{key} "([^"]+)"'
    return series.str.extract(pattern, expand=False)


def load_mane_transcripts(gtf_path: str, genes_of_interest: set[str]) -> pd.DataFrame:
    """
    Read the GTF in chunks, keeping only MANE Select transcript rows for genes
    of interest. Returns a DataFrame with one row per gene (first MANE Select
    transcript encountered).
    """
    print(f"Scanning GTF: {gtf_path}", flush=True)

    chunks = []
    found_genes: set[str] = set()
    chunksize = 200_000

    reader = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=GTF_COLS,
        dtype={"start": "int32", "end": "int32"},
        chunksize=chunksize,
        low_memory=False,
    )

    for i, chunk in enumerate(reader):
        # Keep only transcript rows that mention MANE_Select
        chunk = chunk[
            (chunk["feature"] == "transcript") &
            chunk["attributes"].str.contains("MANE_Select", na=False)
        ]

        if chunk.empty:
            continue

        # Extract gene_name and transcript_id
        chunk = chunk.copy()
        chunk["gene_name"]     = extract_attr(chunk["attributes"], "gene_name")
        chunk["transcript_id"] = extract_attr(chunk["attributes"], "transcript_id")

        # Keep only genes we care about that haven't been found yet
        chunk = chunk[
            chunk["gene_name"].isin(genes_of_interest) &
            ~chunk["gene_name"].isin(found_genes)
        ]

        if not chunk.empty:
            # One row per gene (first occurrence)
            chunk = chunk.drop_duplicates(subset="gene_name", keep="first")
            chunks.append(chunk[["chrom", "start", "end", "strand", "gene_name", "transcript_id"]])
            found_genes.update(chunk["gene_name"].tolist())

        print(f"  … chunk {i+1} processed, {len(found_genes)}/{len(genes_of_interest)} genes found", flush=True)

        if found_genes >= genes_of_interest:
            print("  All genes resolved – stopping early.", flush=True)
            break

    if not chunks:
        return pd.DataFrame(columns=["chrom", "start", "end", "strand", "gene_name", "transcript_id"])

    result = pd.concat(chunks, ignore_index=True)
    # Convert GTF 1-based inclusive start to BED 0-based
    result["start"] = result["start"] - 1
    return result


# ── Step 3: write a single combined BED file ──────────────────────────────────
def write_bed_file(df: pd.DataFrame, bed_path: str) -> None:
    os.makedirs(os.path.dirname(bed_path), exist_ok=True)
    df = df.sort_values(["chrom", "start"]).copy()
    df["name"]  = df["gene_name"] + "|" + df["transcript_id"]
    df["score"] = 0
    bed = df[["chrom", "start", "end", "name", "score", "strand"]]
    bed.to_csv(bed_path, sep="\t", header=False, index=False)
    print(f"Wrote {len(df)} gene regions to: {bed_path}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main() -> None:
    # 1. Load gene symbols
    genes = load_gene_symbols(GENE_LIST_CSV)
    print(f"Unique genes in gene list: {len(genes)}")
    genes_set = set(genes)

    # 2. Parse GTF
    mane_df = load_mane_transcripts(GTF_FILE, genes_set)

    # 3. Report any genes without a MANE Select transcript
    missing = genes_set - set(mane_df["gene_name"].tolist())
    if missing:
        print(f"\nWARNING: No MANE Select transcript found for {len(missing)} gene(s):")
        for g in sorted(missing):
            print(f"  {g}")

    # 4. Write combined BED file
    print()
    write_bed_file(mane_df, OUTPUT_BED)


if __name__ == "__main__":
    main()
