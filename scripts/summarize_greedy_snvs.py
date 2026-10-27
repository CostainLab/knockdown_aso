#!/usr/bin/env python3
"""
Concatenate all *.greedy_selected_snvs.tsv files into a
single master table, adding a 'gene' and 'transcript_id' column derived from each
filename (<GENE>__<TRANSCRIPT_ID>.greedy_selected_snvs.tsv).

"""

import os
import glob
import pandas as pd
import pysam
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description="Combine all greedy snv tsv files")
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--dbsnp", required=True)

args = parser.parse_args()

INPUT_DIR = Path(args.input)
OUTPUT_FILE = Path(args.output)
DBSNP = Path(args.dbsnp)

def parse_filename(path: str) -> tuple[str, str]:
    """Return (gene, transcript_id) from <GENE>__<TRANSCRIPT>.greedy_selected_snvs.tsv."""
    basename = os.path.basename(path)
    stem = basename.replace(".greedy_selected_snvs.tsv", "")
    if "__" in stem:
        gene, transcript_id = stem.split("__", 1)
    else:
        gene, transcript_id = stem, ""
    return gene, transcript_id

# MODIFICATION: Added function to pull rsID from NCBI file
# map chromosome numbers to ncbi format
chrom_map = {
    '1': 'NC_000001.11', '2': 'NC_000002.12', '3': 'NC_000003.12',
    '4': 'NC_000004.12', '5': 'NC_000005.10', '6': 'NC_000006.12',
    '7': 'NC_000007.14', '8': 'NC_000008.11', '9': 'NC_000009.12',
    '10': 'NC_000010.11', '11': 'NC_000011.10', '12': 'NC_000012.12',
    '13': 'NC_000013.11', '14': 'NC_000014.9', '15': 'NC_000015.10',
    '16': 'NC_000016.10', '17': 'NC_000017.11', '18': 'NC_000018.10',
    '19': 'NC_000019.10', '20': 'NC_000020.11', '21': 'NC_000021.9',
    '22': 'NC_000022.11', 'X': 'NC_000023.11', 'Y': 'NC_000024.10',
    'MT': 'NC_012920.1', 'M': 'NC_012920.1' # Included both MT and M just in case
}

# open dbSNP file
dbsnp = pysam.VariantFile(DBSNP)

# define function to pull IDs
def get_rsid(chrom, pos, ref, alt):
    # convert pos to int
    pos = int(pos)

    # remove 'chr' prefix from tsv
    clean_chrom = str(chrom).replace('chr', '')
    
    # translate chrom to NCBI terminology
    ncbi_chrom = chrom_map.get(clean_chrom)
    
    # if the chromosome isn't in dictionary, fail
    if not ncbi_chrom:
        return "Invalid_Chr"
        
    try:
        # use NCBI chrom to query dbsnp file
        for rec in dbsnp.fetch(ncbi_chrom, pos - 1, pos):
            if rec.ref == ref and alt in rec.alts:
                return rec.id if rec.id else "NA"
    except ValueError:
        pass 
    
    return "Not_Found"


def main() -> None:
    pattern = os.path.join(INPUT_DIR, "*.greedy_selected_snvs.tsv")
    files = sorted(glob.glob(pattern))

    if not files:
        print(f"No *.greedy_selected_snvs.tsv files found in: {INPUT_DIR}")
        return

    print(f"Found {len(files)} file(s) — reading...", flush=True)

    chunks = []
    for path in files:
        gene, transcript_id = parse_filename(path)
        df = pd.read_csv(path, sep="\t", dtype=str)

        # get rsID
        rsids = []
        for row in df.itertuples(index=False):
            rsids.append(str(get_rsid(row.chrom, row.pos, row.ref, row.alt)))

        df['rsID'] = rsids

        # pull gene and transcript, then append to chunks
        df.insert(0, "gene", gene)
        df.insert(1, "transcript_id", transcript_id)
        chunks.append(df)

    master = pd.concat(chunks, ignore_index=True)

    master.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"Master table: {len(master)} rows across {len(files)} genes")
    print(f"Written to:   {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
