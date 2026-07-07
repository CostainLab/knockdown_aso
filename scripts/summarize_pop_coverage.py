#!/usr/bin/env python3
import pandas as pd
import glob
import sys

def main():
    # grab every tally file
    file_list = glob.glob("coverage/*.sample_coverage_tally.tsv")

    if not file_list:
        print("No files found")
        sys.exit(1)

    print(f"Found {len(file_list)} coverage summaries. Compiling...")

    # read all files into list and stick together vertically
    df_list = [pd.read_csv(f, sep='\t') for f in file_list]
    master_df = pd.concat(df_list, ignore_index=True)

    # pivot dataframe from long to wide
    wide_df = master_df.pivot(
        index='gene',
        columns='super_pop',
        values=['covered', 'not_covered', 'total_samples', 'proportion_covered']
    )

    # flatten stacked column headers
    new_cols = []
    for metric, pop in wide_df.columns:
        new_cols.append(f"{pop}_{metric}")

    wide_df.columns = new_cols

    # bring gene back from index to column
    wide_df = wide_df.reset_index()

    # fill missing values with 0
    wide_df = wide_df.fillna(0)

    # write master table
    out_file = "coverage/master_gene_coverage_by_pop.tsv"
    wide_df.to_csv(out_file, sep='\t', index=False)
    master_df.to_csv("coverage/master_gene_coverage_by_pop_long.tsv", sep='\t', index=False)
    print(f"Success! Master file written to {out_file}")

if __name__ == "__main__":
    main()