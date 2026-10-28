#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path
import sys

def main():
    # Define arguments for script
    parser = argparse.ArgumentParser(
        description="Tally sample coverage proportions per population."
        )
    
    parser.add_argument(
        "--popdata",
        required = True,
        help="1000G population file"
    )

    parser.add_argument(
        "--coverage_files",
        nargs='+',
        required=True,
        help="List of input sample_coverage.tsv files (e.g., coverage/*.sample_coverage.tsv)"
    )
    parser.add_argument(
        "--out_prefix",
        default="coverage/gene_coverage",
        help="Prefix for output master files"
    )

    args = parser.parse_args()

    if not args.coverage_files:
        print("No coverage files provided.")
        sys.exit(1)

    print(f"Loading population data...")

    # Load population data file
    popdata_df = pd.read_csv(
        args.popdata, 
        sep=' ', 
        usecols=['SampleID', 'Population', 'Superpopulation']
    )

    popdata_df = popdata_df.rename(columns={
        'SampleID': 'sample_id', 
        'Population': 'pop',
        'Superpopulation': 'super_pop'
    })

    print(f"Processing {len(args.coverage_files)} coverage files...")
    df_list = []

# Process each coverage file
    for cov_file in args.coverage_files:
        # Extract gene name from filename (Assumes format like "BRCA1__ENST.sample_coverage_tally.tsv")
        gene_name = Path(cov_file).name.split('__')[0]

        # Load individual coverage file
        cov_df = pd.read_csv(cov_file, sep='\t')

        # Merge datasets
        cov_pop_df = pd.merge(cov_df, popdata_df, on='sample_id', how='left')
        cov_pop_df['super_pop'] = cov_pop_df['super_pop'].fillna('Unknown')

        # Group by super_pop and count by coverage status
        tally = cov_pop_df.groupby(['super_pop', 'coverage_status']).size().unstack(fill_value=0)

        # Ensure both columns exist
        for status in ['Covered', 'Not_Covered']:
            if status not in tally.columns:
                tally[status] = 0
        
        # Calculate totals and proportions
        tally['Total_Samples'] = tally['Covered'] + tally['Not_Covered']
        tally['Proportion_Covered'] = tally['Covered'] / tally['Total_Samples']
        tally['Gene'] = gene_name

        # Clean df structure for master compilation
        tally = tally.reset_index()
        final_cols = ['Gene', 'super_pop', 'Covered', 'Not_Covered', 'Total_Samples', 'Proportion_Covered']
        tally = tally[final_cols]
        tally.columns = tally.columns.str.lower()
        
        df_list.append(tally)

    print("Compiling master tables...")
    # Concatenate all processed dataframes
    master_df = pd.concat(df_list, ignore_index=True)

# Pivot dataframe from long to wide
    wide_df = master_df.pivot(
        index='gene',
        columns='super_pop',
        values=['covered', 'not_covered', 'total_samples', 'proportion_covered']
    )

    # Flatten stacked column headers
    new_cols = [f"{pop}_{metric}" for metric, pop in wide_df.columns]
    wide_df.columns = new_cols

    # Bring gene back from index to column, fill NAs
    wide_df = wide_df.reset_index().fillna(0)

    # Write master tables
    long_out = f"{args.out_prefix}_by_pop_long.tsv"
    wide_out = f"{args.out_prefix}_by_pop.tsv"
    
    master_df.to_csv(long_out, sep='\t', index=False)
    wide_df.to_csv(wide_out, sep='\t', index=False)
    
    print(f"Success! Master files written to:\n  - {long_out}\n  - {wide_out}")

if __name__ == "__main__":
    main()