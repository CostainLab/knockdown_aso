#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    # define arguments for script
    parser = argparse.ArgumentParser(
        description="Tally sample coverage proportions per population."
        )

    parser.add_argument(
        "--coverage",
        required = True,
        help="Input sample_coverage.tsv file")
    
    parser.add_argument(
        "--popdata",
        required = True,
        help="1000G population file"
    )

    parser.add_argument(
        "--gene",
        required=True,
        help="Name of gene"
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output tsv file name"
    )

    args = parser.parse_args()

    # load coverage file
    cov_df = pd.read_csv(args.coverage, sep='\t')

    # load population data file
    popdata_df = pd.read_csv(
        args.popdata, 
        sep=' ', 
        usecols=['SampleID', 'Population', 'Superpopulation']
    )
    # rename columns
    popdata_df = popdata_df.rename(columns={
        'SampleID': 'sample_id', 
        'Population': 'pop',
        'Superpopulation': 'super_pop'
    })

    # merge datasets
    cov_pop_df = pd.merge(cov_df, popdata_df, on='sample_id', how='left')
    # if any sample ids were missing from popdata, fills with Unknown
    cov_pop_df['super_pop'] = cov_pop_df['super_pop'].fillna('Unknown')

    # group by super_pop and count by coverage status
    tally = cov_pop_df.groupby(['super_pop', 'coverage_status']).size().unstack(fill_value=0)

    # ensure both columns exist (may have 100% covered or 0% covered gene)
    for status in ['Covered', 'Not_Covered']:
        if status not in tally.columns:
            tally[status] = 0
    
    # calculate totals and proportions
    tally['Total_Samples'] = tally['Covered'] + tally['Not_Covered']
    tally['Proportion_Covered'] = tally['Covered'] / tally['Total_Samples']
    tally['Gene'] = args.gene

    # clean df structure for output
    tally = tally.reset_index()

    # reorder columns
    final_cols = [
        'Gene', 'super_pop', 'Covered', 'Not_Covered', 
        'Total_Samples', 'Proportion_Covered'
    ]
    tally = tally[final_cols]
    tally.columns = tally.columns.str.lower()

    # sort by highest ppln coverate
    tally = tally.sort_values(by='proportion_covered', ascending=False)

    # write to file
    tally.to_csv(args.out, sep='\t', index=False)
    print(f"Coverage summary for {args.gene} written to {args.out}")

if __name__ == "__main__":
    main()
