#!/usr/bin/env python3

import pandas as pd
import argparse

def load_and_rename(file_path, prefix):
    """Load a file, keep the first 7 columns, and rename the last 3 columns to {prefix}_{col}"""
    df = pd.read_csv(file_path, sep='\t')
    df = df.iloc[:, :7]  # Take the first 7 columns
    new_columns = list(df.columns[:4]) + [f"{prefix}_{col}" for col in df.columns[4:]]
    df.columns = new_columns
    return df

def merge_dataframes(dfs):
    """Merge multiple DataFrames based on the first four columns (outer join)"""
    merged_df = dfs[0]
    for df in dfs[1:]:
        key_columns = list(merged_df.columns[:4])  # Columns to merge on: Chrom, Position, Ref, Cons
        merged_df = pd.merge(merged_df, df, on=key_columns, how='outer')
    return merged_df

def reorder_columns(merged_df, prefixes):
    """Reorder columns according to the order of prefixes provided via -n argument"""
    fixed_columns = list(merged_df.columns[:4])  # These are: Chrom, Position, Ref, Cons
    desired_columns = fixed_columns.copy()

    # Append columns corresponding to each prefix in the given order
    for prefix in prefixes:
        for col in merged_df.columns:
            if col.startswith(f"{prefix}_"):
                desired_columns.append(col)

    # Add remaining columns that were not matched (if any)
    remaining_cols = [col for col in merged_df.columns if col not in desired_columns]
    desired_columns += remaining_cols

    return merged_df[desired_columns]

def main():
    parser = argparse.ArgumentParser(description="Merge multiple VCF-like files by renaming and aligning columns.")
    parser.add_argument('-input', nargs='+', required=True, help="Input TSV/VCF files (no header comment lines)")
    parser.add_argument('-out', required=True, help="Output merged TSV file")
    parser.add_argument('-n', nargs='+', required=True, help="Prefixes for each input file (same number as input files)")

    args = parser.parse_args()

    if len(args.input) != len(args.n):
        raise ValueError("The number of input files must match the number of prefixes.")

    # Step 1: Load and rename each input file
    dfs = []
    for file, prefix in zip(args.input, args.n):
        df = load_and_rename(file, prefix)
        dfs.append(df)

    # Step 2: Merge all DataFrames
    merged_df = merge_dataframes(dfs)

    # Step 3: Reorder columns based on the prefix order from -n
    merged_df = reorder_columns(merged_df, args.n)

    # Step 4: Replace missing values with empty strings
    merged_df.fillna('', inplace=True)

    # Step 5: Save the merged result to a TSV file
    merged_df.to_csv(args.out, sep='\t', index=False)

    print(f"Merged file has been saved to {args.out}")

if __name__ == '__main__':
    main()
