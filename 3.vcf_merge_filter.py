#!/usr/bin/env python3

import pandas as pd
import argparse


def remove_blank_rows(df):
    """Remove all rows that are completely empty (all NaN)"""
    return df.dropna(how='all')


def calculate_sums_and_filter(df, threshold):
    """
    Calculate the sum of column pairs:
        - x1: sum of columns 5 and 6 ('RMR-1_Reads1' and 'RMR-1_Reads2')
        - x2: sum of columns 7 and 8 ('RMR-2_Reads1' and 'RMR-2_Reads2')
        - x3: sum of columns 9 and 10 ('RMR-con-1_Reads1' and 'RMR-con-1_Reads2')
        - x4: sum of columns 11 and 12 ('RMR-con-2_Reads1' and 'RMR-con-2_Reads2')

    Only keep rows where x1, x2, x3, and x4 are all greater than the given threshold.
    """

    # List of relevant columns
    columns_to_sum = [
        'RMR-1_Reads1', 'RMR-1_Reads2',
        'RMR-2_Reads1', 'RMR-2_Reads2',
        'RMR-con-1_Reads1', 'RMR-con-1_Reads2',
        'RMR-con-2_Reads1', 'RMR-con-2_Reads2'
    ]

    # Check for missing columns
    missing_columns = set(columns_to_sum) - set(df.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    # Convert to numeric, coerce errors to NaN and fill with 0
    df[columns_to_sum] = df[columns_to_sum].apply(pd.to_numeric, errors='coerce').fillna(0)

    # Compute sums
    df['x1'] = df['RMR-1_Reads1'] + df['RMR-1_Reads2']
    df['x2'] = df['RMR-2_Reads1'] + df['RMR-2_Reads2']
    df['x3'] = df['RMR-con-1_Reads1'] + df['RMR-con-1_Reads2']
    df['x4'] = df['RMR-con-2_Reads1'] + df['RMR-con-2_Reads2']

    # Filter based on threshold
    df = df[(df['x1'] > threshold) & (df['x2'] > threshold) &
            (df['x3'] > threshold) & (df['x4'] > threshold)]

    # Remove temporary columns
    df.drop(columns=['x1', 'x2', 'x3', 'x4'], inplace=True)

    return df


def process_chunk(chunk, output_file, threshold, header=False):
    """Process a chunk: remove blank rows, filter by column sums, and append to output file"""
    chunk = remove_blank_rows(chunk)
    chunk = calculate_sums_and_filter(chunk, threshold)
    chunk.to_csv(output_file, sep='\t', index=False, mode='a', header=header)


def main(input_file, output_file, threshold=100, chunk_size=10 ** 6):
    """Main function: read input in chunks, process each one, and save the result"""
    reader = pd.read_csv(input_file, sep='\t', chunksize=chunk_size)

    first_chunk = True
    for chunk in reader:
        process_chunk(chunk, output_file, threshold, header=first_chunk)
        first_chunk = False

    print(f"Processed file has been saved to {output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Process a tab-separated file by removing blank rows and filtering based on column sums."
    )
    parser.add_argument('-input', required=True, help="Input tab-separated file to be processed")
    parser.add_argument('-out', required=True, help="Output file to save the processed result")
    parser.add_argument('-d', '--threshold', type=int, default=100,
                        help="Threshold value for filtering rows (default: 100)")

    args = parser.parse_args()

    main(args.input, args.out, args.threshold)
