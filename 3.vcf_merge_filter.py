#!/usr/bin/env python3

import pandas as pd
import argparse

def remove_blank_rows(df):
    """Remove all rows that are completely blank"""
    return df.dropna(how='all')

def calculate_sums_and_filter(df, threshold, columns_to_group):
    """
    Perform pairwise summation based on the provided list of column names, and filter rows where all group sums exceed the threshold.

    :param df: Input DataFrame
    :param threshold: The sum of each column group must be greater than this threshold
    :param columns_to_group: List of lists, each sublist contains two column names to be summed
    :return: Filtered DataFrame
    """
    # Flatten all column names
    flat_columns = [col for group in columns_to_group for col in group]

    # Check for missing columns
    missing_columns = set(flat_columns) - set(df.columns)
    if missing_columns:
        raise ValueError(f"Missing columns: {missing_columns}")

    # Convert columns to numeric, set non-convertible values to 0
    df[flat_columns] = df[flat_columns].apply(pd.to_numeric, errors='coerce').fillna(0)

    # Dynamically compute the sum for each group of columns
    for i, (col1, col2) in enumerate(columns_to_group):
        df[f'x{i+1}'] = df[col1] + df[col2]

    # Create a condition where all x columns must be greater than the threshold
    x_cols = [f'x{i+1}' for i in range(len(columns_to_group))]
    condition = (df[x_cols] > threshold).all(axis=1)

    # Apply the filter condition and remove temporary columns
    filtered_df = df[condition]
    filtered_df.drop(columns=x_cols, inplace=True)

    return filtered_df

def process_chunk(chunk, output_file, threshold, columns_to_group, write_header):
    """Process each data chunk and append to the output file"""
    chunk = remove_blank_rows(chunk)
    chunk = calculate_sums_and_filter(chunk, threshold, columns_to_group)
    chunk.to_csv(output_file, sep='\t', index=False, mode='a', header=write_header)

def main(input_file, output_file, threshold=100, column_names=None, chunk_size=10**6):
    """Main function: read in chunks, process, and save results"""
    if column_names is None:
        column_names = []

    # Validate that an even number of column names is provided for pairing
    if len(column_names) % 2 != 0:
        raise ValueError("An even number of column names must be provided for pairing")

    # Group column names into pairs
    columns_to_group = [column_names[i:i+2] for i in range(0, len(column_names), 2)]

    reader = pd.read_csv(input_file, sep='\t', chunksize=chunk_size)

    first_chunk = True
    for chunk in reader:
        process_chunk(chunk, output_file, threshold, columns_to_group, first_chunk)
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
    parser.add_argument('-name', nargs='+', required=True,
                        help="List of column names to process in pairs, e.g.: "
                             "RMR-1_Reads1 RMR-1_Reads2 RMR-2_Reads1 RMR-2_Reads2 ...")

    args = parser.parse_args()

    main(args.input, args.out, args.threshold, args.name)
