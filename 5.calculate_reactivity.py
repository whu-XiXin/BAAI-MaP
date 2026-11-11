# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import argparse
import sys

def pct_to_float(x):
    if pd.isna(x) or x == "" or x == ".":
        return np.nan
    try:
        return float(str(x).rstrip('%')) / 100.0
    except Exception:
        return np.nan

def compute_coefficient(group):
    n = len(group)
    if n < 10:
        return np.nan
    vals = group['reactivity.mutation'].sort_values(ascending=False)
    start = int(np.ceil(0.02 * n))
    end = int(np.ceil(0.10 * n))
    end = max(end, start + 1)
    return vals.iloc[start - 1 : end].mean()

def main():
    parser = argparse.ArgumentParser(
        description="End-to-end reactivity analysis with customizable sample groups",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True, help='Input file (TSV/VCF with GENE, Reads, VarFreq)')
    parser.add_argument('-o', '--output', required=True, help='Output final TSV')
    parser.add_argument('--treatment', nargs='+', required=True, help='Treatment group prefixes (e.g., CDI-2-1 CDI-2-2)')
    parser.add_argument('--control', nargs='+', required=True, help='Control group prefixes (e.g., Total-con-1 Total-con-2)')
    parser.add_argument('--chunksize', type=int, default=50000, help='Chunk size for Fisher test')
    args = parser.parse_args()

    treatment_prefixes = args.treatment
    control_prefixes = args.control

    print(f"Treatment groups: {treatment_prefixes}")
    print(f"Control groups: {control_prefixes}")

    # ========== Step 1: Load data ==========
    print("Loading input data...")
    df = pd.read_csv(args.input, sep='\t', low_memory=False)
    
    # Check GENE column exists
    if 'GENE' not in df.columns:
        raise ValueError("Column 'GENE' is required!")

    df = df[
        df['GENE'].notna() &
        (df['GENE'].str.strip() != '') &
        (df['GENE'].str.strip() != '.')
    ]

    df = df.reset_index(drop=True)

    # Filter genes with >=50 sites
    gene_counts = df['GENE'].value_counts()
    valid_genes = gene_counts[gene_counts >= 50].index
    df = df[df['GENE'].isin(valid_genes)].reset_index(drop=True)
    
    print(f"After GENE filtering: {len(df)} rows")

    # ========== Step 2: Build column names ==========
    varfreq_cols_treat = [f"{p}_VarFreq" for p in treatment_prefixes]
    varfreq_cols_ctrl  = [f"{p}_VarFreq" for p in control_prefixes]
    reads_cols_treat   = [f"{p}_Reads1" for p in treatment_prefixes] + [f"{p}_Reads2" for p in treatment_prefixes]
    reads_cols_ctrl    = [f"{p}_Reads1" for p in control_prefixes] + [f"{p}_Reads2" for p in control_prefixes]

    all_varfreq_cols = varfreq_cols_treat + varfreq_cols_ctrl
    all_reads_cols   = reads_cols_treat + reads_cols_ctrl

    # Validate columns exist
    missing_cols = [col for col in all_varfreq_cols + all_reads_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # ========== Step 3: Compute normalized reactivity ==========
    print("Computing normalized reactivity...")

    # Convert VarFreq
    for col in all_varfreq_cols:
        df[col] = df[col].apply(pct_to_float)

    # Average mutation frequencies
    df['mutation.average'] = df[varfreq_cols_treat].mean(axis=1, skipna=True)
    df['BG.mutation.average'] = df[varfreq_cols_ctrl].mean(axis=1, skipna=True)
    df['reactivity.mutation'] = df['mutation.average'] - df['BG.mutation.average']

    # Coefficient per GENE
    coeff_map = df.groupby('GENE').apply(compute_coefficient).to_dict()
    df['coefficient'] = df['GENE'].map(coeff_map)
    df['Reactivity_value'] = df['reactivity.mutation'] / df['coefficient']
    
    # 替换 +inf 为 10，-inf 为 NaN（后续 fillna(0) 会变成 0）
    df['Reactivity_value'] = (
        df['Reactivity_value']
        .replace(np.inf, 10)          # 正无穷 → 10
        .replace(-np.inf, np.nan)     # 负无穷 → NaN
        .fillna(0)                    # 所有 NaN（包括原生 NaN 和 -inf）→ 0
    )

    # ========== Step 4: Prepare Reads as integers ==========
    for col in all_reads_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)

    # ========== Step 5: Fisher exact test for each treatment ==========
    print("Performing Fisher exact tests...")

    # Initialize result columns
    for i, prefix in enumerate(treatment_prefixes):
        df[f'p_{prefix}_vs_con'] = 1.0
        df[f'significant_{prefix}'] = 'F'

    n_rows = len(df)
    chunk_size = args.chunksize

    for start in range(0, n_rows, chunk_size):
        end = min(start + chunk_size, n_rows)
        for idx in range(start, end):
            try:
                # Sum control reads
                con_ref = sum(df.at[idx, f"{p}_Reads1"] for p in control_prefixes)
                con_var = sum(df.at[idx, f"{p}_Reads2"] for p in control_prefixes)

                # For each treatment
                for prefix in treatment_prefixes:
                    t_ref = df.at[idx, f"{prefix}_Reads1"]
                    t_var = df.at[idx, f"{prefix}_Reads2"]

                    table = [[t_var, t_ref], [con_var, con_ref]]
                    _, p_val = fisher_exact(table, alternative='greater')

                    df.at[idx, f'p_{prefix}_vs_con'] = p_val
                    df.at[idx, f'significant_{prefix}'] = 'T' if p_val < 0.05 else 'F'
            except Exception:
                # Keep defaults on error
                continue

    # ========== Step 6: Reclassify reactivity ==========
    print("Reclassifying reactivity...")

    # -----------------------------
    # 1. Reactivity
    # -----------------------------
    reactivity_raw = df['Reactivity_value'].clip(lower=0) 
    conditions_raw = [
        (reactivity_raw >= 0.7),
        (reactivity_raw >= 0.3) & (reactivity_raw < 0.7),
        (reactivity_raw > 0) & (reactivity_raw < 0.3)
    ]
    choices = ['H', 'M', 'L']
    df['Reactivity'] = np.select(conditions_raw, choices, default='N')

    # -----------------------------
    # 2. Reactivity_Fisher
    # -----------------------------
    sig_cols = [f'significant_{p}' for p in treatment_prefixes]
    mask = (df[sig_cols] == 'T').all(axis=1)
    df['Reactivity_value_Fisher'] = np.where(mask, df['Reactivity_value'], 0.0)
    df['Reactivity_value_Fisher'] = df['Reactivity_value_Fisher'].clip(lower=0)

    conditions_fisher = [
        (df['Reactivity_value_Fisher'] >= 0.7),
        (df['Reactivity_value_Fisher'] >= 0.3) & (df['Reactivity_value_Fisher'] < 0.7),
        (df['Reactivity_value_Fisher'] > 0) & (df['Reactivity_value_Fisher'] < 0.3)
    ]
    df['Reactivity_Fisher'] = np.select(conditions_fisher, choices, default='N')

    # -----------------------------
    # 3. Enforce ratio < 1.5 → 'N'
    # -----------------------------
    df['ratio'] = np.where(
        df['BG.mutation.average'] == 0,
        np.inf,
        df['mutation.average'] / df['BG.mutation.average']
    )
    low_ratio_mask = df['ratio'] < 1.5
    df.loc[low_ratio_mask, 'Reactivity'] = 'N'
    df.loc[low_ratio_mask, 'Reactivity_Fisher'] = 'N'
    df.drop(columns=['ratio'], inplace=True)

    # ========== Save ==========
    print(f"Saving to: {args.output}")
    df.to_csv(args.output, sep='\t', index=False, na_rep='NA')
    print("Done!")

if __name__ == '__main__':
    main()
