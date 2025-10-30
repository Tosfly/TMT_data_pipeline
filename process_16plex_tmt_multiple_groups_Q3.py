#!/usr/bin/env python3
"""
Enhanced 16-plex TMT multi-group processor

Additions while preserving original behaviors:
1) Keep creating a 'Genename' column from the Description (GN= tag).
2) Remove rows where the 'accession' column contains 'Reverse_' or 'contaminant_' (case-insensitive).
3) Exclude any TMT channel not assigned a Group in Sample_info.csv from the final output.
4) For each post hoc comparison 'A vs B', compute mean(A)/mean(B) and mean(B)/mean(A) as 'A/B' and 'B/A'.
5) Add a flag to name the output file: -o / --output.
6) Add a flag --spec-count-min to drop rows with 'spec count' < threshold.
7) Add a flag --q3-normalize to apply Q3 normalization after data cleanup.
"""

import argparse
from itertools import combinations
import numpy as np
import pandas as pd
from scipy import stats




def remove_avg_columns(df):
    """Drop any columns whose name contains 'avg' (case-insensitive),
    e.g., 'avg int m/z_129.13779'. """
    cols_to_drop = [c for c in df.columns if 'avg' in str(c).lower()]
    if cols_to_drop:
        df = df.drop(columns=cols_to_drop)
    return df


def q3_normalize(df, sample_cols):
    """Q3 normalize each sample to the median Q3 across samples."""
    if not sample_cols:
        return df
    
    # Extract sample data
    df_samples = df[sample_cols].copy()
    
    # Convert to numeric
    df_samples = df_samples.apply(pd.to_numeric, errors='coerce')
    
    # Treat zeros as missing for Q3 calc
    q3 = df_samples.replace(0, np.nan).quantile(0.75, axis=0, numeric_only=True)
    target = np.nanmedian(q3.values)
    
    # Handle edge cases
    if np.isnan(target) or target == 0:
        print("Warning: Q3 normalization skipped (invalid target value)")
        return df
    
    scale = target / q3
    # Replace inf/nan in scale with 1.0 (no scaling for those samples)
    scale = scale.replace([np.inf, -np.inf], 1.0).fillna(1.0)
    
    # Apply scaling
    df[sample_cols] = df_samples * scale
    
    print(f"Q3 normalization applied to {len(sample_cols)} samples (target Q3: {target:.2f})")
    return df


def GN_select(names):
    """Extract gene names from description column as in original script."""
    new = []
    for i in names:
        if isinstance(i, str) and ('GN=' in i):
            s = i.split('GN=')[1]
            e = s.find(' ')
            new.append(s[:e] if e > -1 else s)
        else:
            new.append(i)
    return new


def add_genename_column(df):
    # Match original column name casing: 'Genename'
    if 'Genename' not in df.columns:
        # Try common description fields
        for desc_col in ['Description', 'description', 'Protein Descriptions', 'Protein.Descriptions']:
            if desc_col in df.columns:
                df.insert(0, 'Genename', GN_select(df[desc_col].tolist()))
                break
        else:
            # Fallback: if a GeneNames column already exists, reuse it
            if 'GeneNames' in df.columns:
                df.insert(0, 'Genename', df['GeneNames'].astype(str))
            else:
                df.insert(0, 'Genename', 'NA')
    return df


def drop_reverse_and_contaminants(df):
    """Remove rows with 'Reverse_' or 'contaminant_' in 'accession' column, case-insensitive."""
    acc_col = None
    for cand in ['accession', 'Accession', 'Protein IDs', 'Protein.IDs']:
        if cand in df.columns:
            acc_col = cand
            break
    if acc_col:
        acc = df[acc_col].astype(str)
        mask = ~(acc.str.contains('Reverse_', case=False, na=False) | acc.str.contains('contaminant_', case=False, na=False))
        df = df[mask].copy()
    return df


def build_column_mapping(df, sample_info):
    """Map proList intensity columns to 'Sample_<id>_<group>_<sex>', keep only assigned groups."""
    info = sample_info.copy()
    info['Group'] = info['Group'].astype(str).str.strip()
    info = info[(info['Group'].notna()) & (info['Group'] != '')]

    # Identify raw TMT intensity columns
    raw_tmt_cols = [c for c in df.columns if 'norm int' in c]

    col_mapping = {}
    unmapped = []
    for c in raw_tmt_cols:
        unmapped.append(c)

    for _, row in info.iterrows():
        tmt_channel = str(row['TMT_channels'])
        sample_id = str(row['Sample_ids']).rstrip('.0')
        group = row['Group']
        sex = str(row.get('Sex', 'Unknown'))
        target = f"Sample_{sample_id}_{group}_{sex}"
        # Find the original column that contains this TMT channel
        for col in raw_tmt_cols:
            if tmt_channel in col:
                col_mapping[col] = target
                if col in unmapped:
                    unmapped.remove(col)
                break

    return col_mapping, unmapped, info


def compute_statistics_and_ratios(df, sample_info):
    sample_cols = [c for c in df.columns if c.startswith('Sample_')]
    if not sample_cols:
        return df

    groups = list(pd.unique(sample_info['Group']))
    group_cols = {g: [c for c in sample_cols if f"_{g}_" in c] for g in groups}

    # One-way ANOVA
    anova_pvals = []
    for idx in df.index:
        group_data = []
        for g in groups:
            vals = pd.to_numeric(df.loc[idx, group_cols[g]], errors='coerce').values if group_cols[g] else np.array([])
            vals = vals[~np.isnan(vals)]
            if len(vals) > 0:
                group_data.append(vals)
        if len(group_data) >= 2:
            try:
                _, p = stats.f_oneway(*group_data)
            except Exception:
                p = np.nan
        else:
            p = np.nan
        anova_pvals.append(p)
    df['One-Way ANOVA'] = anova_pvals

    # Post hoc and ratios
    for g1, g2 in combinations(groups, 2):
        comp_col = f"{g1} vs {g2}"
        pvals, r12_list, r21_list = [], [], []
        for idx in df.index:
            v1 = pd.to_numeric(df.loc[idx, group_cols[g1]], errors='coerce').values if group_cols[g1] else np.array([])
            v2 = pd.to_numeric(df.loc[idx, group_cols[g2]], errors='coerce').values if group_cols[g2] else np.array([])
            v1 = v1[~np.isnan(v1)]
            v2 = v2[~np.isnan(v2)]

            if len(v1) > 0 and len(v2) > 0:
                try:
                    _, p = stats.ttest_ind(v1, v2, equal_var=False)
                except Exception:
                    p = np.nan
            else:
                p = np.nan
            pvals.append(p)

            m1 = np.nan if len(v1) == 0 else float(np.mean(v1))
            m2 = np.nan if len(v2) == 0 else float(np.mean(v2))
            r12 = (m1 / m2) if (m1 == m1 and m2 == m2 and m2 != 0) else np.nan
            r21 = (m2 / m1) if (m1 == m1 and m2 == m2 and m1 != 0) else np.nan
            r12_list.append(r12)
            r21_list.append(r21)

        df[comp_col] = pvals
        insert_idx = df.columns.get_loc(comp_col) + 1
        df.insert(insert_idx, f"{g1}/{g2}", r12_list)
        df.insert(insert_idx + 1, f"{g2}/{g1}", r21_list)

    return df


def main():
    parser = argparse.ArgumentParser(description="Process 16-plex TMT data across multiple groups.")
    parser.add_argument('--prolist', required=True, help='Path to proList.csv')
    parser.add_argument('--sampleinfo', required=True, help='Path to Sample_info.csv')
    parser.add_argument('-o', '--output', default='processed_TMT_data.csv', help='Output CSV filename')
    parser.add_argument('--spec-count-min', type=float, default=None,
                        help="Drop rows with 'spec count' less than this threshold")
    parser.add_argument('--q3-normalize', action='store_true',
                        help='Apply Q3 normalization after data cleanup and before statistics')
    args = parser.parse_args()

    # Load proList
    df = pd.read_csv(args.prolist)
    # Remove any average/intensity summary columns not in Sample_info
    df = remove_avg_columns(df)

    # Keep original behaviors
    df = add_genename_column(df)
    df = drop_reverse_and_contaminants(df)

    # Optional spec count filter
    if args['spec_count_min'] if isinstance(args, dict) else args.spec_count_min is not None:
        thresh = args.spec_count_min if not isinstance(args, dict) else args['spec_count_min']
        # Find spec count column variants
        spec_col = None
        for cand in ['spec count', 'Spec count', 'Spec Count', 'Spectral Count', 'Spectral.Count']:
            if cand in df.columns:
                spec_col = cand
                break
        if spec_col is not None:
            df = df[pd.to_numeric(df[spec_col], errors='coerce').fillna(0) >= thresh].copy()

    # Build mapping and drop unmapped intensity columns
    sample_info = pd.read_csv(args.sampleinfo)
    col_mapping, unmapped, info_used = build_column_mapping(df, sample_info)
    if col_mapping:
        df = df.rename(columns=col_mapping)
    if unmapped:
        df = df.drop(columns=unmapped)

    # Optional Q3 normalization (after data cleanup, before statistics)
    if args.q3_normalize:
        sample_cols = [c for c in df.columns if c.startswith('Sample_')]
        df = q3_normalize(df, sample_cols)

    # Compute stats and ratios
    df = compute_statistics_and_ratios(df, info_used)

    # Final touch
    df = df.fillna('NA')
    df.to_csv(args.output, index=False)

    # Basic summary
    groups = list(pd.unique(info_used['Group']))
    print(f"Saved: {args.output}")
    print(f"Rows: {len(df)} | Groups: {', '.join(groups)}")

if __name__ == "__main__":
    main()
