#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations
import argparse
import sys

def GN_select(names):
    """Extract gene names from description column"""
    new = []
    for i in names:
        if i.find('GN=') > 0:
            s = i.split('GN=')[1]
            e = s.find(' ')
            new.append(s[:e])
        else:
            new.append(i)
    return new

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process 16plex TMT data')
    parser.add_argument('--prolist', required=True, help='Path to proList.csv')
    parser.add_argument('--sampleinfo', required=True, help='Path to Sample_info.csv')
    parser.add_argument('--output', default='processed_TMT_data.csv', help='Output filename (default: processed_TMT_data.csv)')
    
    args = parser.parse_args()
    
    # 1. Load the main input file
    print(f"Loading {args.prolist}...")
    df = pd.read_csv(args.prolist)
    
    # 2. Remove columns containing "avg"
    print("Removing avg columns...")
    avg_cols = [col for col in df.columns if 'avg' in col]
    df = df.drop(columns=avg_cols)
    
    # 3. Remove contaminant and reverse rows
    print("Removing contaminant and reverse rows...")
    df = df[~df['accession'].str.contains('contaminant_', na=False)]
    df = df[~df['accession'].str.contains('Reverse_', na=False)]
    
    # 4. Add GeneName column
    print("Extracting gene names...")
    gene_names = GN_select(df['description'].tolist())
    df.insert(0, 'GeneName', gene_names)
    
    # 5. Load sample info and rename TMT columns
    print(f"Loading {args.sampleinfo}...")
    sample_info = pd.read_csv(args.sampleinfo)
    
    # Create mapping dictionary for column renaming
    col_mapping = {}
    for _, row in sample_info.iterrows():
        # Find the corresponding column in proList
        tmt_channel = row['TMT_channels']
        # Look for columns containing the TMT channel value
        for col in df.columns:
            if 'norm int' in col and tmt_channel in col:
                new_name = f"Sample_{row['Sample_ids']}_{row['Group']}_{row['Sex']}"
                col_mapping[col] = new_name
                break
    
    # Rename columns
    print("Renaming TMT columns...")
    df = df.rename(columns=col_mapping)
    
    # 6. Get renamed sample columns for statistics
    sample_cols = [col for col in df.columns if col.startswith('Sample_')]
    
    # Extract groups from sample columns
    groups = sample_info['Group'].unique()
    
    # Create group-to-columns mapping
    group_cols = {}
    for group in groups:
        group_cols[group] = [col for col in sample_cols if f"_{group}_" in col]
    
    # 7. Calculate One-Way ANOVA
    print("Calculating One-Way ANOVA...")
    anova_pvals = []
    
    for idx in df.index:
        # Get values for each group
        group_data = []
        for group in groups:
            values = df.loc[idx, group_cols[group]].values
            # Convert to numeric and remove NaN values
            values = pd.to_numeric(values, errors='coerce')
            values = values[~np.isnan(values)]
            if len(values) > 0:
                group_data.append(values)
        
        # Perform ANOVA if we have at least 2 groups with data
        if len(group_data) >= 2:
            try:
                f_stat, p_val = stats.f_oneway(*group_data)
                anova_pvals.append(p_val)
            except:
                anova_pvals.append(np.nan)
        else:
            anova_pvals.append(np.nan)
    
    df['One-Way ANOVA'] = anova_pvals
    
    # 8. Calculate post hoc p-values for each pair of groups
    print("Calculating post hoc comparisons...")
    
    # Generate all possible pairs of groups
    group_pairs = list(combinations(groups, 2))
    
    for group1, group2 in group_pairs:
        col_name = f"{group1} vs {group2}"
        posthoc_pvals = []
        
        for idx in df.index:
            # Get values for each group
            values1 = df.loc[idx, group_cols[group1]].values
            values2 = df.loc[idx, group_cols[group2]].values
            
            # Convert to numeric and remove NaN values
            values1 = pd.to_numeric(values1, errors='coerce')
            values2 = pd.to_numeric(values2, errors='coerce')
            values1 = values1[~np.isnan(values1)]
            values2 = values2[~np.isnan(values2)]
            
            # Perform t-test if both groups have data
            if len(values1) > 0 and len(values2) > 0:
                try:
                    t_stat, p_val = stats.ttest_ind(values1, values2)
                    posthoc_pvals.append(p_val)
                except:
                    posthoc_pvals.append(np.nan)
            else:
                posthoc_pvals.append(np.nan)
        
        df[col_name] = posthoc_pvals
    
    # 9. Fill empty cells with 'NA'
    print("Filling empty cells with 'NA'...")
    df = df.fillna('NA')
    
    # 10. Save the processed data
    print(f"Saving processed data to {args.output}...")
    df.to_csv(args.output, index=False)
    
    print("Processing complete!")
    print(f"Original rows: {len(pd.read_csv(args.prolist))}")
    print(f"Processed rows: {len(df)}")
    print(f"Groups found: {', '.join(groups)}")
    print(f"Post hoc comparisons: {', '.join([f'{p[0]} vs {p[1]}' for p in group_pairs])}")

if __name__ == "__main__":
    main()
