#!/usr/bin/env python3
"""
process_16plex_tmt_multiple_groups_Volc_plots.py

Extended processing of 16‑plex TMT datasets.

Key additions
-------------
1.  Sample_info.csv **no longer** contains a “Sex” column.  Intensity columns are renamed
    to:  `Sample_<SampleID>_<Group>`.

2.  Calculates one‑way ANOVA p‑values across **all** groups per protein.

3.  For every pair of groups:
      * two‑sample (Welch) t‑test p‑value
      * –log10(p)  (for volcano)
      * mean ratio  (mean(group1)/mean(group2))
      * log2 ratio

4.  Produces:
      * `<prefix>_all_info.csv`       –   complete table
      * `<prefix>_vol_<g1>_vs_<g2>.csv`  – volcano‑ready subset
      * `<prefix>_vol_<g1>_vs_<g2>.png`  – volcano plot

Usage
-----
python process_16plex_tmt_multiple_groups_modified.py \
       --prolist proList.csv \
       --sampleinfo Sample_info.csv \
       --output  result_prefix

Dependencies
------------
pandas, numpy, scipy, matplotlib

"""
import argparse, os, sys
import pandas as pd, numpy as np
from itertools import combinations
from scipy import stats
import matplotlib.pyplot as plt

def GN_select(names):
    out = []
    for x in names:
        if 'GN=' in x:
            s = x.split('GN=')[1]
            e = s.find(' ')
            out.append(s[:e])
        else:
            out.append(x)
    return out

def build_column_mapping(df, sample_info):
    """
    Creates {old_name : new_name} dict and group→columns mapping.
    """
    col_map = {}
    group_cols = {}
    for _, row in sample_info.iterrows():
        chan   = str(row['TMT_channels'])
        sample = str(row['Sample_ids'])
        group  = str(row['Group'])
        new_name = f"Sample_{sample}_{group}"
        # locate channel‑specific column inside df
        cand = [c for c in df.columns if chan in c and 'norm int' in c.lower()]
        if not cand:
            cand = [c for c in df.columns if chan in c]   # fallback
        if not cand:
            sys.stderr.write(f"⚠️  Channel {chan} not found – skipping\n")
            continue
        old = cand[0]
        col_map[old] = new_name
        group_cols.setdefault(group, []).append(new_name)
    return col_map, group_cols

def calc_anova(df, group_cols):
    pvals = []
    for idx in df.index:
        samples = []
        for cols in group_cols.values():
            vals = pd.to_numeric(df.loc[idx, cols], errors='coerce')
            vals = vals.dropna()
            if len(vals):
                samples.append(vals.values)
        # need ≥2 groups with data
        if sum(len(s) > 0 for s in samples) >= 2:
            try:
                _, p = stats.f_oneway(*samples)
            except Exception:
                p = np.nan
        else:
            p = np.nan
        pvals.append(p)
    return pvals

def calc_pairwise(df, group_cols):
    group_pairs = list(combinations(group_cols.keys(), 2))
    for g1, g2 in group_pairs:
        col_p      = f"{g1} vs {g2} p"
        col_logp   = f"-log10({g1} vs {g2})"
        col_ratio  = f"mean({g1})/mean({g2})"
        col_log2fc = f"log2FC {g1}/{g2}"
        pvals, logp, ratio, log2fc = [], [], [], []
        for idx in df.index:
            v1 = pd.to_numeric(df.loc[idx, group_cols[g1]], errors='coerce').dropna()
            v2 = pd.to_numeric(df.loc[idx, group_cols[g2]], errors='coerce').dropna()
            # Welch's t‑test if both groups have data
            if len(v1) and len(v2):
                try:
                    _, p = stats.ttest_ind(v1, v2, equal_var=False, nan_policy='omit')
                except Exception:
                    p = np.nan
                r  = v1.mean() / v2.mean() if v2.mean() != 0 else np.nan
            else:
                p, r = np.nan, np.nan
            pvals.append(p)
            logp.append(-np.log10(p) if p and p>0 else np.nan)
            ratio.append(r)
            log2fc.append(np.log2(r) if r and r>0 else np.nan)
        df[col_p]      = pvals
        df[col_logp]   = logp
        df[col_ratio]  = ratio
        df[col_log2fc] = log2fc
    return group_pairs

def make_volcano_files(df, group_pairs, prefix):
    keep_cols_base = ['GeneName', 'description']
    for g1, g2 in group_pairs:
        pcol   = f"{g1} vs {g2} p"
        logp   = f"-log10({g1} vs {g2})"
        logfc  = f"log2FC {g1}/{g2}"
        subset = df[keep_cols_base + [pcol, logp, logfc]].copy()
        fout   = f"{prefix}_vol_{g1}_vs_{g2}.csv"
        subset.to_csv(fout, index=False)
        # volcano plot
        png    = f"{prefix}_vol_{g1}_vs_{g2}.png"
        x = subset[logfc]
        y = subset[logp]
        plt.figure(figsize=(6,6))
        plt.scatter(x, y, s=8, alpha=0.6)
        plt.axvline(0, linestyle='--', linewidth=1)
        plt.xlabel(f"log2 mean({g1}) / mean({g2})")
        plt.ylabel(f"-log10 p‑value")
        plt.title(f"Volcano: {g1} vs {g2}")
        plt.tight_layout()
        plt.savefig(png, dpi=300)
        plt.close()

def main():
    ap = argparse.ArgumentParser(description="Process 16‑plex TMT data with ANOVA / volcano")
    ap.add_argument("--prolist", required=True, help="proList input CSV")
    ap.add_argument("--sampleinfo", required=True, help="Sample_info CSV (TMT_channels, Sample_ids, Group)")
    ap.add_argument("--output", default="processed", help="Output prefix (default: processed)")
    args = ap.parse_args()

    print("Loading data…")
    df = pd.read_csv(args.prolist)
    # drop 'avg' columns
    df = df.drop(columns=[c for c in df.columns if 'avg' in c.lower()], errors='ignore')

    # gene names
    if 'description' in df.columns and 'GeneName' not in df.columns:
        df.insert(0, 'GeneName', GN_select(df['description'].astype(str)))

    sample_info = pd.read_csv(args.sampleinfo)
    col_map, group_cols = build_column_mapping(df, sample_info)
    print(f"Renaming {len(col_map)} intensity columns…")
    df = df.rename(columns=col_map)
    sample_cols = list(col_map.values())

    print("Calculating one‑way ANOVA p‑values…")
    df['One‑Way ANOVA'] = calc_anova(df, group_cols)

    print("Pairwise post‑hoc & ratios…")
    group_pairs = calc_pairwise(df, group_cols)

    # save full table
    prefix = os.path.splitext(args.output)[0]
    all_info_path = f"{prefix}_all_info.csv"
    df.to_csv(all_info_path, index=False)
    print(f"Saved {all_info_path}")

    print("Generating volcano csv/plots…")
    make_volcano_files(df, group_pairs, prefix)
    print("Done.")

if __name__ == "__main__":
    main()
