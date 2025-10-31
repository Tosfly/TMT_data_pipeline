#!/usr/bin/env python3
from __future__ import annotations
"""
t-SNE from wide matrix with headers like 'Sample_556_Control'.

Features are rows. Samples are columns named with a pattern that includes an ID and a group,
for example: 'Sample_556_Control' where group='Control' and ID='556'.

Main steps:
1) Read CSV.
2) Parse sample headers -> sample_id, group.
3) Optional Q3 normalization (per sample).
4) Select top % features using one-way ANOVA across groups (highest F).
5) Run t-SNE and make a scatter plot with labels (IDs) and colors by group.

Outputs:
- PNG plot
- CSV of t-SNE coordinates
- TXT list of selected features (when feature selection applied)

Dependencies:
  pip install pandas numpy scipy scikit-learn matplotlib
"""

import argparse
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import f_oneway
from sklearn.manifold import TSNE
from matplotlib.cm import get_cmap


def parse_args():
    p = argparse.ArgumentParser(description="t-SNE plot from matrix with sample headers like 'Sample_556_Control'.")
    p.add_argument("csv", type=str, help="Input CSV file. Rows are features, columns are samples.")
    p.add_argument("--index-col", type=str, default=None,
                   help="Name of the feature ID column. If omitted, the first column will be used if it is non-numeric.")
    p.add_argument("--q3", action="store_true", help="Apply Q3 normalization per sample.")
    p.add_argument("--top-percent", type=float, default=100.0,
                   help="Top percent of features to keep by ANOVA F-stat (0-100). Default 100 means use all features.")
    p.add_argument("--perplexity", type=float, default=30.0, help="t-SNE perplexity.")
    p.add_argument("--n-iter", type=int, default=1000, help="t-SNE iterations.")
    p.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    p.add_argument("--out-prefix", type=str, default="tsne",
                   help="Prefix for outputs: <prefix>_plot.png, <prefix>_coords.csv, <prefix>_features.txt")
    p.add_argument("--separate-plots", action="store_true",
                   help="Create separate comparison plots for each group pair in addition to the main plot")
    return p.parse_args()


HEADER_REGEXES = [
    # Match pattern like "Sample_568_Experimental control" - capture ID and everything after it as group
    re.compile(r"^Sample[_\- ]+(\d+)[_\- ]+(.+)$", re.IGNORECASE),
    # More flexible: any text with numeric ID followed by underscore/hyphen and group name (with spaces)
    re.compile(r".*?(\d+)[_\- ]+(.+)$"),
]


def parse_header(name: str):
    """Return (sample_id:str, group:str) if matched, else None."""
    for rx in HEADER_REGEXES:
        m = rx.match(name)
        if m:
            sid, grp = m.group(1), m.group(2)
            return sid, grp
    return None


def choose_index_column(df: pd.DataFrame, index_col_opt: str | None):
    if index_col_opt:
        if index_col_opt not in df.columns:
            sys.exit(f"Index column '{index_col_opt}' not found in CSV.")
        df = df.set_index(index_col_opt)
        return df
    # Auto: if first column is mostly non-numeric, use it as index
    first_col = df.columns[0]
    col_vals = df[first_col]
    numeric_ratio = pd.to_numeric(col_vals, errors="coerce").notna().mean()
    if numeric_ratio < 0.5:
        df = df.set_index(first_col)
    return df


def q3_normalize(df_samples: pd.DataFrame) -> pd.DataFrame:
    """Q3 normalize each sample to the median Q3 across samples."""
    # Treat zeros as missing for Q3 calc
    q3 = df_samples.replace(0, np.nan).quantile(0.75, axis=0, numeric_only=True)
    target = np.nanmedian(q3.values)
    
    # Handle edge cases
    if np.isnan(target) or target == 0:
        print("Warning: Q3 normalization skipped (invalid target value)")
        return df_samples
    
    scale = target / q3
    # Replace inf/nan in scale with 1.0 (no scaling for those samples)
    scale = scale.replace([np.inf, -np.inf], 1.0).fillna(1.0)
    
    return df_samples * scale


def anova_f_vector(df: pd.DataFrame, groups: dict[str, list[str]]) -> pd.Series:
    """Compute one-way ANOVA F-stat for each feature across groups."""
    f_vals = []
    index = df.index
    # Pre-assemble arrays of column blocks for speed
    group_cols = {g: df[cols].values for g, cols in groups.items()}
    for i in range(len(index)):
        arrays = []
        for g, mat in group_cols.items():
            vals = mat[i, :]
            vals = vals[~np.isnan(vals)]
            if vals.size >= 2:
                arrays.append(vals)
        if len(arrays) >= 2:
            try:
                f, _ = f_oneway(*arrays)
                # Handle NaN F-statistic
                if np.isnan(f):
                    f = 0.0
            except Exception:
                f = 0.0
        else:
            f = 0.0
        f_vals.append(f)
    return pd.Series(f_vals, index=index, name="F")


def pick_top_features_by_f(df: pd.DataFrame, f_series: pd.Series, top_percent: float):
    top_percent = float(top_percent)
    if top_percent >= 100.0:
        return df, f_series
    if top_percent <= 0.0:
        sys.exit("top-percent must be > 0.")
    k = max(1, int(np.ceil(len(f_series) * (top_percent / 100.0))))
    keep = f_series.sort_values(ascending=False).head(k).index
    return df.loc[keep], f_series.loc[keep]


def make_color_map(group_names):
    """Create a color map for groups using tab20 for better color variety."""
    try:
        # Use tab20 for up to 20 distinct colors
        base = get_cmap("tab20")
        n_colors = 20
    except:
        # Fallback to tab10 if tab10 not available
        base = get_cmap("tab10")
        n_colors = 10
    
    uniq = list(group_names)
    colors = {g: base(i % n_colors) for i, g in enumerate(uniq)}
    return colors


def make_marker_map(group_names):
    """Create a marker map for groups to visually distinguish them."""
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h', 'H', '+', 'x', 'X', 'd']
    uniq = list(group_names)
    marker_map = {g: markers[i % len(markers)] for i, g in enumerate(uniq)}
    return marker_map


def create_comparison_plots(plot_df, color_map, text_offset, out_prefix):
    """Create separate plots for each pair of groups."""
    from itertools import combinations
    
    groups = sorted(plot_df["group"].unique())
    pairs = list(combinations(groups, 2))
    
    print(f"Creating {len(pairs)} comparison plots...")
    
    for g1, g2 in pairs:
        subset = plot_df[plot_df["group"].isin([g1, g2])]
        
        plt.figure(figsize=(7, 6), dpi=150)
        
        for g, sub in subset.groupby("group"):
            plt.scatter(sub["x"], sub["y"], s=80, alpha=0.85,
                       edgecolors="black", linewidths=0.5,
                       c=[color_map[g]] * len(sub))
            for _, r in sub.iterrows():
                # Remove "Sample_" prefix for cleaner labels
                label = r["col"].replace("Sample_", "").replace("sample_", "")
                plt.text(r["x"], r["y"] + text_offset, label, 
                        fontsize=6, ha="center", va="bottom")
        
        plt.xlabel("t-SNE 1", fontsize=10)
        plt.ylabel("t-SNE 2", fontsize=10)
        plt.title(f"t-SNE: {g1} vs {g2}", fontsize=12, fontweight='bold')
        # Legend removed
        plt.tight_layout()
        
        # Safe filename
        safe_g1 = g1.replace(" ", "_").replace("/", "_")
        safe_g2 = g2.replace(" ", "_").replace("/", "_")
        out_file = f"{out_prefix}_compare_{safe_g1}_vs_{safe_g2}.png"
        plt.savefig(out_file, bbox_inches="tight")
        plt.close()
        print(f"  Saved {out_file}")



def main():
    args = parse_args()
    in_path = Path(args.csv)
    if not in_path.exists():
        sys.exit(f"Input file not found: {in_path}")

    df = pd.read_csv(in_path)
    df = choose_index_column(df, args.index_col)

    # Identify sample columns and parse headers
    sample_cols = []
    sample_meta = []  # list of dicts with name, id, group
    for c in df.columns:
        parsed = parse_header(c)
        if parsed is not None:
            sid, grp = parsed
            sample_cols.append(c)
            sample_meta.append({"col": c, "id": str(sid), "group": str(grp)})

    if len(sample_cols) == 0:
        sys.exit("No sample columns matched the expected header pattern. Example: 'Sample_556_Control'.")

    meta_df = pd.DataFrame(sample_meta)
    groups = {g: meta_df.loc[meta_df["group"] == g, "col"].tolist()
              for g in sorted(meta_df["group"].unique())}

    print("Detected groups and sample counts:")
    for g, cols in groups.items():
        print(f"  {g}: {len(cols)}")

    # Extract sample data matrix (features x samples)
    X = df[sample_cols].copy()

    # Convert to numeric and handle missing
    X = X.apply(pd.to_numeric, errors="coerce")
    
    # Basic imputation by feature median
    # For features where median is NaN (all values missing), use 0
    def impute_row(row):
        med = row.median()
        if pd.isna(med):
            return row.fillna(0)
        return row.fillna(med)
    
    X = X.apply(impute_row, axis=1)
    
    # Remove features that are still all NaN or all zeros (uninformative)
    row_valid = X.notna().any(axis=1) & (X != 0).any(axis=1)
    if not row_valid.all():
        n_removed = (~row_valid).sum()
        print(f"Removing {n_removed} uninformative features (all NaN or all zeros)")
        X = X[row_valid]
    
    # Final safety check: ensure no NaN values remain
    if X.isna().any().any():
        print("Warning: NaN values still present after imputation, filling with 0")
        X = X.fillna(0)

    # Optional Q3 normalization
    if args.q3:
        X = q3_normalize(X)

    # One-way ANOVA to rank features
    f_series = anova_f_vector(X, groups)
    X_sel, f_sel = pick_top_features_by_f(X, f_series, args.top_percent)

    # Prepare t-SNE input (samples x features)
    X_tsne = X_sel.T.values
    n_samples = X_tsne.shape[0]

    # Adjust perplexity if needed
    max_perp = max(5.0, (n_samples - 1) / 3.0)
    perp = min(args.perplexity, max_perp)
    if perp < 5.0:
        perp = 5.0

    print(f"t-SNE with perplexity={perp:.1f}, n_iter={args.n_iter}, samples={n_samples}, features={X_sel.shape[0]}")

    tsne = TSNE(n_components=2, perplexity=perp, n_iter=args.n_iter, random_state=args.seed, init="pca")
    coords = tsne.fit_transform(X_tsne)

    # Build a clean plotting frame
    plot_df = meta_df.copy()
    plot_df["x"] = coords[:, 0]
    plot_df["y"] = coords[:, 1]

    # Colors only - using circles for all groups
    color_map = make_color_map(plot_df["group"].unique())

    # Plot
    plt.figure(figsize=(10, 8), dpi=150)
    
    # Calculate adaptive text offset based on y-axis range
    y_range = plot_df["y"].max() - plot_df["y"].min()
    text_offset = y_range * 0.02  # 2% of y-range
    
    for g, sub in plot_df.groupby("group"):
        plt.scatter(sub["x"], sub["y"], s=80, alpha=0.85, 
                    edgecolors="black", linewidths=0.5,
                    c=[color_map[g]] * len(sub))
        # Add sample labels slightly above each point
        for _, r in sub.iterrows():
            # Remove "Sample_" prefix for cleaner labels
            label = r["col"].replace("Sample_", "").replace("sample_", "")
            plt.text(r["x"], r["y"] + text_offset, label, fontsize=6, ha="center", va="bottom")

    plt.xlabel("t-SNE 1", fontsize=11)
    plt.ylabel("t-SNE 2", fontsize=11)
    plt.title("t-SNE Analysis by Group", fontsize=13, fontweight='bold', pad=15)
    # Legend removed as per user request
    plt.tight_layout()

    out_png = f"{args.out_prefix}_plot.png"
    out_csv = f"{args.out_prefix}_coords.csv"
    plt.savefig(out_png, bbox_inches="tight")
    plt.close()

    print(f"Saved main plot to {out_png}")

    # Create separate comparison plots if requested
    if args.separate_plots:
        create_comparison_plots(plot_df, color_map, text_offset, args.out_prefix)

    # Save coordinates
    plot_df[["col", "id", "group", "x", "y"]].to_csv(out_csv, index=False)

    # Save selected features if feature selection applied
    if args.top_percent < 100.0:
        out_feats = f"{args.out_prefix}_features.txt"
        X_sel.index.to_series().to_csv(out_feats, index=False, header=False)
        print(f"Saved selected features to {out_feats}")

    print(f"Saved coordinates to {out_csv}")
    print("Done!")


if __name__ == "__main__":
    main()
