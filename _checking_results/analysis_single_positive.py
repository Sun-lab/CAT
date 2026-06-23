#!/usr/bin/env python
"""
Summarize the two-layer cancer-reactive definition for a single study.

For a given study folder it reports, separately for CD4 and CD8:
  1. The per-cluster reactive fraction (the quantity used to call a cluster
     "cancer-reactive": fraction of per-cell-positive cells in the cluster,
     threshold > 0.4), plus the resulting cluster-level call.
  2. The 2x2 cross-tabulation of the two component flags
        cancer_reactive_per_cell  (the "overall" signature call)
        cancer_reactive_by_cluster
     so the "single-positive" cells (high overall but not in a reactive
     cluster, and in a reactive cluster but not high overall) are counted
     explicitly, alongside the final AND-based call.

Input: the per-cell tables written by process_output.ipynb, i.e.
       <study>/cell_meta_data_CD8.csv and <study>/cell_meta_data_CD4.csv
       (these carry the columns cluster, cancer_reactive_per_cell,
        cancer_reactive_by_cluster).

Usage:
    python analysis_single_positive.py Zheng_2021
    python analysis_single_positive.py Zheng_2021 --write-table
    python analysis_single_positive.py Chen_2024 --threshold 0.4 --write-table
"""
import argparse
import os
import sys
import pandas as pd

CLUSTER_THRESHOLD = 0.4  # a cluster is reactive if reactive_fraction > this


def load_cell_meta(path):
    """Read a cell_meta_data_*.csv regardless of single- or multi-index header."""
    df = pd.read_csv(path, low_memory=False)
    # coerce the two flags to real booleans (csv may store as True/False strings)
    for col in ["cancer_reactive_per_cell", "cancer_reactive_by_cluster"]:
        if col not in df.columns:
            raise KeyError(
                f"{path} is missing column '{col}'. "
                "Re-run process_output.ipynb to regenerate it."
            )
        df[col] = df[col].astype(str).str.strip().str.lower().map(
            {"true": True, "false": False, "1": True, "0": False}
        )
    return df


def resolve_cluster_col(df, celltype):
    """Cluster label column is 'cluster' in most studies, 'cluster_<CT>' in others."""
    for cand in ("cluster", f"cluster_{celltype}"):
        if cand in df.columns:
            return cand
    raise KeyError(
        f"No cluster column found (looked for 'cluster' / 'cluster_{celltype}'). "
        f"Available: {list(df.columns)}"
    )


def per_cluster_fractions(df, threshold, cluster_col="cluster"):
    """Reactive fraction per cluster + the cluster-level reactive call.

    cluster_is_reactive is taken from the STORED cancer_reactive_by_cluster flag
    (the authoritative call made on the full data in process_output.ipynb), so it
    stays consistent with the published selection even when this file is a TCR-
    filtered subset. reactive_fraction is recomputed on the available cells for
    description; reactive_fraction_gt_threshold flags where the two disagree
    (borderline clusters near the threshold).
    """
    g = (
        df.groupby(cluster_col)
        .agg(
            n_cells=("cancer_reactive_per_cell", "size"),
            n_pos=("cancer_reactive_per_cell", "sum"),
            cluster_is_reactive=("cancer_reactive_by_cluster", "first"),
        )
    )
    g["reactive_fraction"] = g["n_pos"] / g["n_cells"]
    g["reactive_fraction_gt_threshold"] = g["reactive_fraction"] > threshold
    g.index.name = "cluster_id"
    g = g[["n_cells", "n_pos", "reactive_fraction",
           "cluster_is_reactive", "reactive_fraction_gt_threshold"]]
    return g.sort_values("reactive_fraction", ascending=False)


def single_positive_crosstab(df):
    """2x2 of per-cell vs by-cluster flags, with the AND result summarized."""
    ct = pd.crosstab(
        df["cancer_reactive_per_cell"].rename("overall_high"),
        df["cancer_reactive_by_cluster"].rename("in_reactive_cluster"),
    )
    # named cells of the 2x2
    def cell(a, b):
        try:
            return int(ct.loc[a, b])
        except KeyError:
            return 0

    summary = {
        "both (final cancer_reactive=True)": cell(True, True),
        "single-positive: overall-high only (per_cell & ~cluster)": cell(True, False),
        "single-positive: in-cluster only (cluster & ~per_cell)": cell(False, True),
        "neither": cell(False, False),
        "total cells": int(ct.to_numpy().sum()),
    }
    return ct, summary


def write_supp_tables(study, celltype, frac, summary, outdir):
    """Write a publication-ready supplementary table (one CSV per cell type)."""
    os.makedirs(outdir, exist_ok=True)

    # Sheet 1: per-cluster reactive fraction + cluster call
    frac_out = frac.reset_index()
    frac_path = os.path.join(outdir, f"{study}_{celltype}_cluster_fractions.csv")
    frac_out.to_csv(frac_path, index=False)

    # Sheet 2: single-positive cross-tab summary (long format, tidy for a table)
    sp = pd.DataFrame(
        [{"category": k, "n_cells": v} for k, v in summary.items()
         if k != "total cells"]
    )
    tot = summary["total cells"]
    sp["pct_cells"] = (100 * sp["n_cells"] / tot).round(2) if tot else 0.0
    sp.insert(0, "cell_type", celltype)
    sp.insert(0, "study", study)
    sp_path = os.path.join(outdir, f"{study}_{celltype}_single_positive_summary.csv")
    sp.to_csv(sp_path, index=False)

    print(f"    [written] {frac_path}")
    print(f"    [written] {sp_path}")


def report_celltype(study, celltype, threshold, outdir=None, write_table=False,
                    datadir=None, suffix="_cleaned"):
    label = os.path.basename(os.path.normpath(study))  # for filenames/titles
    base = os.path.join(datadir, label) if datadir else study
    path = os.path.join(base, f"cell_meta_data_{celltype}{suffix}.csv")
    if not os.path.exists(path):
        print(f"[skip] {path} not found")
        return
    print("=" * 70)
    print(f"{label}  —  {celltype}")
    print("=" * 70)
    df = load_cell_meta(path)
    cluster_col = resolve_cluster_col(df, celltype)

    print(f"\n[1] Per-cluster reactive fraction (cluster reactive if > {threshold}):")
    frac = per_cluster_fractions(df, threshold, cluster_col)
    with pd.option_context("display.float_format", "{:.3f}".format):
        print(frac.to_string())
    n_react_clusters = int(frac["cluster_is_reactive"].sum())
    print(f"\n    -> {n_react_clusters}/{len(frac)} clusters called cancer-reactive")

    print("\n[2] overall-high (per_cell) x in-reactive-cluster cross-tab:")
    ct, summary = single_positive_crosstab(df)
    print(ct.to_string())
    print("\n    Counts:")
    for k, v in summary.items():
        print(f"      {k}: {v}")
    tot = summary["total cells"]
    sp1 = summary["single-positive: overall-high only (per_cell & ~cluster)"]
    sp2 = summary["single-positive: in-cluster only (cluster & ~per_cell)"]
    if tot:
        print(
            f"\n    Single-positive cells = {sp1 + sp2} "
            f"({100 * (sp1 + sp2) / tot:.1f}% of all cells) "
            "are excluded by the AND criterion."
        )

    if write_table:
        print()
        write_supp_tables(label, celltype, frac, summary, outdir)
    print()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("study", help="study folder, e.g. Zheng_2021")
    ap.add_argument("--threshold", type=float, default=CLUSTER_THRESHOLD,
                    help="cluster reactive-fraction threshold (default 0.4)")
    ap.add_argument("--celltypes", nargs="+", default=["CD8", "CD4"])
    ap.add_argument("--write-table", action="store_true",
                    help="write supplementary CSV tables")
    ap.add_argument("--outdir", default=None,
                    help="output dir for tables "
                         "(default <study>/single_positive_analysis)")
    ap.add_argument("--datadir", default=None,
                    help="base dir holding the study subfolders "
                         "(default: the study path itself)")
    ap.add_argument("--suffix", default="_cleaned",
                    help="cell_meta_data file suffix (default '_cleaned')")
    args = ap.parse_args()

    outdir = args.outdir or os.path.join(args.study, "single_positive_analysis")

    for ct in args.celltypes:
        report_celltype(args.study, ct, args.threshold, outdir=outdir,
                        write_table=args.write_table,
                        datadir=args.datadir, suffix=args.suffix)


if __name__ == "__main__":
    main()
