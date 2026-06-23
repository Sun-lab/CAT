#!/usr/bin/env python
"""
Run the single-positive / cluster-fraction analysis across the production
studies and combine the results into ONE file each.

Reuses the logic in analysis_single_positive.py. For every study x cell type
that has a <study>/cell_meta_data_<CT>.csv, it appends:
  - the per-cluster reactive-fraction rows  -> combined_cluster_fractions.csv
  - the single-positive summary rows        -> combined_single_positive_summary.csv
Studies whose input file is missing are reported and skipped (not fabricated).

Usage:
    python combine_single_positive.py
    python combine_single_positive.py --studies Zheng_2021 Chen_2024 --threshold 0.4
"""
import argparse
import os
import sys
import pandas as pd

from analysis_single_positive import (
    CLUSTER_THRESHOLD,
    load_cell_meta,
    per_cluster_fractions,
    resolve_cluster_col,
    single_positive_crosstab,
)

DEFAULT_STUDIES = ["Zheng_2021", "Chen_2024", "Chow_2023", "Liu_2022", "Liu_2025"]
DEFAULT_CELLTYPES = ["CD8", "CD4"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--studies", nargs="+", default=DEFAULT_STUDIES)
    ap.add_argument("--celltypes", nargs="+", default=DEFAULT_CELLTYPES)
    ap.add_argument("--threshold", type=float, default=CLUSTER_THRESHOLD)
    ap.add_argument("--datadir", default=".",
                    help="base dir holding the study subfolders")
    ap.add_argument("--suffix", default="_cleaned",
                    help="cell_meta_data file suffix (default '_cleaned')")
    ap.add_argument("--outdir", default="combined_single_positive_analysis")
    args = ap.parse_args()

    frac_rows, summary_rows, missing = [], [], []

    for study in args.studies:
        label = os.path.basename(os.path.normpath(study))
        for ct in args.celltypes:
            path = os.path.join(args.datadir, label,
                                f"cell_meta_data_{ct}{args.suffix}.csv")
            if not os.path.exists(path):
                missing.append(path)
                continue
            df = load_cell_meta(path)
            cluster_col = resolve_cluster_col(df, ct)

            # per-cluster fractions
            frac = per_cluster_fractions(df, args.threshold, cluster_col).reset_index()
            frac.insert(0, "cell_type", ct)
            frac.insert(0, "study", label)
            frac_rows.append(frac)

            # single-positive summary
            _, summary = single_positive_crosstab(df)
            tot = summary["total cells"]
            sp = pd.DataFrame(
                [{"category": k, "n_cells": v} for k, v in summary.items()
                 if k != "total cells"]
            )
            sp["pct_cells"] = (100 * sp["n_cells"] / tot).round(2) if tot else 0.0
            sp.insert(0, "n_total_cells", tot)
            sp.insert(0, "cell_type", ct)
            sp.insert(0, "study", label)
            summary_rows.append(sp)

            print(f"[ok] {label} {ct}: {len(frac)} clusters, {tot} cells")

    if missing:
        print("\n[missing inputs - skipped] "
              "(run process_output.ipynb or copy these files in):")
        for m in missing:
            print(f"   - {m}")

    if not frac_rows:
        sys.exit("\nNo input tables found; nothing to combine.")

    os.makedirs(args.outdir, exist_ok=True)
    frac_all = pd.concat(frac_rows, ignore_index=True)
    summary_all = pd.concat(summary_rows, ignore_index=True)
    frac_path = os.path.join(args.outdir, "combined_cluster_fractions.csv")
    summary_path = os.path.join(args.outdir, "combined_single_positive_summary.csv")
    frac_all.to_csv(frac_path, index=False)
    summary_all.to_csv(summary_path, index=False)

    print(f"\n[written] {frac_path}  ({len(frac_all)} rows)")
    print(f"[written] {summary_path}  ({len(summary_all)} rows)")


if __name__ == "__main__":
    main()
