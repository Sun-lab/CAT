#!/usr/bin/env python
"""
Clonality of cells filtered by IEDB / McPAS (both-chain non-human-antigen match),
stratified by whether the TCR is shared with VDJdb (score>=2) or database-only.

Combines IEDB + McPAS into one table per cell type (CD8, CD4) and writes both
tables to a Word document.

Reuses the database loaders / normalization from zz3_db_compare.py.
Output: ~/research/CAT/_revision/db_comparison/clonality_filtered_cells.docx
"""
import os
import numpy as np
import pandas as pd
from docx import Document
from docx.shared import Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH

import zz3_db_compare as z

OUTDIR = os.path.expanduser("~/research/CAT/_revision/db_comparison")
DOCX   = os.path.join(OUTDIR, "clonality_filtered_cells.docx")


def both_cells(cr, alpha, beta):
    """Per-cell boolean: both chains match the database (exact CDR3+V)."""
    ak = set(alpha["cdr3"] + "|" + alpha["v_gene"])
    bk = set(beta["cdr3"] + "|" + beta["v_gene"])
    a = (cr["a_cdr3"] + "|" + cr["a_v"]).isin(ak).to_numpy()
    b = (cr["b_cdr3"] + "|" + cr["b_v"]).isin(bk).to_numpy()
    return a & b


def clonality_table(cell_type, dbs):
    """Combined IEDB + McPAS clonality table for one cell type."""
    cr = z.load_cr(cell_type)
    cr["_clone"] = list(zip(cr.TRA_cdr3, cr.TRA_v_gene, cr.TRB_cdr3, cr.TRB_v_gene))
    size_map = cr["_clone"].value_counts().to_dict()
    vdj_clones = set(cr.loc[both_cells(cr, *dbs["VDJdb"]), "_clone"])

    rows = []
    for dbname in ("IEDB", "McPAS"):
        db_clones = set(cr.loc[both_cells(cr, *dbs[dbname]), "_clone"])
        groups = [("shared with VDJdb", db_clones & vdj_clones),
                  ("DB-only (not in VDJdb)", db_clones - vdj_clones),
                  ("all", db_clones)]
        for label, clones in groups:
            sz = np.array([size_map[c] for c in clones]) if clones else np.array([0])
            rows.append({
                "database": dbname, "group": label,
                "n_TCRs": len(clones), "n_cells": int(sz.sum()),
                "mean_cells_per_TCR": round(float(sz.mean()), 1),
                "median": int(np.median(sz)), "max": int(sz.max()),
                "pct_singleton_TCRs": round(100 * float(np.mean(sz == 1)), 1),
            })
    return pd.DataFrame(rows)


def add_df_table(doc, df, heading):
    doc.add_heading(heading, level=2)
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Light Grid Accent 1"
    for j, col in enumerate(df.columns):
        cell = t.rows[0].cells[j]
        cell.text = col
        for p in cell.paragraphs:
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            for r in p.runs:
                r.bold = True; r.font.size = Pt(9)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for j, col in enumerate(df.columns):
            v = row[col]
            cells[j].text = f"{v:,}" if isinstance(v, (int, np.integer)) else str(v)
            for p in cells[j].paragraphs:
                if j >= 2:
                    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                for r in p.runs:
                    r.font.size = Pt(9)


def main():
    dbs = {"VDJdb": z.db_vdjdb(), "McPAS": z.db_mcpas(), "IEDB": z.db_iedb()}
    tables = {ct: clonality_table(ct, dbs) for ct in ("CD8", "CD4")}

    doc = Document()
    doc.styles["Normal"].font.name = "Calibri"
    doc.styles["Normal"].font.size = Pt(10)
    doc.add_heading("Clonality of cells filtered by IEDB / McPAS, "
                    "stratified by sharing with VDJdb", level=1)
    doc.add_paragraph(
        "Cells are filtered when both TCR chains match a non-human-antigen TCR "
        "(exact CDR3 + V gene). For each database the filtered clonotypes are split "
        "into those also flagged by VDJdb (score ≥ 2) vs database-only, and "
        "summarized by clone size (cells per TCR).")
    for ct in ("CD8", "CD4"):
        tables[ct].to_csv(os.path.join(OUTDIR, f"clonality_filtered_{ct}.csv"), index=False)
        add_df_table(doc, tables[ct], f"{ct}")
        print(f"\n=== {ct} ===")
        print(tables[ct].to_string(index=False))
    doc.save(DOCX)
    print(f"\n[written] {DOCX}")


if __name__ == "__main__":
    main()
