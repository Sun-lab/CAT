#!/usr/bin/env python
"""
Reviewer comment 3.3 — would adding McPAS-TCR and/or IEDB to the VDJdb-based
"non-human antigen" removal change the cancer-reactive (CR) TCR set?

Mirrors the VDJdb filter in zz2_construct_training_data.Rmd WITHOUT modifying the
pipeline. For each cell with a paired alpha/beta TCR, each chain is matched
against a database of human TCRs recognizing a non-human (microbial/viral)
antigen; a cell is flagged for removal only if BOTH chains match (the zz2 rule
non_human_antigen == "Both"). We then report, per database, how many additional
cancer-reactive cells would be removed beyond VDJdb.

To compare databases fairly despite different CDR3 / V-gene conventions, all
sources (CR cells + databases) are normalized identically:
  - V gene: uppercase, TCR->TR, drop allele (*..) and any ", ..." list tail
  - CDR3 : uppercase, drop a leading 'C' and trailing 'F'/'W' (conserved
           junction residues) so junction vs core encodings align
Two match modes: exact (CDR3 core + V gene) and cdr3-only (lenient upper bound).
"""
import os
import re
import pandas as pd

HOME   = os.path.expanduser("~")
CR_DIR = f"{HOME}/research/CAT/data"
VDJDB  = f"{HOME}/research/GitHub/cancer_reactive_TCR/data_VDJdb/vdjdb-2025-02-21/vdjdb.slim.txt"
MCPAS  = f"{HOME}/research/GitHub/cancer_reactive_TCR/data_McPAS/McPAS-TCR.csv.gz"
IEDB   = f"{HOME}/research/GitHub/cancer_reactive_TCR/data_IEDB/tcr_full_v3.zip"
OUTDIR = f"{HOME}/research/CAT/_revision/db_comparison"
os.makedirs(OUTDIR, exist_ok=True)

CLONO = ["TRA_cdr3", "TRA_v_gene", "TRB_cdr3", "TRB_v_gene"]


def norm_v(s):
    if pd.isna(s): return ""
    s = str(s).replace("\xa0", "").strip().upper().replace("TCR", "TR")
    s = s.split(",")[0]            # first gene of a list
    s = s.split("*")[0]           # drop allele
    return s.strip()


def norm_cdr3(s):
    if pd.isna(s): return ""
    s = str(s).replace("\xa0", "").strip().upper()
    if not s or not s.isalpha(): return s if s else ""
    if s.startswith("C"): s = s[1:]
    if s.endswith(("F", "W")): s = s[:-1]
    return s


def prep_db(alpha, beta):
    """alpha/beta: DataFrames with raw columns ['cdr3','v_gene']; normalize."""
    out = []
    for d in (alpha, beta):
        d = d.copy()
        d["cdr3"] = d["cdr3"].map(norm_cdr3)
        d["v_gene"] = d["v_gene"].map(norm_v)
        d = d[(d["cdr3"] != "")]
        out.append(d.drop_duplicates())
    return out


def load_cr(cell_type):
    df = pd.read_csv(f"{CR_DIR}/{cell_type}_combined.tsv.gz", sep="\t", low_memory=False)
    for c in CLONO:
        df[c] = df[c].fillna("").astype(str)
    df = df[(df["TRA_cdr3"] != "") & (df["TRB_cdr3"] != "")].copy()
    df["cancer_reactive"] = (df["cancer_reactive_per_cell"].astype(bool)
                             & df["study_specific_CR_by_cluster"].astype(bool))
    df["a_cdr3"] = df["TRA_cdr3"].map(norm_cdr3); df["a_v"] = df["TRA_v_gene"].map(norm_v)
    df["b_cdr3"] = df["TRB_cdr3"].map(norm_cdr3); df["b_v"] = df["TRB_v_gene"].map(norm_v)
    return df


def db_vdjdb():
    v = pd.read_csv(VDJDB, sep="\t", low_memory=False)
    v = v[(v["species"] == "HomoSapiens") & (v["vdjdb.score"] >= 2)].copy()
    v["v_gene"] = v["v.segm"]
    v = v[v["antigen.species"] != "HomoSapiens"]
    a = v[v["gene"] == "TRA"][["cdr3", "v_gene"]]
    b = v[v["gene"] == "TRB"][["cdr3", "v_gene"]]
    return prep_db(a, b)


def db_mcpas():
    m = pd.read_csv(MCPAS, encoding="latin-1", low_memory=False)
    m = m[(m["Species"] == "Human") & (m["Category"] == "Pathogens")]
    a = m[["CDR3.alpha.aa", "TRAV"]].rename(columns={"CDR3.alpha.aa": "cdr3", "TRAV": "v_gene"})
    b = m[["CDR3.beta.aa", "TRBV"]].rename(columns={"CDR3.beta.aa": "cdr3", "TRBV": "v_gene"})
    return prep_db(a, b)


def db_iedb():
    # IEDB constant points to the persisted zip containing tcr_full_v3.csv
    read_arg = IEDB
    if IEDB.endswith(".zip"):
        import zipfile, io
        with zipfile.ZipFile(IEDB) as z:
            name = [n for n in z.namelist() if n.endswith(".csv")][0]
            read_arg = io.BytesIO(z.read(name))
    df = pd.read_csv(read_arg, header=[0, 1], low_memory=False)
    df.columns = [" / ".join([str(x), str(y)]) for x, y in df.columns]
    C = df.columns
    src = C[8]
    # non-human antigen = source organism present and not Homo sapiens
    nonhuman = df[src].notna() & (~df[src].astype(str).str.startswith("Homo sapiens"))
    d = df[nonhuman]
    # chain1 = alpha (TRA), chain2 = beta (TRB); prefer curated then calculated
    def pick(row, cur, calc):
        v = row[cur]
        return v if (pd.notna(v) and str(v).strip() != "") else row[calc]
    a_cdr3 = d.apply(lambda r: pick(r, C[23], C[24]), axis=1)
    a_v    = d.apply(lambda r: pick(r, C[15], C[16]), axis=1)
    b_cdr3 = d.apply(lambda r: pick(r, C[54], C[55]), axis=1)
    b_v    = d.apply(lambda r: pick(r, C[46], C[47]), axis=1)
    a = pd.DataFrame({"cdr3": a_cdr3, "v_gene": a_v})
    b = pd.DataFrame({"cdr3": b_cdr3, "v_gene": b_v})
    return prep_db(a, b)


def chain_hits(df, alpha, beta, mode):
    """Return (a_hit, b_hit) boolean arrays: does each chain match the database."""
    if mode == "exact":
        ak = set(zip(alpha["cdr3"], alpha["v_gene"])); bk = set(zip(beta["cdr3"], beta["v_gene"]))
        a_hit = df.apply(lambda r: (r["a_cdr3"], r["a_v"]) in ak, axis=1)
        b_hit = df.apply(lambda r: (r["b_cdr3"], r["b_v"]) in bk, axis=1)
    else:
        ak = set(alpha["cdr3"]); bk = set(beta["cdr3"])
        a_hit = df["a_cdr3"].isin(ak); b_hit = df["b_cdr3"].isin(bk)
    return a_hit.to_numpy(), b_hit.to_numpy()


def flag_both(df, alpha, beta, mode):
    a_hit, b_hit = chain_hits(df, alpha, beta, mode)
    return a_hit & b_hit


def nonhuman_category(a_hit, b_hit):
    """Reproduce zz2 add_non_human_antigen(): Both / TRA / TRB / None."""
    import numpy as np
    cat = np.full(len(a_hit), "None", dtype=object)
    cat[a_hit & b_hit]  = "Both"
    cat[a_hit & ~b_hit] = "TRA"
    cat[~a_hit & b_hit] = "TRB"
    return cat


def make_iedb_plot(dbs):
    """Reproduce the zz2 'draw_bar_plot' panels using IEDB annotations:
    proportion of cancer-reactive cells in each non-human-antigen category.
    Four panels = {CD8, CD4} x {cancer_reactive_per_cell, cancer_reactive}.
    """
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    alpha, beta = dbs["IEDB"]
    order = ["Both", "TRB", "TRA", "None"]            # same order as zz2
    cr_cols = ["cancer_reactive_per_cell", "cancer_reactive"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    prop_rows = []

    for i, ct in enumerate(["CD8", "CD4"]):
        cr = load_cr(ct)
        cr["cancer_reactive_per_cell"] = cr["cancer_reactive_per_cell"].astype(bool)
        a_hit, b_hit = chain_hits(cr, alpha, beta, "exact")   # exact = pipeline-faithful
        cat = nonhuman_category(a_hit, b_hit)
        for j, col in enumerate(cr_cols):
            ax = axes[i, j]
            sub = pd.DataFrame({"cat": cat, "cr": cr[col].to_numpy().astype(bool)})
            grp = sub.groupby("cat")["cr"].agg(total="size", n_reactive="sum")
            grp["prop"] = grp["n_reactive"] / grp["total"]
            grp = grp.reindex(order)
            ax.bar(order, grp["prop"].values, color="steelblue")
            for x, (p, n) in enumerate(zip(grp["prop"].values, grp["total"].values)):
                if pd.notna(p):
                    ax.text(x, p + 0.01, f"n={int(n):,}", ha="center", va="bottom", fontsize=7)
            ax.set_ylim(0, 0.4)
            ax.set_xlabel("Antigen Category (IEDB, non-human)")
            ax.set_ylabel("Prop of Cancer Reactive")
            ax.set_title(f"{ct}: {col.replace('_',' ')}")
            for cat_name in order:
                r = grp.loc[cat_name]
                prop_rows.append({"cell_type": ct, "cr_definition": col,
                                  "antigen_category": cat_name,
                                  "n_total": int(r["total"]) if pd.notna(r["total"]) else 0,
                                  "n_cancer_reactive": int(r["n_reactive"]) if pd.notna(r["n_reactive"]) else 0,
                                  "proportion_reactive": round(float(r["prop"]), 4) if pd.notna(r["prop"]) else None})

    fig.suptitle("Proportion of cancer-reactive T cells by IEDB non-human-antigen category",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out_png = f"{OUTDIR}/iedb_prop_cancer_reactive_by_antigen.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    pd.DataFrame(prop_rows).to_csv(
        f"{OUTDIR}/iedb_prop_cancer_reactive_by_antigen.csv", index=False)
    print(f"[written] {out_png}")
    print(f"[written] {OUTDIR}/iedb_prop_cancer_reactive_by_antigen.csv")


def run(cell_type, dbs):
    print("=" * 74); print(cell_type); print("=" * 74)
    cr = load_cr(cell_type)
    crpos = cr["cancer_reactive"].to_numpy(); n_cr = int(crpos.sum())
    print(f"CR set: {len(cr)} paired-TCR cells, {n_cr} cancer-reactive")

    rows = []
    for mode in ("exact", "cdr3"):
        flags = {name: flag_both(cr, a, b, mode) for name, (a, b) in dbs.items()}
        vdj = flags["VDJdb"]
        print(f"\n--- {mode} match ---")
        for name, fl in flags.items():
            removed_cr = int((fl & crpos).sum())
            add_cr = int((fl & crpos & ~vdj).sum()) if name != "VDJdb" else 0
            tag = "" if name == "VDJdb" else f"; ADDITIONAL beyond VDJdb: {add_cr} ({100*add_cr/n_cr:.3f}%)"
            print(f"  {name:7s}: removes {removed_cr} CR cells ({100*removed_cr/n_cr:.3f}%){tag}")
            rows.append({"cell_type": cell_type, "match_mode": mode, "database": name,
                         "n_cancer_reactive": n_cr, "CR_cells_removed": removed_cr,
                         "additional_vs_VDJdb": add_cr})
        # union of all three
        union = vdj.copy()
        for fl in flags.values(): union = union | fl
        add_union = int((union & crpos & ~vdj).sum())
        print(f"  UNION(McPAS+IEDB) adds {add_union} CR cells ({100*add_union/n_cr:.3f}%) beyond VDJdb")
        rows.append({"cell_type": cell_type, "match_mode": mode, "database": "McPAS+IEDB union",
                     "n_cancer_reactive": n_cr, "CR_cells_removed": int((union & crpos).sum()),
                     "additional_vs_VDJdb": add_union})
    return rows


if __name__ == "__main__":
    print("Loading databases (non-human-antigen human TCRs)...")
    dbs = {"VDJdb": db_vdjdb(), "McPAS": db_mcpas(), "IEDB": db_iedb()}
    for n, (a, b) in dbs.items():
        print(f"  {n}: {len(a)} TRA, {len(b)} TRB normalized entries")
    all_rows = []
    for ct in ("CD8", "CD4"):
        all_rows += run(ct, dbs)
    pd.DataFrame(all_rows).to_csv(f"{OUTDIR}/db_comparison_summary.csv", index=False)
    print("\n[written]", f"{OUTDIR}/db_comparison_summary.csv")

    print("\nGenerating IEDB proportion-cancer-reactive plot (zz2-style)...")
    make_iedb_plot(dbs)
