"""Functional category enrichment of the cancer-reactive (CR) gene signatures.

Run from this folder:  python enrichment.py
Outputs: signature_genes.tsv, enrichment_results.tsv, enrichment_dotplot.pdf

Test: one-sided hypergeometric test of each signature union (CD8/CD4 x
positive/negative) against MSigDB GO BP, Reactome and Hallmark gene sets.
Background: genes expressed in >= 5% of cells (the genes in DE_{CD8,CD4}.csv).
FDR: Benjamini-Hochberg within each gene set and library.
"""
import pickle, textwrap, urllib.request
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import hypergeom, false_discovery_control

here = Path(__file__).resolve().parent
data, de_dir = here.parent, here.parents[1] / "classify_cancer_reactive_T/DE_results"
msigdb = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/"
libs = {"Hallmark": "h.all", "Reactome": "c2.cp.reactome", "GO BP": "c5.go.bp"}
min_size, max_size, fdr_cut, n_top = 10, 500, 0.05, 5


def read_gmt(prefix):
    f = here / "gmt" / f"{prefix}.v2024.1.Hs.symbols.gmt"
    if not f.exists():
        f.parent.mkdir(exist_ok=True)
        urllib.request.urlretrieve(msigdb + f.name, f)
    return {x[0]: set(x[2:]) for x in (l.rstrip("\n").split("\t") for l in open(f))}


gmt = {lib: read_gmt(prefix) for lib, prefix in libs.items()}

# ---- 1. signature gene table, ordered by the number of signatures per gene ----
genes = []
for ct in ["CD8", "CD4"]:
    sig = pickle.load(open(data / f"signatures_{ct}.pkl", "rb"))
    used = pd.read_csv(data / f"scale_{ct}.tsv", sep="\t")["signature"]
    for name, gs in sig.items():
        if used.str.endswith(name).any():  # keep the signatures used for scoring
            direction = "negative" if ("neg" in name or "virus" in name) else "positive"
            genes += [(ct, direction, g, name) for g in gs]
genes = pd.DataFrame(genes, columns=["cell_type", "direction", "gene", "signature"])
tab = (genes.groupby(["cell_type", "direction", "gene"])["signature"]
       .agg(n_signatures="nunique", signatures=lambda s: ";".join(sorted(set(s))))
       .reset_index())
for ct in ["CD8", "CD4"]:  # flag genes in the enrichment background
    bg = set(pd.read_csv(de_dir / f"DE_{ct}.csv")["names"])
    tab.loc[tab.cell_type == ct, "in_background"] = tab.gene[tab.cell_type == ct].isin(bg)
tab = tab.sort_values(["cell_type", "direction", "n_signatures", "gene"],
                      ascending=[False, False, False, True])
tab.to_csv(here / "signature_genes.tsv", sep="\t", index=False)

# ---- 2. hypergeometric enrichment ----
res = []
for (ct, direction), d in tab[tab.in_background.astype(bool)].groupby(["cell_type", "direction"], sort=False):
    bg = set(pd.read_csv(de_dir / f"DE_{ct}.csv")["names"])
    query, N = set(d.gene), len(bg)
    for lib, terms in gmt.items():
        r = pd.DataFrame([(t, len(query & m), len(m & bg), ",".join(sorted(query & m)))
                          for t, m in terms.items()], columns=["term", "k", "K", "genes"])
        r = r[r.K.between(min_size, max_size)].copy()
        r["p"] = hypergeom.sf(r.k - 1, N, r.K, len(query))
        r["FDR"] = false_discovery_control(r.p)
        r["fold"] = (r.k / len(query)) / (r.K / N)
        res.append(r.assign(gene_set=f"{ct} {direction}", library=lib, n=len(query), N=N))
res = pd.concat(res).sort_values(["gene_set", "library", "p"])
cols = ["gene_set", "library", "term", "k", "K", "n", "N", "fold", "p", "FDR", "genes"]
sig_res = res.loc[res.FDR < fdr_cut, cols]
sig_res.to_csv(here / "enrichment_results.tsv", sep="\t", index=False, float_format="%.4g")

# ---- 3. dot plot: top terms per gene set and library ----
sets = ["CD8 positive", "CD8 negative", "CD4 positive", "CD4 negative"]
bins, labels = [0, 1e-6, 1e-3, fdr_cut], ["< 1e-6", "1e-6 to 1e-3", "1e-3 to 0.05"]
colors = dict(zip(labels, ["#08306b", "#4292c6", "#9ecae1"]))
size = lambda fold: 20 * np.minimum(fold, 10)  # fold enrichment truncated at 10
top = {lib: list(dict.fromkeys(t for s in sets for t in d[d.gene_set == s].head(n_top).term))
       for lib, d in sig_res.groupby("library")}
fig = plt.figure(figsize=(11, 6))
gs = fig.add_gridspec(2, 2, height_ratios=[len(top["Hallmark"]), len(top["Reactome"])],
                      left=0.215, right=0.86, top=0.95, bottom=0.08, wspace=0.95, hspace=0.14)
axes = {"Hallmark": fig.add_subplot(gs[0, 0]), "Reactome": fig.add_subplot(gs[1, 0]),
        "GO BP": fig.add_subplot(gs[:, 1])}
for lib, ax in axes.items():
    terms = top[lib]
    d = sig_res[(sig_res.library == lib) & sig_res.term.isin(terms)]
    fdr_bin = pd.cut(d.FDR, bins, labels=labels).astype(str)
    ax.scatter(d.gene_set.map(sets.index), d.term.map(terms.index), s=size(d.fold),
               c=fdr_bin.map(colors), edgecolors="white", linewidths=0.5, zorder=3)
    names = [textwrap.fill(t.split("_", 1)[1].replace("_", " "), 40) for t in terms]
    ax.set_yticks(range(len(terms)), names, fontsize=6, linespacing=1.0)
    ax.set_xticks(range(len(sets)), [s.replace(" ", "\n") for s in sets], fontsize=8)
    ax.set_ylim(len(terms) - 0.5, -0.5)
    ax.set_xlim(-0.5, len(sets) - 0.5)
    ax.set_title(lib, fontsize=9, loc="left")
    ax.grid(color="0.92", linewidth=0.5, zorder=0)
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_color("0.8")
axes["Hallmark"].tick_params(labelbottom=False)
handles = [plt.scatter([], [], s=60, color=c, label=l) for l, c in colors.items()]
fig.legend(handles=handles, title="FDR", loc="upper left", bbox_to_anchor=(0.87, 0.95),
           frameon=False, fontsize=7, title_fontsize=8)
handles = [plt.scatter([], [], s=size(f), color="0.6", label=l) for f, l in [(2, "2"), (5, "5"), (10, "≥ 10")]]
fig.legend(handles=handles, title="Fold enrichment", loc="upper left", bbox_to_anchor=(0.87, 0.72),
           frameon=False, fontsize=7, title_fontsize=8, labelspacing=1.2)
fig.savefig(here / "enrichment_dotplot.pdf")
