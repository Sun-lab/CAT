"""
Apply the trained CD8 careTCR autoencoder+classifier to an external test
dataset (e.g. Caushi_2021) and produce per-cell predictions.

Matches the preprocessing chain in step1a_prepare_CD8.ipynb +
step4_classify_careT_CD8.ipynb exactly:
  (1) restrict test genes to the training common-gene space (genes2use),
      inserting zero columns for any common gene missing from the test data
  (2) compute per-cell total_counts on that common-gene matrix
      (this is the same definition step1a uses via sc.pp.calculate_qc_metrics)
  (3) library-size normalize to target_sum = 3000, then log1p
  (4) subset to the top-1000 DE genes used as model input, in the same order
  (5) load model/CD8.keras and predict
"""

import os
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import issparse, csr_matrix
from sklearn.metrics import roc_curve, auc
from keras.models import load_model
import matplotlib.pyplot as plt

from scipy.sparse import hstack

from collections import Counter

# -----------------------------------------------------------------------------
# Paths -- specify the cell type to use, either CD4 or CD8
# -----------------------------------------------------------------------------
cell_type = "CD8"

# -----------------------------------------------------------------------------
# Paths -- edit these for your environment
# -----------------------------------------------------------------------------
COMMON_GENES_FILE = "../data/common_genes.txt"
GENE_1000_FILE    = f"./model/{cell_type}_top1000_genes.txt"
MODEL_FILE        = f"./model/{cell_type}.keras"                # from step4
TEST_DATA_DIR     = "../data/Caushi_2021/for_python"
OUT_DIR           = "predictions/Caushi_2021_CD8"

TARGET_SUM        = 3000        # must match step4

os.makedirs(OUT_DIR, exist_ok=True)


# -----------------------------------------------------------------------------
# 1. Load the training common-gene list and the top-1000 DE gene list
# -----------------------------------------------------------------------------
genes2use = pd.read_csv(COMMON_GENES_FILE, header=None)[0].tolist()
print(f"Number of common genes (training space): {len(genes2use)}")

# load the genes in the GENE_1000_FILE
top_genes = pd.read_csv(GENE_1000_FILE, header=None)[0].tolist()
len(top_genes)

# Sanity check: every top gene should live in the common-gene space
missing_from_common = [g for g in top_genes if g not in set(genes2use)]
if missing_from_common:
    print(f"WARNING: {len(missing_from_common)} top genes are not in common_genes.txt "
          f"-- check that you're using the right files. Examples: {missing_from_common[:5]}")

# -----------------------------------------------------------------------------
# 2. Load the test data (Caushi_2021 example)
# -----------------------------------------------------------------------------
counts  = mmread(f"{TEST_DATA_DIR}/counts.mtx.gz").T.tocsr()   # cells x genes
genes_t = pd.read_csv(f"{TEST_DATA_DIR}/genes.tsv.gz",
                      header=None, sep="\t")[0].values
barcodes = pd.read_csv(f"{TEST_DATA_DIR}/barcodes.tsv.gz",
                       header=None, sep="\t")[0].values
metadata = pd.read_csv(f"{TEST_DATA_DIR}/metadata.csv.gz",
                       index_col=0, low_memory=False)

assert all(metadata.index == barcodes), "Metadata index does not match barcodes order!"

adata = anndata.AnnData(
    X=counts,
    obs=metadata,
    var=pd.DataFrame(index=genes_t),
)
adata.var["gene_name"] = adata.var.index.tolist()
adata.obs["barcode"] = adata.obs.index.tolist()
print(f"Raw test AnnData: {adata.shape}")


# -----------------------------------------------------------------------------
# 3. Reindex to the training common-gene space (genes2use)
#
#    Important: any common gene missing from the test data is inserted as a
#    zero column so the gene order and dimensionality match training exactly.
# -----------------------------------------------------------------------------
existing_genes = set(adata.var_names)
present  = [g for g in genes2use if g in existing_genes]
missing  = [g for g in genes2use if g not in existing_genes]
print(f"Common genes present in test:  {len(present)}")
print(f"Common genes missing from test: {len(missing)}  (inserted as zeros)")

# see how many missing genes are among the top-1000 DE genes used as model input
missing_top = [g for g in top_genes if g in missing]
print(f"Missing common genes among top-1000 DE: {len(missing_top)} (these will be zeroed out in the model input)")

# Subset to present genes (preserve genes2use order among the present ones)
adata_pres = adata[:, present].copy()

# Pad with zero columns for missing common genes, then reorder to genes2use
if missing:
    n_cells = adata_pres.n_obs
    zero_block = csr_matrix((n_cells, len(missing)), dtype=adata_pres.X.dtype)
    X_full = hstack([adata_pres.X, zero_block]).tocsr()
    var_full = pd.DataFrame(index=list(adata_pres.var_names) + missing)
    adata_common = anndata.AnnData(X=X_full, obs=adata_pres.obs.copy(), var=var_full)
    # Reorder columns into genes2use order
    adata_common = adata_common[:, genes2use].copy()
else:
    adata_common = adata_pres[:, genes2use].copy()

assert list(adata_common.var_names) == list(genes2use), "Gene-order mismatch after reindex!"
print(f"Test AnnData on training common-gene space: {adata_common.shape}")


# -----------------------------------------------------------------------------
# 4. Compute total_counts on the common-gene matrix (matches step1a/step4)
# -----------------------------------------------------------------------------
sc.pp.calculate_qc_metrics(adata_common, percent_top=None,
                           log1p=False, inplace=True)
print("total_counts summary:")
print(adata_common.obs["total_counts"].describe())


# -----------------------------------------------------------------------------
# 5. Library-size normalize to TARGET_SUM, then log1p
#    (Same code path as step4: scale = target_sum / total_counts, no median.)
# -----------------------------------------------------------------------------
d = adata_common.obs["total_counts"].to_numpy()
scale = np.divide(TARGET_SUM, d, out=np.zeros_like(d, dtype=float), where=d > 0)

if issparse(adata_common.X):
    X_norm = adata_common.X.tocsr().multiply(scale[:, None]).tocsr()
else:
    X_norm = adata_common.X * scale[:, None]

# summary the read depth after normalization
total_counts_norm = X_norm.sum(axis=1).A1 if issparse(X_norm) else X_norm.sum(axis=1)
print("total_counts summary after normalization:")
print(pd.Series(total_counts_norm).describe())

adata_common.X = X_norm
sc.pp.log1p(adata_common)

print(adata_common.X.shape)
print(adata_common.X[:2, :4])

# -----------------------------------------------------------------------------
# 6. Subset to the top-1000 DE genes (model input order)
# -----------------------------------------------------------------------------
adata_top = adata_common[:, top_genes].copy()
assert list(adata_top.var_names) == top_genes, "Top-gene order mismatch!"

X_test = adata_top.X.toarray() if issparse(adata_top.X) else np.asarray(adata_top.X)
print(f"Model input matrix: {X_test.shape}")


# -----------------------------------------------------------------------------
# 7. Load the trained model and predict
# -----------------------------------------------------------------------------
model = load_model(MODEL_FILE)
out = model.predict(X_test)
pred = out["classification"].flatten()
print(f"Predictions: n={len(pred)}, mean={pred.mean():.3f}, "
      f"min={pred.min():.3f}, max={pred.max():.3f}")


# -----------------------------------------------------------------------------
# 8. Save predictions and a probability histogram
# -----------------------------------------------------------------------------
pred_df = adata_top.obs.copy()
pred_df["pred_prob"] = pred
pred_path = os.path.join(OUT_DIR, "predictions.csv")
pred_df.to_csv(pred_path)
print(f"Saved predictions to: {pred_path}")

plt.figure(figsize=(5, 3))
plt.hist(pred, bins=50, color="skyblue", edgecolor="black")
plt.xlabel("Predicted probability (cancer-reactive)")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pred_hist.png"), dpi=200)
plt.close()


# -----------------------------------------------------------------------------
# 9. Optional: evaluate against ground-truth labels if present
# -----------------------------------------------------------------------------

# take cells with antigen label as MANA as cancer-reactive
# and 'InfluenzaA' and 'EBV' as non-reactive
Counter(pred_df["antigen"])

pos_scores = []
neg_scores = []

for antg, score in zip(pred_df["antigen"], pred_df["pred_prob"]):
    if antg == "MANA":
        pos_scores += [score]
    elif antg in ["InfluenzaA", "EBV"]:
        neg_scores += [score]

both_scores = pos_scores + neg_scores
both_labels = [1] * len(pos_scores) + [0] * len(neg_scores)

fpr, tpr, _ = roc_curve(both_labels, both_scores)
roc_auc = auc(fpr, tpr)
print(f"ROC-AUC on test set Caushi et al. 2021: {roc_auc:.4f}")

plt.figure(figsize=(4, 4))
plt.plot(fpr, tpr, label=f"AUC = {roc_auc:.3f}")
plt.plot([0, 1], [0, 1], "k--", linewidth=1)
plt.xlabel("False positive rate")
plt.ylabel("True positive rate")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "roc_curve.png"), dpi=200)
plt.close()
