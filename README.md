# CAT

Cancer Associated T cells — identification of cancer-reactive (CR) T cells from
single-cell RNA-seq, and neural-network (NN) models that predict CR status from
gene expression.

There are two ways to use this repository:

- **Annotate your own T cells** with the trained models (CAT-NN). Go to
  [Annotate your own T cells with CAT-NN](#annotate-your-own-t-cells-with-cat-nn).
  This needs only the trained models and two gene lists, which are all included
  here. No training data and no re-training are needed.
- **Reproduce the atlas and the model training** of the paper. See
  [Repository structure](#repository-structure) and
  [Top-level workflow scripts](#top-level-workflow-scripts-zz).

## Repository structure

**Study folders** — each processes one published scRNA-seq dataset: restrict to the
common gene set, normalize, score the cancer-reactive gene signatures (ssGSEA),
cluster, select cancer-reactive cells, and integrate TCR data. The typical entry
points in each folder are `process_all_data.py` (signature scoring) followed by
`process_output.ipynb` (clustering, CR-cell selection, TCR metadata).
- **Used in the production training set:** `Zheng_2021`, `Chen_2024`, `Chow_2023`,
  `Liu_2022`, `Liu_2025`.
- **Additional datasets** (processed, not used in the final models):
  `Chin_2025`, `Heiser_2023`, `Li_2025`, `Upadhye_2023`, `Wang_2024`.

**`classify_cancer_reactive_T/`** — training, leave-one-study-out cross-validation,
and prediction of the CD4/CD8 cancer-reactive NN classifiers (gene-selection
cascade → DE top-1000 genes → autoencoder+classifier). Has its own README.

**`data/`** — shared resources: `common_genes.txt` (16,323 genes common to all
studies), `signatures_CD4.pkl` / `signatures_CD8.pkl` (cancer-reactive gene
signatures used for ssGSEA), and `scale_CD4.tsv` / `scale_CD8.tsv` (scaling
parameters). Large generated tables (e.g. `*_combined.tsv.gz`) are written to the
working data directory.

**`_checking_results/`** — auxiliary QC / exploration scripts (e.g.
`analysis_single_positive.py`, `combine_single_positive.py`, `make_word_table.py`,
`explore_files.ipynb`).

**Other:** `environment.yml` (conda environment), `LICENSE`, `_summary.xlsx`.

## Top-level workflow scripts (`zz*`)

These run **after** the per-study folders and tie the datasets together. Run in order:
They are used to build the atlas and the training data. **They are not needed to annotate new data.**

- **`zz0_get_common_gene_set.ipynb`** — derive the set of genes shared across all
  studies (the 16,323-gene common space used throughout, saved to
  `data/common_genes.txt`).
- **`zz1_prepare_data.Rmd`** (`.html` = rendered report) — combine the per-study
  cleaned cell + TCR tables into harmonized `CD4_combined.tsv.gz` /
  `CD8_combined.tsv.gz`, standardizing TCR column names and cell metadata across
  studies.
- **`zz2_construct_training_data.Rmd`** (`.html` = rendered report) — build the
  refined CR labels for NN training: define `cancer_reactive` = per-cell signature
  call **AND** CR-cluster membership; keep clonally expanded CR TCRs; remove TCRs
  matching **non-human antigens** via VDJdb. Outputs `CD4_updated.tsv(.gz)` /
  `CD8_updated.tsv(.gz)`.
- **`zz3_db_compare.py`** and **`zz3_db_compare.ipynb`** — database comparison for
  the non-human-antigen TCR-removal step: how many TCRs/cells **VDJdb vs McPAS vs
  IEDB** would each filter, how they overlap, and whether adding McPAS/IEDB would
  change the cancer-reactive set. The `.py` is an importable module + command-line
  entry point; the `.ipynb` runs the same analysis with all tables and figures
  embedded inline for viewing. (`zz4` imports the `.py`, so keep both.)
- **`zz4_clonality_tables.py`** — summarize the clonality (clone sizes) of cells
  filtered by IEDB / McPAS, stratified by whether the TCR is shared with VDJdb,
  and export the CD8/CD4 tables to a Word document. Imports `zz3_db_compare`.

<!--
# Output files of the per-study pipelines:

**cell_meta_data_CD4.csv, cell_meta_data_CD8.csv**: contains scores for each cell, basic cell information, whether a cell is cancer reactive or not, TCR (could have multiple) associated with the cell. These files can be multi-indexed if TCR info is available for each cell.

**cell_meta_data_CD4_cleaned.csv, cell_meta_data_CD8_cleaned.csv**: contains scores for each cell, basic cell information, whether a cell is cancer reactive or not, and one alpha chain and one beta chain from the TCR data (not multi-indexed). If more than 1 alpha/beta chain available for a single cell, we only take the one with highest UMI.

**tcr_meta_data_CD4_cells.csv, tcr_meta_data_CD8_cells.csv**: TCR data for alpha and beta chains including v gene, j gene, cdr3, and cell level information (average score, num_cells have this TCR, num_cancer_reactive_cells have this TCR).

## Run it using local conda environment

```
conda create --name CAT python=3.10.8
conda activate CAT
conda install -c conda-forge pandas numpy scanpy anndata seaborn matplotlib
conda install -c conda-forge python-igraph leidenalg
```

First run ```zz0_get_common_gene_set.ipynb``` to get the common set of genes to be used, and then run the codes within each study-specific folder. Typically to run ```process_all_data.py``` followed by ```process_output.ipynb```. Finally run ```zz1_prepare_data.Rmd``` and ```zz2_construct_training_data.Rmd``` in this folder.
-->

## Annotate your own T cells with CAT-NN

CAT-NN is a neural network that predicts whether a T cell is tumor-reactive
(called cancer-reactive, CR, elsewhere in this repository) from its gene
expression. The output is a score between 0 and 1 for each cell. There is one
model for CD8+ T cells and one for CD4+ T cells.

### What you need

**Your data:** single-cell RNA-seq data of CD8+ T cells or CD4+ T cells, as
**raw UMI counts** (not normalized) with **gene symbols**. Select the T cells and
separate CD8+ from CD4+ T cells before using CAT-NN, and use the model that
matches the cell type.

**Files in this repository:**

| File | Content |
| --- | --- |
| `classify_cancer_reactive_T/model/CD8.keras`, `CD4.keras` | Trained CAT-NN models |
| `classify_cancer_reactive_T/model/CD8_top1000_genes.txt`, `CD4_top1000_genes.txt` | The 1,000 input genes of each model, in the order used by the model |
| `data/common_genes.txt` | The 16,323 genes shared by the five training studies, used to calculate the read depth of each cell |

**Environment:**

```
conda env create -f environment.yml
conda activate CAT
```

The models were saved with Keras 3 (TensorFlow 2.21). An older version of Keras
may fail to load them.

### Steps

The input of CAT-NN has to be prepared in the same way as the training data:

1. Calculate the read depth of each cell as its total count over the 16,323
   common genes. A common gene that is absent from your data counts as 0.
2. Extract the 1,000 input genes in the order used by the model. An input gene
   that is absent from your data is set to 0.
3. Normalize the counts to 3,000 per cell using the read depth of step 1, and
   apply log(1 + x).
4. Apply the model.

The function below does these four steps.

```python
import numpy as np
import pandas as pd
from scipy import sparse
from keras.models import load_model

def predict_cat_nn(adata, cell_type="CD8", cat_dir="."):
    """Return the CAT-NN score (0-1) of each cell.

    adata     : AnnData of CD8+ (or CD4+) T cells with raw UMI counts in adata.X
                and unique gene symbols in adata.var_names
    cell_type : "CD8" or "CD4"
    cat_dir   : path to the CAT repository
    """
    model_dir = f"{cat_dir}/classify_cancer_reactive_T/model"
    common = pd.read_csv(f"{cat_dir}/data/common_genes.txt", header=None)[0]
    top = pd.read_csv(f"{model_dir}/{cell_type}_top1000_genes.txt", header=None)[0]
    model = load_model(f"{model_dir}/{cell_type}.keras", compile=False)

    X = sparse.csr_matrix(adata.X)
    genes = pd.Index(adata.var_names)
    # 1. read depth = total count of the 16,323 common genes (absent genes count as 0)
    depth = np.asarray(X[:, genes.isin(common)].sum(axis=1)).ravel()
    # 2. the 1,000 input genes in the order used by the model; absent genes are set to 0
    pos = genes.get_indexer(top)
    X_top = np.zeros((adata.n_obs, len(top)), dtype=np.float32)
    X_top[:, pos >= 0] = X[:, pos[pos >= 0]].toarray()
    print(f"{(pos < 0).sum()} of the {len(top)} input genes are absent from the data")
    # 3. normalize to 3,000 counts per cell, then log(1 + x)
    X_top = np.log1p(X_top * (3000 / np.maximum(depth, 1))[:, None])
    return model.predict(X_top, verbose=0)["classification"].ravel()
```

Example:

```python
import scanpy as sc

adata = sc.read_h5ad("my_CD8_T_cells.h5ad")   # raw UMI counts
adata.var_names_make_unique()
adata.obs["CAT_NN_score"] = predict_cat_nn(adata, cell_type="CD8", cat_dir="path/to/CAT")
adata.obs["tumor_reactive"] = adata.obs["CAT_NN_score"] > 0.5
```

### Output

- A cell is annotated as tumor-reactive if its score is larger than 0.5, the
  cutoff used in the paper.
- The score of a TCR clone is the average score of the cells of the clone.

### Notes

- Use raw counts. If your data are already normalized, the scores are not valid.
- The read depth is the total over the common genes, not over all genes of your
  data and not over the 1,000 input genes.
- The function prints the number of input genes that are absent from your data.
  If many are absent (for example, with a targeted gene panel), the scores are
  less reliable.
- In the training data, cells with fewer than 500 or more than 20,000 UMIs were
  removed.

### A complete example

[classify_cancer_reactive_T/step5_predict_on_Caushi_data_CD8.py](https://github.com/Sun-lab/CAT/blob/main/classify_cancer_reactive_T/step5_predict_on_Caushi_data_CD8.py)
applies the CD8 model to an independent dataset (Caushi et al. 2021), following
the same four steps, and evaluates the scores against the antigen specificity
of the T cells. The data of this example are not included in this repository:
the script reads a count matrix, gene names, barcodes and cell metadata from
`data/Caushi_2021/for_python`.

For how the models were trained, see
[classify_cancer_reactive_T/README.md](classify_cancer_reactive_T/README.md).
