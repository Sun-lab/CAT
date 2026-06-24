# CAT

Cancer Associated T cells — identification of cancer-reactive (CR) T cells from
single-cell RNA-seq, and neural-network (NN) models that predict CR status from
gene expression.

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

## NN models for whether a T cell is cancer-associated or not

- **classify_cancer_reactive_T**

NN models were trained for predicting whether a T cell is cancer-associated or not based on gene experssion. For each cell, the output is a numerical score between 0 and 1. One model was trained for CD4 and one model was for CD8.

### Build the conda environment

```
conda env create -f environment.yml
conda activate CAT
```

### Make prediction on test data

The test data should be T cells with gene expression count. The NN model for the corresponding cell type should be used. For example, for CD8 T cells, the NN model trained for CD8 should be used.

To prepare the test data for making predictions on, these steps are needed:

1. restrict test genes to a set of genes, which are the genes on which the model development data library-size normalize was run. For the genes that do not exist in the test data, sert zero columns.
2. compute per-cell total_counts on the set of genes in point 1.
3. library-size normalize to target_sum = 3000, then log1p
4. subset to the top-1000 DE genes used as model input, arranged in the same order as in model input

and these steps are all included in an example tutorial for using the CD8 NN model to make predictions on independent test data:

[step5_predict_on_Caushi_data_CD8.py](https://github.com/Sun-lab/CAT/blob/main/classify_cancer_reactive_T/step5_predict_on_Caushi_data_CD8.py)
