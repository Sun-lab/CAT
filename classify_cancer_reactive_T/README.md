# `classify_cancer_reactive_T` — neural-network classifier for cancer-reactive (CR) T cells

This folder trains and evaluates the neural-network (NN) models that predict whether a
CD8⁺ or CD4⁺ T cell is **cancer-reactive (CR)** from its transcriptome, and applies the
trained models to external datasets.

The CR training labels themselves are produced **upstream** (per-study processing in the
`../<Study>/` folders, then combined and refined in `../zz1_prepare_data.Rmd` and
`../zz2_construct_training_data.Rmd`). This folder takes those labels as input and handles
gene selection, model training, leave-one-study-out cross-validation (LOSO-CV), and
prediction.

Five studies are used in production: **Zheng_2021, Chen_2024, Chow_2023, Liu_2022,
Liu_2025** (Chow_2023 contributes CD8 only; it has no usable CD4 CR set, so CD4 LOSO-CV
skips it).

---

## 1. Gene selection — a three-stage cascade

The NN does **not** see the genome. Input features are reduced in three stages, and it is
important to keep them distinct:

| Stage | Genes | How it is obtained | Where in code |
|------|-------|--------------------|---------------|
| **(a) Common genes** | **16,323** | Intersection of the detected genes across all 5 studies | `step1a_prepare_CD8.ipynb` / `step1b_prepare_CD4.ipynb` (`set.intersection(...)`) |
| **(b) Prevalence-filtered genes** | **6,256 (CD8)** / **6,241 (CD4)** | Of the common genes, keep those expressed (non-zero count) in **≥ 5% of pooled cells** | `step1a/step1b`, code block 23 (`keep_genes = gene_expr_frac >= 0.05`) |
| **(c) Top-1,000 DE genes (model input)** | **1,000** | Rank genes by Wilcoxon differential expression (CR=1 vs non-CR=0) and take the top 1,000 — selected **within the 6,256/6,241 filtered set**, not the 16,323 common genes | `step2_classify_careT.ipynb` (per CV fold) and `step4_classify_careT_{CD8,CD4}.ipynb` (final model) |

Notes on stage (b): block 22 prints the sweep used to choose the 5% cutoff — for CD8,
≥1%→9,510, ≥2%→8,514, ≥3%→7,656, **≥5%→6,256** genes; for CD4, ≥1%→9,507, ≥2%→8,465,
≥3%→7,592, **≥5%→6,241** genes. The filtered matrix is saved as
`../data/{CD8,CD4}_combined_filtered.h5ad` and reused by all downstream steps. The
prevalence filter is computed **once on all five studies**, so the 6,256/6,241-gene
universe is the same in every CV fold; only the DE ranking in stage (c) is recomputed per
fold.

Notes on stage (c): in LOSO-CV the DE ranking uses only the four **training** studies
(`adata_sub = adata[adata.obs["study"] != test_study]`), so feature selection never sees
the held-out study. In the final model the DE ranking pools **all five** studies. The
gene lists are saved (`DE_results/…` and `model/{CD8,CD4}_top1000_genes.txt`).

---

## 2. Pipeline and purpose of each file

Run order: **step1 → step2 (+ step2 runner) → step3 → step4 → step5**.

### Notebooks / scripts

| File | Purpose |
|------|---------|
| `step1a_prepare_CD8.ipynb` | Build the CD8 training matrix: load per-study raw CD8 counts, intersect to the **16,323** common genes, apply the **≥5% prevalence filter → 6,256 genes**, attach CR labels, and save `../data/CD8_combined_filtered.h5ad`. **(Defines gene stages a–b.)** |
| `step1b_prepare_CD4.ipynb` | Same as above for CD4 → **6,241 genes**, saved as `../data/CD4_combined_filtered.h5ad`. |
| `step2_classify_careT.ipynb` | **LOSO-CV template** (parameterized by `test_study`). Loads the filtered h5ad, runs DE on the four training studies, selects the **top-1,000 genes for that fold**, trains the autoencoder+classifier, predicts on the held-out study, and writes the fold's predictions/figures. Not run directly — executed by the runner below. |
| `step2_run_classify_careT.py` | **LOSO-CV driver.** Uses `papermill` to execute `step2_classify_careT.ipynb` once per (study × cell type), producing the `CD8_<study>_…ipynb` / `CD4_<study>_…ipynb` executed notebooks. Skips CD4 for Chow_2023. |
| `CD8_<study>_lr_0.001_bs_32.ipynb`, `CD4_<study>_lr_0.001_bs_32.ipynb` | **Executed CV-fold notebooks** (one per held-out study), the rendered output of the runner. Each holds that fold's DE, training curves, and held-out ROC/predictions. |
| `step3_update_CD8_anno.ipynb`, `step3_update_CD4_anno.ipynb` | **Aggregate the LOSO-CV predictions.** Collect each fold's `predictions_test_data.txt` (each cell scored by the model that did *not* see its study), merge them back onto the cells, and write `../data/{CD8,CD4}_with_NN_predictions.tsv`. These out-of-sample scores (threshold 0.5) give the honest, final CR annotation. |
| `step4_classify_careT_CD8.ipynb`, `step4_classify_careT_CD4.ipynb` | **Train the final, deployable model** on all five studies. DE/top-1,000 selection is pooled across studies; saves `model/{CD8,CD4}.keras` and `model/{CD8,CD4}_top1000_genes.txt`. |
| `step5_predict_on_Caushi_data_CD8.py` | **Apply a trained model to a new/external dataset** (example: Caushi 2021). Reindexes test data to the common-gene space, normalizes (target depth 3,000, log1p), subsets to the saved top-1,000 genes (in order), loads `model/CD8.keras`, and writes per-cell predictions + ROC. |

### Output subfolders

| Folder | Contents |
|--------|----------|
| `model/` | Final deployable artifacts: `CD8.keras`, `CD4.keras`, the `{CD8,CD4}_top1000_genes.txt` input-gene lists, and `{CD8,CD4}_learn_rate_…/` training-curve plots/tables. |
| `DE_results/` | Wilcoxon DE tables. `DE_{CD8,CD4}.csv` = pooled (final model); `DE_{CD8,CD4}_test_study_<study>.csv` = per LOSO-CV fold. The top 1,000 rows (by p-value) are the model-input genes for that setting. |
| `classification_result/` | Per CV-fold outputs (`{CT}_<study>_…/`): held-out `predictions_test_data.txt`, `roc_curve.pdf`, box plots, training curves. Consumed by `step3`. |
| `figures/` | LOSO-CV logistic-regression ROC plots, DE volcano plots (pooled and per fold), and `umapfigures/` (UMAPs colored by pos/neg score, study, tissue). |

---

## 3. Model architecture (summary)

The classifier is a joint **autoencoder + classification head** taking the 1,000-gene
log-normalized expression vector as input, trained with combined reconstruction +
classification loss. Default hyperparameters (encoded in folder names): learning rate
0.001, batch size 32, classification/reconstruction loss weights 1.0/1.0, up to 20 epochs.
A T cell is annotated CR if its prediction score (0–1) exceeds **0.5**. See the main-text
Methods for full details.

---

## 4. Reproduce

```bash
conda activate CAT                  # see ../environment.yml

# 1) Build filtered training matrices (gene stages a–b)
#    run step1a_prepare_CD8.ipynb and step1b_prepare_CD4.ipynb

# 2) Leave-one-study-out cross-validation (per-fold gene stage c + training)
python step2_run_classify_careT.py

# 3) Aggregate out-of-sample CV predictions into final CR annotation
#    run step3_update_CD8_anno.ipynb and step3_update_CD4_anno.ipynb

# 4) Train the final deployable models (pooled gene stage c)
#    run step4_classify_careT_CD8.ipynb and step4_classify_careT_CD4.ipynb

# 5) Apply a trained model to external data
python step5_predict_on_Caushi_data_CD8.py
```
