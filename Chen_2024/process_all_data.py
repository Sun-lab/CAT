# %%
import os
import pandas as pd
import scanpy as sc
import numpy as np

# %%
# set target read depth for normalization
target_rd = 3000

# Define thresholds for filtering
min_rd = 500
max_rd = 20000

features = pd.read_csv("GSE236581_features.tsv.gz", compression='gzip', sep='\t', header=None)
barcodes = pd.read_csv("GSE236581_barcodes.tsv.gz", compression='gzip', header=None)
features.head(3)
barcodes.head(3)

from scipy.io import mmread

import gzip
with gzip.open("GSE236581_counts.mtx.gz", 'rb') as f:
    counts = mmread(f).tocsc()  # or .tocsr()
print(counts.shape)

import anndata
features.columns = ["gene_id", "gene_name", "feature_type"]
barcodes.columns = ["barcode"]

counts = counts.transpose()  

# Create the AnnData object
adata = anndata.AnnData(X=counts)

# AnnData “obs” (rows) will correspond to cells
adata.obs_names = barcodes["barcode"].values

# AnnData “var” (columns) will correspond to genes
adata.var_names = features["gene_id"].values

# Store additional gene annotations
adata.var["gene_name"] = features["gene_name"].values
adata.var["feature_type"] = features["feature_type"].values

adata

#%%
import pickle

with open("signatures_CD8.pkl", "rb") as f:
    sigs_CD8 = pickle.load(f)

print({k: len(v) for k, v in sigs_CD8.items()})

with open("signatures_CD4.pkl", "rb") as f:
    sigs_CD4 = pickle.load(f)

print({k: len(v) for k, v in sigs_CD4.items()})

del sigs_CD4['Jansen_TermDiff_73g']
del sigs_CD8['Lowery_neg_99g']

metadata = pd.read_csv("GSE236581_CRC-ICB_metadata.txt.gz", compression='gzip', sep=' ')
metadata.shape
metadata.head(3)

#%%
def average_genes(genes_list, df):
    valid_genes = [g for g in genes_list if g in df.columns]
    if len(valid_genes) == 0:
        return np.repeat(np.nan, df.shape[0])
    sub_df = df[valid_genes]
    return sub_df.mean(axis=1)  # mean over columns => per-cell average

def match_genes_in_sets(sigs_CD8, se2):    
    # Dictionary to hold matched genes for each set
    gs = {}
    # Loop over each key in the sigs_CD8 dictionary
    for s1 in sigs_CD8.keys():
        # 'sigs_CD8[s1]' is a list of genes
        # We'll find which of these genes appear in se2.index
        matched_genes = []
        
        for gene in sigs_CD8[s1]:
            if gene in se2['gene_name'].tolist():
                matched_genes.append(gene)
        
        # Print match info
        print(f"{s1}: {len(matched_genes)}/{len(sigs_CD8[s1])} genes are found.")
        
        # Store in gs
        gs[s1] = matched_genes
    return

import gseapy as gp
final_df_CD8 = pd.DataFrame()
final_df_CD4 = pd.DataFrame()
t_cell = {"CD4":[
 'c01_CD4_Tn_CCR7',
 'c02_CD4_Tn_SELL',
 'c03_CD4_Tn_NR4A2',
 'c04_CD4_Tcm_ANXA1',
 'c05_CD4_Tcm_GPR183',
 'c06_CD4_Trm_HSPA1A',
 'c07_CD4_Th17_CTSH',
 'c08_CD4_Tfh_CXCL13_IL6ST',
 'c09_CD4_Th1_CXCL13_HAVCR2',
 'c10_CD4_Temra_GZMB',
 'c11_CD4_Treg_FOXP3',
 'c12_CD4_Treg_KLRB1',
 'c13_CD4_Treg_TNFRSF9',
 'c14_CD4_MT'], 'CD8':[
 'c15_CD8_Tn_CCR7',
 'c16_CD8_Tn_SELL',
 'c17_CD8_Tcm_GPR183',
 'c18_CD8_Tcm_ANXA1',
 'c19_CD8_Tem_CMC1',
 'c20_CD8_Tem_GZMK',
 'c21_CD8_Trm_XCL1',
 'c22_CD8_Trm_HSPA1B',
 'c23_CD8_Tex_LAYN',
 'c24_CD8_Temra_CX3CR1',
 'c25_CD8_Temra_TYROBP',
 'c26_CD8_MAIT_KLRB1',
 'c27_CD8_MAIT_SLC4A10',
 'c28_CD8_IEL_CD160',
]}

sc.pp.calculate_qc_metrics(adata, inplace=True)
adata.obs = adata.obs.join(metadata)

# Summarize percentiles and min/max for adata
print("Percentiles (25%, 50%, 75%) of total_counts for adata:", 
    np.percentile(adata.obs['total_counts'].values, [25, 50, 75]))
print("Min and max of total_counts for adata:", 
    np.min(adata.obs['total_counts'].values), 
    np.max(adata.obs['total_counts'].values))

# Check percentage of cells with counts < min_rd or > max_rd
low_count = (adata.obs['total_counts'] < min_rd).sum()
high_count = (adata.obs['total_counts'] > max_rd).sum()
total_cells = adata.shape[0]
print(f"Percentage of cells with counts < {min_rd}: {100 * low_count / total_cells:.2f}%")
print(f"Percentage of cells with counts > {max_rd}: {100 * high_count / total_cells:.2f}%")

#%%

# Filter cells with total_counts < min_rd or > max_rd
adata = adata[(adata.obs['total_counts'] >= min_rd) & (adata.obs['total_counts'] <= max_rd)]

adata_all_CD4 = adata[adata.obs["SubCellType"].isin(t_cell['CD4'])]
adata_all_CD8 = adata[adata.obs["SubCellType"].isin(t_cell['CD8'])]


adata_all_CD4
adata_all_CD8

print(adata_all_CD8.X[:3, :3].toarray())

import scipy.io
import subprocess

scipy.io.mmwrite("Chen_2024_CD4_counts.mtx", adata_all_CD4.X)
scipy.io.mmwrite("Chen_2024_CD8_counts.mtx", adata_all_CD8.X)

subprocess.run(["gzip", "-f", "Chen_2024_CD4_counts.mtx"])
subprocess.run(["gzip", "-f", "Chen_2024_CD8_counts.mtx"])

del adata  # Free memory

#%%
sc.pp.normalize_total(adata_all_CD4, target_sum=target_rd)
sc.pp.log1p(adata_all_CD4)  # log transform

sc.pp.normalize_total(adata_all_CD8, target_sum=target_rd)
sc.pp.log1p(adata_all_CD8)  # log transform

#%%
for lib in t_cell['CD8']:
    print(lib)
    adata = adata_all_CD8[adata_all_CD8.obs["SubCellType"]==lib]
    print(adata)

    sets_ave = ["Hanada_pos_27g", "Oliveira_virus_26g", "Hanada_neg_5g"]
    # Identify the rest for ssGSEA
    ssgsea_sets = {k: v for k, v in sigs_CD8.items() if k not in sets_ave}
    # gseapy.ssgsea expects a DataFrame where rows = genes, columns = samples.
    # So let's invert adata back: adata.X => n_cells x n_genes
    # We need genes x cells. We'll build a small data frame:
    
    all_cells = np.array(adata.obs.index)
    n_batches = 10
    batches = np.array_split(all_cells, n_batches)
    print([len(batch) for batch in batches])
    # Check that the sum of batch lengths equals the total number of cells
    assert sum(len(batch) for batch in batches) == len(all_cells), \
        f"Batching error: sum={sum(len(batch) for batch in batches)}, expected={len(all_cells)}"

    # Store results
    ssgsea_all_results = []

    ave_Hanada_pos_27g = []
    ave_Hanada_neg_5g  = []
    ave_Oliveira_virus_26g = []

    for i, batch_cells in enumerate(batches):
        print(f"Processing batch {i+1} / {n_batches} ...")

        expr_batch = pd.DataFrame(
            adata[batch_cells].X.transpose().toarray(),  # genes x cells
            index=adata.var['gene_name'].tolist(),   # Genes
            columns=batch_cells
        )
        result = gp.ssgsea(data=expr_batch,
                        gene_sets=ssgsea_sets,
                        sample_norm=False,
                        outdir=None,
                        min_size=1,
                        max_size=20000)
        result.run()

        ssgsea_all_results.append(result.res2d)

        expr_batch = pd.DataFrame(
            adata[batch_cells].X.toarray(), 
            index=batch_cells, 
            columns=adata.var['gene_name'].tolist()
        )

        aa1 = average_genes(sigs_CD8["Hanada_pos_27g"], expr_batch)
        ave_Hanada_pos_27g = np.concatenate([ave_Hanada_pos_27g, aa1.values])

        aa2 = average_genes(sigs_CD8["Hanada_neg_5g"], expr_batch)
        ave_Hanada_neg_5g = np.concatenate([ave_Hanada_neg_5g, aa2.values])

        aa3 = average_genes(sigs_CD8["Oliveira_virus_26g"], expr_batch)
        ave_Oliveira_virus_26g = np.concatenate([ave_Oliveira_virus_26g, aa3.values])

    # Combine all results
    combined_ssgsea = pd.concat(ssgsea_all_results, axis=0)

    # Transpose to cell x set
    df_wide = combined_ssgsea.pivot(
        index='Name',    # Each sample’s name/ID
        columns='Term',  # The gene set name
        values='NES'      # Which score you want to spread out as columns
    )
    # rename columns
    df_wide.columns = [f"CD8_{col}" for col in df_wide.columns]

    # Now combine them with ssgsea_scores
    signature_df = pd.DataFrame(index=adata.obs.index)
    signature_df["CD8_ave_Hanada_pos_27g"] = ave_Hanada_pos_27g
    signature_df["CD8_ave_Hanada_neg_5g"]  = ave_Hanada_neg_5g
    signature_df["CD8_ave_Oliveira_virus_26g"] = ave_Oliveira_virus_26g

    df_combined = df_wide.join(signature_df, how="inner")
    final_df_CD8 = pd.concat([final_df_CD8, df_combined], ignore_index=False)

for lib in t_cell['CD4']:
    print(lib)
    adata = adata_all_CD4[adata_all_CD4.obs["SubCellType"]==lib]
    print(adata)

    sets_ave_CD4 = ["Hanada_pos_9g", "Hanada_neg_4g"]
    # Identify the rest for ssGSEA
    ssgsea_sets = {k: v for k, v in sigs_CD4.items() if k not in sets_ave_CD4}
    # gseapy.ssgsea expects a DataFrame where rows = genes, columns = samples.
    # So let's invert adata back: adata.X => n_cells x n_genes
    # We need genes x cells. We'll build a small data frame:

    all_cells = np.array(adata.obs.index)
    n_batches = 10
    batches = np.array_split(all_cells, n_batches)

    # Store results
    ssgsea_all_results = []

    ave_Hanada_pos_9g = []
    ave_Hanada_neg_4g = []

    for i, batch_cells in enumerate(batches):
        print(f"Processing batch {i+1} / {n_batches} ...")

        # Extract gene expression matrix for this batch
        expr_batch = pd.DataFrame(
            adata[batch_cells].X.transpose().toarray(),  # genes x cells
            index=adata.var['gene_name'].tolist(),
            columns=batch_cells
        )

        # Run ssGSEA
        result = gp.ssgsea(data=expr_batch,
                        gene_sets=ssgsea_sets,
                        sample_norm=False,
                        outdir=None,
                        min_size=1,
                        max_size=20000)
        result.run()

        # Store result
        ssgsea_all_results.append(result.res2d)

        expr_batch = pd.DataFrame(
            adata[batch_cells].X.toarray(), 
            index=batch_cells, 
            columns=adata.var['gene_name'].tolist())

        aa1 = average_genes(sigs_CD4["Hanada_pos_9g"], expr_batch)
        ave_Hanada_pos_9g = np.concatenate([ave_Hanada_pos_9g, aa1.values])   
        aa2 = average_genes(sigs_CD4["Hanada_neg_4g"], expr_batch)
        ave_Hanada_neg_4g = np.concatenate([ave_Hanada_neg_4g, aa2.values])

    combined_ssgsea = pd.concat(ssgsea_all_results, axis=0)
    combined_ssgsea

    # Transpose to cell x set
    df_wide = combined_ssgsea.pivot(
        index='Name',    # Each sample’s name/ID
        columns='Term',  # The gene set name
        values='NES'      # Which score you want to spread out as columns
    )

    # rename columns
    df_wide.columns = [f"CD4_{col}" for col in df_wide.columns]

    signature_df = pd.DataFrame(index=adata.obs.index)
    signature_df["CD4_ave_Hanada_pos_9g"] = ave_Hanada_pos_9g
    signature_df["CD4_ave_Hanada_neg_4g"] = ave_Hanada_neg_4g

    df_combined = df_wide.join(signature_df, how="inner")
    final_df_CD4 = pd.concat([final_df_CD4, df_combined], ignore_index=False)

final_df_CD4.to_csv('output_CD4.csv',index = True)
final_df_CD8.to_csv('output_CD8.csv',index = True)
