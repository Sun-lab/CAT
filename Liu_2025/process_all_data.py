# %%
# Dowload files: 

# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243013/suppl/GSE243013_NSCLC_immune_scRNA_counts.mtx.gz

# %%
import gzip
import pandas as pd
import os
import re
import pandas as pd
import numpy as np
import scanpy as sc
import gseapy as gp
import seaborn as sns
import matplotlib.pyplot as plt
import anndata
import shutil
import pickle
import scipy.io

target_rd = 3000
min_rd = 500
max_rd = 20000

# %%
metadata = pd.read_csv("GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")
print(metadata.shape)

# Select CD8T cells
cd8t_cells = metadata[
    (metadata['major_cell_type'] == 'T/NK cell') &
    (metadata['sub_cell_type'].str.startswith('CD8T'))
]

# Select CD4T cells
cd4t_cells = metadata[
    (metadata['major_cell_type'] == 'T/NK cell') &
    (metadata['sub_cell_type'].str.startswith('CD4T'))
]

print(f"CD8T cells: {cd8t_cells.shape[0]}")
print(f"CD4T cells: {cd4t_cells.shape[0]}")

# %%
print(metadata.head())

# %%
# Tabulate major_cell_type in metadata
major_counts = metadata['major_cell_type'].value_counts()
print("Major cell type counts:")
print(major_counts)

# Tabulate sub_cell_type among 'T/NK cell'
tnk_sub_counts = metadata.loc[metadata['major_cell_type'] == 'T/NK cell', 'sub_cell_type'].value_counts()
print("\nSub cell type counts within 'T/NK cell':")
print(tnk_sub_counts)

# %%
barcodes = pd.read_csv("GSE243013_barcodes.csv.gz")
print(barcodes.shape)
# Check if the barcodes and metadata.cellID are identical and in the same order
same_order = (barcodes['barcode'].values == metadata['cellID'].values).all()
print(f"Barcodes and metadata.cellID are identical and in the same order: {same_order}")

# %%
features = pd.read_csv("GSE243013_genes.csv.gz")
print(features.shape)
features.head()

# %%
X = scipy.io.mmread("GSE243013_NSCLC_immune_scRNA_counts.mtx.gz").tocsr()
print(X.shape)

# %%

# Set feature names (row names)
var = pd.DataFrame(index=features['geneSymbol'].tolist())
obs = pd.DataFrame(index=barcodes['barcode'].tolist())

# Create AnnData
adata = anndata.AnnData(X=X, obs=obs, var=var)
adata

adata.var["gene_name"] = features["geneSymbol"].values

# %%
with open("signatures_CD8.pkl", "rb") as f:
    sigs_CD8 = pickle.load(f)

print({k: len(v) for k, v in sigs_CD8.items()})

with open("signatures_CD4.pkl", "rb") as f:
    sigs_CD4 = pickle.load(f)

print({k: len(v) for k, v in sigs_CD4.items()})

del sigs_CD4['Jansen_TermDiff_73g']
del sigs_CD8['Lowery_neg_99g']

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

# %%
final_df = pd.DataFrame()
t_cell = {"CD4":[
'CD4T_Tem_GZMA', 
'CD4T_Tfh_CXCL13', 
'CD4T_Th1-like_CXCL13', 
'CD4T_Tm_ANXA1', 
'CD4T_Tm_XCL1', 
'CD4T_Tn_CCR7', 
'CD4T_Treg_CCR8', 
'CD4T_Treg_FOXP3', 
'CD4T_Treg_MKI67'], 
"CD8":[
'CD8T_ISG15', 
'CD8T_MAIT_KLRB1', 
'CD8T_NK-like_FGFBP2', 
'CD8T_prf_MKI67', 
'CD8T_Tem_GZMK+GZMH+', 
'CD8T_Tem_GZMK+NR4A1+', 
'CD8T_terminal_Tex_LAYN', 
'CD8T_Tex_CXCL13', 
'CD8T_Tm_IL7R', 
'CD8T_Trm_ZNF683'
]}
t_cell

# %%
adata.obs = adata.obs.join(metadata.set_index('cellID'))
adata.obs

sc.pp.calculate_qc_metrics(adata, inplace=True)

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

# Filter cells with total_counts < min_rd or > max_rd
adata = adata[(adata.obs['total_counts'] >= min_rd) & (adata.obs['total_counts'] <= max_rd)]

# %%
adata_all_CD4 = adata[adata.obs["sub_cell_type"].isin(t_cell['CD4'])]
adata_all_CD8 = adata[adata.obs["sub_cell_type"].isin(t_cell['CD8'])]

print(f"CD4T cells: {adata_all_CD4.shape[0]}")
print(f"CD8T cells: {adata_all_CD8.shape[0]}")

print(adata_all_CD8.X[:3, :3].toarray())

#%% Save the filtered data
import scipy.io
import subprocess

scipy.io.mmwrite("Liu_2025_CD4_counts.mtx", adata_all_CD4.X)
scipy.io.mmwrite("Liu_2025_CD8_counts.mtx", adata_all_CD8.X)

subprocess.run(["gzip", "-f", "Liu_2025_CD4_counts.mtx"])
subprocess.run(["gzip", "-f", "Liu_2025_CD8_counts.mtx"])

del adata  # Free memory

#%% 
sc.pp.normalize_total(adata_all_CD4, target_sum=target_rd)
sc.pp.log1p(adata_all_CD4)  # log transform

sc.pp.normalize_total(adata_all_CD8, target_sum=target_rd)
sc.pp.log1p(adata_all_CD8)  # log transform


# %%
final_df_CD8 = pd.DataFrame()
final_df_CD4 = pd.DataFrame()

for lib in t_cell['CD8']:
    print(lib)
    adata = adata_all_CD8[adata_all_CD8.obs["sub_cell_type"]==lib]
    print(adata)
    print(adata.var.head())

    sets_ave = ["Hanada_pos_27g", "Oliveira_virus_26g", "Hanada_neg_5g"]
    # Identify the rest for ssGSEA
    ssgsea_sets = {k: v for k, v in sigs_CD8.items() if k not in sets_ave}

    # gseapy.ssgsea expects a DataFrame where rows = genes, columns = samples.
    # So let's invert adata back: adata.X => n_cells x n_genes
    # We need genes x cells. We'll build a small data frame:

    all_cells = np.array(adata.obs.index)
    n_batches = 10
    batches = np.array_split(all_cells, n_batches)

    # Store results
    ssgsea_all_results = []

    ave_Hanada_pos_27g = []
    ave_Hanada_neg_5g  = []
    ave_Oliveira_virus_26g = []

    for i, batch_cells in enumerate(batches):
        print(f"Processing batch {i+1} / {n_batches} ...")

        # Extract gene expression matrix for this batch
        expr_batch = pd.DataFrame(
            adata[batch_cells].X.transpose().toarray(),  # genes x cells
            index=adata.var.index,
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
            columns=adata.var.index)

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

    # Create a DataFrame for the average signatures
    signature_df = pd.DataFrame(index=adata.obs.index)
    signature_df["CD8_ave_Hanada_pos_27g"] = ave_Hanada_pos_27g
    signature_df["CD8_ave_Hanada_neg_5g"]  = ave_Hanada_neg_5g
    signature_df["CD8_ave_Oliveira_virus_26g"] = ave_Oliveira_virus_26g

    df_combined = df_wide.join(signature_df, how="inner")
    final_df_CD8 = pd.concat([final_df_CD8, df_combined], ignore_index=False)


for lib in t_cell['CD4']:
    print(lib)
    adata = adata_all_CD4[adata_all_CD4.obs["sub_cell_type"]==lib]
    print(adata)
    print(adata.var.head())

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
            index=adata.var.index,
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
            columns=adata.var.index)

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

    # Create a DataFrame for the average signatures
    signature_df = pd.DataFrame(index=adata.obs.index)
    signature_df["CD4_ave_Hanada_pos_9g"] = ave_Hanada_pos_9g
    signature_df["CD4_ave_Hanada_neg_4g"] = ave_Hanada_neg_4g

    df_combined = df_wide.join(signature_df, how="inner")
    final_df_CD4 = pd.concat([final_df_CD4, df_combined], ignore_index=False)


# Save the final DataFrames to CSV files
final_df_CD4.to_csv('output_CD4.csv',index = True)
final_df_CD8.to_csv('output_CD8.csv',index = True)
