#%%
import pandas as pd
import scanpy as sc
import numpy as np
import os

target_rd = 3000
min_rd = 500
max_rd = 20000

#%%

metadata_df = pd.read_csv("GSE179994_Tcell.metadata.tsv.gz", sep="\t", compression="gzip")
print(metadata_df.shape)
metadata_df.head()

# Load the sparse matrix
adata = sc.read_mtx("GSE179994_all.Tcell.rawCounts.mtx")
print(adata.shape)
adata.obs_names = [line.strip() for line in open("barcodes.txt")]
adata.var_names = [line.strip() for line in open("features.txt")]   

metadata_df = metadata_df.set_index('cellid')

adata.obs.index.name = 'cellid'

adata.obs = adata.obs.join(metadata_df, how='left')

import pickle

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
            if gene in se2.index.tolist():
                matched_genes.append(gene)
        
        # Print match info
        print(f"{s1}: {len(matched_genes)}/{len(sigs_CD8[s1])} genes are found.")
        
        # Store in gs
        gs[s1] = matched_genes
    return

import gseapy as gp
final_df_CD8 = pd.DataFrame()
final_df_CD4 = pd.DataFrame()

t_cell_types = [
'CD4',
'CD8'
]
print(adata.obs["celltype"].value_counts(dropna=False))

adata = adata[adata.obs["celltype"].isin(t_cell_types)]
sc.pp.calculate_qc_metrics(adata, inplace=True)

#%%
# filter adata by read depth
print(adata.shape)
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
print(adata.shape)

adata_all_CD4 = adata[adata.obs["celltype"]=='CD4']
adata_all_CD8 = adata[adata.obs["celltype"]=='CD8']

print(adata_all_CD8.X[:3, :3].toarray())

#%% Save the filtered data

import scipy.io
import subprocess

os.makedirs("Liu_2022_CD4", exist_ok=True)
os.makedirs("Liu_2022_CD8", exist_ok=True)

scipy.io.mmwrite("Liu_2022_CD4/matrix.mtx", adata_all_CD4.X)
scipy.io.mmwrite("Liu_2022_CD8/matrix.mtx", adata_all_CD8.X)

subprocess.run(["gzip", "-f", "Liu_2022_CD4/matrix.mtx"])
subprocess.run(["gzip", "-f", "Liu_2022_CD8/matrix.mtx"])

# Save gene names (variables)
adata_all_CD4.var_names.to_series().to_csv("Liu_2022_CD4/genes.tsv", 
                                           sep='\t', index=False, header=False)
adata_all_CD8.var_names.to_series().to_csv("Liu_2022_CD8/genes.tsv", 
                                           sep='\t', index=False, header=False)

# Save cell names (observations)
adata_all_CD4.obs_names.to_series().to_csv("Liu_2022_CD4/barcodes.tsv", 
                                           sep='\t', index=False, header=False)
adata_all_CD8.obs_names.to_series().to_csv("Liu_2022_CD8/barcodes.tsv", 
                                           sep='\t', index=False, header=False)

del adata  # Free memory

#%% 
sc.pp.normalize_total(adata_all_CD4, target_sum=target_rd)
sc.pp.log1p(adata_all_CD4)  # log transform

sc.pp.normalize_total(adata_all_CD8, target_sum=target_rd)
sc.pp.log1p(adata_all_CD8)  # log transform

#%% CD8
adata = adata_all_CD8.copy()
print(adata)
sets_ave = ["Hanada_pos_27g", "Oliveira_virus_26g", "Hanada_neg_5g"]
# Identify the rest for ssGSEA
ssgsea_sets = {k: v for k, v in sigs_CD8.items() if k not in sets_ave}
# gseapy.ssgsea expects a DataFrame where rows = genes, columns = samples.
# So let's invert adata back: adata.X => n_cells x n_genes
# We need genes x cells. We'll build a small data frame:

expr_for_ssgsea = pd.DataFrame(
    adata.X.transpose().toarray(),
    index=adata.var.index.tolist(),   # Genes
    columns=adata.obs.index
)
ssgsea_results = gp.ssgsea(data=expr_for_ssgsea,
                               gene_sets=ssgsea_sets,
                               sample_norm=False, # We already normalized above
                               outdir=None,        # don't write to disk
                               min_size=1,
                               max_size=20000)
ssgsea_results.run()
ssgsea_scores = ssgsea_results.res2d # shape = #sets x #cells
# Transpose to cell x set
df_wide = ssgsea_scores.pivot(
    index='Name',    # Each sample’s name/ID
    columns='Term',  # The gene set name
    values='NES'      # Which score you want to spread out as columns
)
# rename columns
df_wide.columns = [f"CD8_{col}" for col in df_wide.columns]
# ----------------------------------------------------------------
# Compute average expression for certain sets
# ----------------------------------------------------------------
# We'll get the log-transformed data from adata.X: (n_cells, n_genes)
# We'll define a small helper to average genes from 'gs'

adata_df = pd.DataFrame(adata.X.toarray(), 
                        index=adata.obs_names, 
                        columns=adata.var.index.tolist())

ave_Hanada_pos_27g = average_genes(sigs_CD8["Hanada_pos_27g"], adata_df)
ave_Hanada_neg_5g  = average_genes(sigs_CD8["Hanada_neg_5g"],  adata_df)
ave_Oliveira_virus_26g = average_genes(sigs_CD8["Oliveira_virus_26g"], adata_df)

# Now combine them with ssgsea_scores
# The original code merges them in a single data.frame
# We'll do the same:
signature_df = pd.DataFrame(index=adata.obs.index)
signature_df["CD8_ave_Hanada_pos_27g"] = ave_Hanada_pos_27g.values
signature_df["CD8_ave_Hanada_neg_5g"]  = ave_Hanada_neg_5g.values
signature_df["CD8_ave_Oliveira_virus_26g"] = ave_Oliveira_virus_26g.values

df_combined = df_wide.join(signature_df, how="inner")
final_df_CD8 = pd.concat([final_df_CD8, df_combined], ignore_index=False)

del(adata)
del(adata_all_CD8)

#%% CD4

adata = adata_all_CD4.copy()
del(adata_all_CD4)
print(adata)
sets_ave_CD4 = ["Hanada_pos_9g", "Hanada_neg_4g"]
# Identify the rest for ssGSEA
ssgsea_sets = {k: v for k, v in sigs_CD4.items() if k not in sets_ave_CD4}
# gseapy.ssgsea expects a DataFrame where rows = genes, columns = samples.
# So let's invert adata back: adata.X => n_cells x n_genes
# We need genes x cells. We'll build a small data frame:

expr_for_ssgsea = pd.DataFrame(
    adata.X.transpose().toarray(),
    index=adata.var.index.tolist(),   # Genes
    columns=adata.obs.index
)
ssgsea_results = gp.ssgsea(data=expr_for_ssgsea,
                               gene_sets=ssgsea_sets,
                               sample_norm=False, # We already normalized above
                               outdir=None,        # don't write to disk
                               min_size=1,
                               max_size=20000)
ssgsea_results.run()
ssgsea_scores = ssgsea_results.res2d # shape = #sets x #cells

# Transpose to cell x set
df_wide = ssgsea_scores.pivot(
    index='Name',    # Each sample’s name/ID
    columns='Term',  # The gene set name
    values='NES'      # Which score you want to spread out as columns
)
# rename columns
df_wide.columns = [f"CD4_{col}" for col in df_wide.columns]
# ----------------------------------------------------------------
# Compute average expression for certain sets
# ----------------------------------------------------------------
# We'll get the log-transformed data from adata.X: (n_cells, n_genes)
# We'll define a small helper to average genes from 'gs'

adata_df = pd.DataFrame(adata.X.toarray(), 
                        index=adata.obs_names, 
                        columns=adata.var.index.tolist())

ave_Hanada_pos_9g = average_genes(sigs_CD4["Hanada_pos_9g"], adata_df)
ave_Hanada_neg_4g = average_genes(sigs_CD4["Hanada_neg_4g"], adata_df)
    
# Now combine them with ssgsea_scores
# The original code merges them in a single data.frame
# We'll do the same:
signature_df = pd.DataFrame(index=adata.obs.index)
signature_df["CD4_ave_Hanada_pos_9g"] = ave_Hanada_pos_9g.values
signature_df["CD4_ave_Hanada_neg_4g"] = ave_Hanada_neg_4g.values

df_combined = df_wide.join(signature_df, how="inner")
final_df_CD4 = pd.concat([final_df_CD4, df_combined], ignore_index=False)

final_df_CD4.to_csv('output_CD4.csv',index = True)
final_df_CD8.to_csv('output_CD8.csv',index = True)