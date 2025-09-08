# %%
print("started")

import os

import gseapy as gp
import pandas as pd
import scanpy as sc
import numpy as np
import pickle
from scipy.io import mmread
# %%
target_rd = 3000
min_rd = 500
max_rd = 20000

data_folder = "data_from_RDS"

with open("signatures_CD8.pkl","rb") as f:
    sigs_CD8 = pickle.load(f)
del sigs_CD8['Lowery_neg_99g']
print({k: len(v) for k, v in sigs_CD8.items()})

# Read the matrix
with open(os.path.join(data_folder,"rna_counts.mtx"), "r") as f:
    for _ in range(10):
        print(f.readline().strip())
rna_matrix = mmread(os.path.join(data_folder,"rna_counts.mtx")).T.tocsr()  # Transpose to cells × genes

print("read rna counts matrix")

# %%

# Read gene and cell names
genes = pd.read_csv(os.path.join(data_folder,"genes_rna.txt"), header=None)[0].tolist()
cells = pd.read_csv(os.path.join(data_folder,"cells_rna.txt"), header=None)[0].tolist()

# Create AnnData object
adata = sc.AnnData(rna_matrix)
adata.var_names = genes
adata.obs_names = cells



# Read the ADT file
adt = pd.read_csv(os.path.join(data_folder,"adt_values.csv"), index_col=0)

# Print all antibody/protein tags (row names)
print("Antibody tags found:")
print(adt.index.tolist())

# %%
meta = pd.read_csv("data_from_RDS/cell_metadata.csv", index_col=0)

#%%
print(meta['timepoint'].value_counts())

# %%

adata.var_names_make_unique()
adata.var["gene_name"] = adata.var_names

sc.pp.calculate_qc_metrics(adata, inplace=True)

# Check percentage of cells with counts < min_rd or > max_rd
low_count = (adata.obs['total_counts'] < min_rd).sum()
high_count = (adata.obs['total_counts'] > max_rd).sum()
total_cells = adata.shape[0]
print(f"Percentage of cells with counts < {min_rd}: {100 * low_count / total_cells:.2f}%")
print(f"Percentage of cells with counts > {max_rd}: {100 * high_count / total_cells:.2f}%")

# Filter cells with total_counts < min_rd or > max_rd
adata = adata[(adata.obs['total_counts'] >= min_rd) & (adata.obs['total_counts'] <= max_rd)]

sc.pp.normalize_total(adata, target_sum=target_rd)
sc.pp.log1p(adata)

print(adata)
print("Final min and max of total_counts for adata:", 
    np.min(adata.obs['total_counts'].values), 
    np.max(adata.obs['total_counts'].values))

# %% pasted from example

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

final_df = pd.DataFrame()

sets_ave = ["Hanada_pos_27g", "Oliveira_virus_26g", "Hanada_neg_5g"]
# Identify the rest for ssGSEA
ssgsea_sets = {k: v for k, v in sigs_CD8.items() if k not in sets_ave}
# gseapy.ssgsea expects a DataFrame where rows = genes, columns = samples.
# So let's invert adata back: adata.X => n_cells x n_genes
# We need genes x cells. We'll build a small data frame:

# %%

all_cells = np.array(adata.obs.index)
n_batches = 100
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
final_df = pd.concat([final_df, df_combined], ignore_index=False)

final_df.to_csv('output.csv',index = True)

# %%

for s1 in sigs_CD8.keys():
    print(len(sigs_CD8[s1]), s1)
# %%
