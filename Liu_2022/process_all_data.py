import pandas as pd
metadata_df = pd.read_csv("GSE179994_Tcell.metadata.tsv.gz", sep="\t", compression="gzip")
import scanpy as sc

# Load the sparse matrix
adata = sc.read_mtx("GSE179994_all.Tcell.rawCounts.mtx")

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
final_df = pd.DataFrame()
t_cell_types = [
'CD4',
'CD8'
]
adata_all = adata[adata.obs["celltype"].isin(t_cell_types)]

for lib in adata_all.obs['celltype'].unique():
    adata = adata_all[adata_all.obs['celltype']==lib]
    print(adata)
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)  # log transform

    match_genes_in_sets(sigs_CD8, adata.var)
    match_genes_in_sets(sigs_CD4, adata.var)
        
    sets_ave = ["Hanada_pos_27g", "Oliveira_virus_26g", "Hanada_neg_5g"]
    sets_ave_CD4 = ["Hanada_pos_9g", "Hanada_neg_4g"]
    # Identify the rest for ssGSEA
    ssgsea_sets = {k: v for k, v in sigs_CD8.items() if k not in sets_ave}
    ssgsea_sets_CD4 = {k: v for k, v in sigs_CD4.items() if k not in sets_ave}
    # gseapy.ssgsea expects a DataFrame where rows = genes, columns = samples.
    # So let's invert adata back: adata.X => n_cells x n_genes
    # We need genes x cells. We'll build a small data frame:
    
    expr_for_ssgsea = pd.DataFrame(
        adata.X.transpose().toarray(),
        index=adata.var.index.tolist(),   # Genes
        columns=adata.obs.index
    )
    
    # gseapy ssgsea
    # Each gene set is a list of genes. We'll create the required format for gseapy.
    # gseapy requires: gene_sets = { 'set_name': [gene1, gene2, ... ], ... }
    
    # The result of ssGSEA is a DataFrame with sets on rows, samples on columns.
    # Then we can transpose to get sample x sets. 
    
    # Large gene lists can be somewhat slow for ssGSEA. 
    # Adjust min_size, max_size, if needed to skip too small/large sets:
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
    df_wide.columns = [f"CD8{col}" for col in df_wide.columns]
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
    signature_df["ave_Hanada_pos_27g"] = ave_Hanada_pos_27g.values
    signature_df["ave_Hanada_neg_5g"]  = ave_Hanada_neg_5g.values
    signature_df["ave_Oliveira_virus_26g"] = ave_Oliveira_virus_26g.values

    df_combined = df_wide.join(signature_df, how="inner")

    ssgsea_results_CD4 = gp.ssgsea(data=expr_for_ssgsea,
                               gene_sets=ssgsea_sets_CD4,
                               sample_norm=False, # We already normalized above
                               outdir=None,        # don't write to disk
                               min_size=1,
                               max_size=20000)
    ssgsea_results_CD4.run()
    ssgsea_scores_CD4 = ssgsea_results_CD4.res2d # shape = #sets x #cells
    # Transpose to cell x set
    df_wide_CD4 = ssgsea_scores_CD4.pivot(
        index='Name',    # Each sample’s name/ID
        columns='Term',  # The gene set name
        values='NES'      # Which score you want to spread out as columns
    )
    # rename columns 
    df_wide_CD4.columns = [f"CD4{col}" for col in df_wide_CD4.columns]
    # ----------------------------------------------------------------
    # Compute average expression for certain sets
    # ----------------------------------------------------------------    

    ave_Hanada_pos_9g = average_genes(sigs_CD4["Hanada_pos_9g"], adata_df)
    ave_Hanada_neg_4g = average_genes(sigs_CD4["Hanada_neg_4g"], adata_df)
    
    # Now combine them with ssgsea_scores
    # The original code merges them in a single data.frame
    # We'll do the same:
    signature_df_CD4 = pd.DataFrame(index=adata.obs.index)
    signature_df_CD4["Hanada_pos_9g"] = ave_Hanada_pos_9g.values
    signature_df_CD4["Hanada_neg_4g"] = ave_Hanada_neg_4g.values

    df_combined_CD4 = df_wide_CD4.join(signature_df_CD4, how="inner")
    df_CD4_CD8 = df_combined.join(df_combined_CD4, how="inner")
    
    final_df = pd.concat([final_df, df_CD4_CD8], ignore_index=False)
final_df.to_csv('cell_scores.csv',index = True)
