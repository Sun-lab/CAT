import os
import pandas as pd
import scanpy as sc

adata_all = sc.read_h5ad("GSE253173_single_cell_DREAM.h5ad")
print(adata_all)

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
            if gene in se2['gene'].tolist():
                matched_genes.append(gene)
        
        # Print match info
        print(f"{s1}: {len(matched_genes)}/{len(sigs_CD8[s1])} genes are found.")
        
        # Store in gs
        gs[s1] = matched_genes
    return

import gseapy as gp
import numpy as np

final_df_CD8 = pd.DataFrame()
final_df_CD4 = pd.DataFrame()
t_cell = {'CD4':[
'CD4 CTL',
'CD4 Naive',
'CD4 Proliferating',
'CD4 TCM',
'CD4 TEM', 'Treg'], 'CD8':[
'CD8 Naive',
'CD8 Proliferating',
'CD8 TCM',
'CD8 TEM']
}

sc.pp.calculate_qc_metrics(adata_all, inplace=True)

adata_all_CD4 = adata_all[adata_all.obs["predicted.celltype.l2"].isin(t_cell['CD4'])].copy()
adata_all_CD8 = adata_all[adata_all.obs["predicted.celltype.l2"].isin(t_cell['CD8'])].copy()

median_depth = np.median(adata_all_CD4.obs['total_counts'].values)
sc.pp.normalize_total(adata_all_CD4, target_sum=median_depth)
sc.pp.log1p(adata_all_CD4)  # log transform

median_depth = np.median(adata_all_CD8.obs['total_counts'].values)
sc.pp.normalize_total(adata_all_CD8, target_sum=median_depth)
sc.pp.log1p(adata_all_CD8)  # log transform

# pull off the barcode prefix (everything before the first “-”)
barcodes = [bc.rsplit('-', 1)[0] for bc in adata_all_CD8.obs.index]

# grab the library name column
libs     = adata_all_CD8.obs['LibraryName']

# build your new names:  "<barcode>-<LibraryName>"
new_index = [barcodes[i] + '-' + libs[i] for i in range(len( adata_all_CD8.obs.index))]

adata_all_CD8.obs.index=new_index

# pull off the barcode prefix (everything before the first “-”)
barcodes = [bc.rsplit('-', 1)[0] for bc in adata_all_CD4.obs.index]

# grab the library name column
libs     = adata_all_CD4.obs['LibraryName']

# build your new names:  "<barcode>-<LibraryName>"
new_index = [barcodes[i] + '-' + libs[i] for i in range(len(adata_all_CD4.obs.index))]

adata_all_CD4.obs.index=new_index


for lib in t_cell['CD8']:
    print(lib)
    adata = adata_all_CD8[adata_all_CD8.obs["predicted.celltype.l2"]==lib]
    print(adata)
    sets_ave = ["Hanada_pos_27g", "Oliveira_virus_26g", "Hanada_neg_5g"]
    # Identify the rest for ssGSEA
    ssgsea_sets = {k: v for k, v in sigs_CD8.items() if k not in sets_ave}
    # gseapy.ssgsea expects a DataFrame where rows = genes, columns = samples.
    # So let's invert adata back: adata.X => n_cells x n_genes
    # We need genes x cells. We'll build a small data frame:
    
    expr_for_ssgsea = pd.DataFrame(
        adata.X.transpose().toarray(),
        index=adata.var['gene'].tolist(),   # Genes
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
    df_wide.columns = [f"CD8{col}" for col in df_wide.columns]
    # ----------------------------------------------------------------
    # Compute average expression for certain sets
    # ----------------------------------------------------------------
    # We'll get the log-transformed data from adata.X: (n_cells, n_genes)
    # We'll define a small helper to average genes from 'gs'
    
    adata_df = pd.DataFrame(adata.X.toarray(), 
                            index=adata.obs_names, 
                            columns=adata.var['gene'].tolist())
    
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
    final_df_CD8 = pd.concat([final_df_CD8, df_combined], ignore_index=False)
for lib in t_cell['CD4']:
    print(lib)
    adata = adata_all_CD4[adata_all_CD4.obs["predicted.celltype.l2"]==lib]
    print(adata)    
    sets_ave_CD4 = ["Hanada_pos_9g", "Hanada_neg_4g"]
    # Identify the rest for ssGSEA
    ssgsea_sets = {k: v for k, v in sigs_CD4.items() if k not in sets_ave_CD4}
    # gseapy.ssgsea expects a DataFrame where rows = genes, columns = samples.
    # So let's invert adata back: adata.X => n_cells x n_genes
    # We need genes x cells. We'll build a small data frame:
    
    expr_for_ssgsea = pd.DataFrame(
        adata.X.transpose().toarray(),
        index=adata.var['gene'].tolist(),   # Genes
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
    df_wide.columns = [f"CD4{col}" for col in df_wide.columns]
    # ----------------------------------------------------------------
    # Compute average expression for certain sets
    # ----------------------------------------------------------------
    # We'll get the log-transformed data from adata.X: (n_cells, n_genes)
    # We'll define a small helper to average genes from 'gs'
    
    adata_df = pd.DataFrame(adata.X.toarray(), 
                            index=adata.obs_names, 
                            columns=adata.var['gene'].tolist())
    
    ave_Hanada_pos_9g = average_genes(sigs_CD4["Hanada_pos_9g"], adata_df)
    ave_Hanada_neg_4g = average_genes(sigs_CD4["Hanada_neg_4g"], adata_df)
        
    # Now combine them with ssgsea_scores
    # The original code merges them in a single data.frame
    # We'll do the same:
    signature_df = pd.DataFrame(index=adata.obs.index)
    signature_df["Hanada_pos_9g"] = ave_Hanada_pos_9g.values
    signature_df["Hanada_neg_4g"] = ave_Hanada_neg_4g.values

    df_combined = df_wide.join(signature_df, how="inner")
    final_df_CD4 = pd.concat([final_df_CD4, df_combined], ignore_index=False)
final_df_CD4.to_csv('output_CD4.csv',index = True)
final_df_CD8.to_csv('output_CD8.csv',index = True)

