import os
import pandas as pd

folder_path = './GSE252432_RAW'
dataframes = {}

# Loop through all .csv.gz files in the folder
for file_name in os.listdir(folder_path):
    if file_name.endswith('.csv.gz'):
        full_path = os.path.join(folder_path, file_name)
        df = pd.read_csv(full_path, compression='gzip')
        dataframes[file_name] = df  # Store in dict with filename as key

# Example: show one of them
dataframes['GSM8001307_KH_48_TCR.csv.gz'].head()

import scanpy as sc
adata_all = sc.read_h5ad("GSE253173_single_cell_DREAM.h5ad")
print(adata_all)

import pickle

# Load the dictionary from the file
with open('tcr_file_map.pkl', 'rb') as f:
    tcr_file_map = pickle.load(f)

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
            if gene in se2['gene'].tolist():
                matched_genes.append(gene)
        
        # Print match info
        print(f"{s1}: {len(matched_genes)}/{len(sigs_CD8[s1])} genes are found.")
        
        # Store in gs
        gs[s1] = matched_genes
    return

import gseapy as gp
final_df = pd.DataFrame()
t_cell_types = [
'CD4 CTL',
'CD4 Naive',
'CD4 Proliferating',
'CD4 TCM',
'CD4 TEM',
'CD8 Naive',
'CD8 Proliferating',
'CD8 TCM',
'CD8 TEM',
'Treg'
]
adata_all = adata_all[adata_all.obs["predicted.celltype.l2"].isin(t_cell_types)]

for lib in adata_all.obs['LibraryName'].unique():
    print(lib)
    if lib in tcr_file_map:
        adata = adata_all[adata_all.obs['LibraryName']==lib]
        print(adata)
        
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)  # log transform
        
        new_index = [bc.rsplit('-', 1)[0] + f"-{lib}" for bc in adata.obs.index]
        adata.obs.index = new_index
    
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
            index=adata.var['gene'].tolist(),   # Genes
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
        
        joined_df = adata.obs.join(df_CD4_CD8, how='inner')
        
        tcr_df = dataframes[tcr_file_map[lib]]
        new_index = [bc.rsplit('-', 1)[0] + f"-{lib}" for bc in tcr_df.barcode]
        tcr_df.index = new_index

        tcr_df['TCR_entry_number'] = tcr_df.groupby(tcr_df.index).cumcount() + 1
        
        # Reset index to allow merging without losing info
        tcr_df_reset = tcr_df.reset_index()
        joined_df_reset = joined_df.reset_index().rename(columns={'index': 'cell_id'})
        
        # Merge using a left join: keep all cells from joined_df, bring in matching TCR entries
        merged = joined_df_reset.merge(
            tcr_df_reset.rename(columns={'index': 'cell_id'}), 
            on='cell_id', 
            how='left'
        )
        # merged[merged['TCR_entry_number'].notna()]
        # Set multi-index if you want to preserve original + TCR number
        
        merged.set_index(['cell_id', 'TCR_entry_number'], inplace=True)
        # merged[merged.index.get_level_values('TCR_entry_number').notna()]
        final_df = pd.concat([final_df, merged], ignore_index=False)

final_df.to_csv('cell_meta_data.csv',index = True)