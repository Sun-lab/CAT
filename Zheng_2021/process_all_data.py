import os
import glob
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

# set target read depth for normalization
target_rd = 3000

# Define thresholds for filtering
min_rd = 500
max_rd = 20000


DATA_DIR = "GSE156728"

def read_table(path, **kwargs):
    """Convenience for gzipped TSV with rownames in col 0."""
    return pd.read_csv(path, sep="\t", index_col=0, compression="gzip", **kwargs)
pd.set_option('display.max_columns', None)
meta_fp = os.path.join(DATA_DIR, "GSE156728_metadata.txt.gz")
meta = read_table(meta_fp)
meta.shape
meta.head()

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
l = ['BC','BCL','ESCA','MM','PACA','RC','THCA','UCEC','OV','FTC']

for i in l:
    df_path = glob.glob(os.path.join(DATA_DIR, f"*{i}_10X.CD8.counts.txt.gz"))[0]
    df = read_table(df_path)
    adata = sc.AnnData(df.T)
    # Check whether adata.obs.index is a subset of meta['cellID']
    assert(set(adata.obs.index).issubset(set(meta.index)))
    
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.obs = adata.obs.join(meta)
    
    new_index = [bc + f"-{i}" for bc in adata.obs.index]
    adata.obs.index = new_index
    
    median_depth = np.median(adata.obs['total_counts'].values)
    print(f"Median depth for {i}: {median_depth}")

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

    sc.pp.normalize_total(adata, target_sum=target_rd)
    sc.pp.log1p(adata)  # log transform

    print(adata)
    sets_ave = ["Hanada_pos_27g", "Oliveira_virus_26g", "Hanada_neg_5g"]

    # Identify the rest for ssGSEA
    ssgsea_sets = {k: v for k, v in sigs_CD8.items() if k not in sets_ave}
    # gseapy.ssgsea expects a DataFrame where rows = genes, columns = samples.
    # So let's invert adata back: adata.X => n_cells x n_genes
    # We need genes x cells. We'll build a small data frame:
    
    expr_for_ssgsea = pd.DataFrame(
        adata.X.T,
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
    adata_df = pd.DataFrame(adata.X, 
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

for i in l:
    df_path = glob.glob(os.path.join(DATA_DIR, f"*{i}_10X.CD4.counts.txt.gz"))[0]
    df = read_table(df_path)
    adata = sc.AnnData(df.T)
    print(adata)
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    
    assert(set(adata.obs.index).issubset(set(meta.index)))
    
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.obs = adata.obs.join(meta)
    
    new_index = [bc + f"-{i}" for bc in adata.obs.index]
    adata.obs.index = new_index
    
    median_depth = np.median(adata.obs['total_counts'].values)
    print(f"Median depth for {i}: {median_depth}")

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

    sc.pp.normalize_total(adata, target_sum=target_rd)
    sc.pp.log1p(adata)  # log transform

    
    print(adata)
    sets_ave_CD4 = ["Hanada_pos_9g", "Hanada_neg_4g"]
    # Identify the rest for ssGSEA
    ssgsea_sets = {k: v for k, v in sigs_CD4.items() if k not in sets_ave_CD4}
    # gseapy.ssgsea expects a DataFrame where rows = genes, columns = samples.
    # So let's invert adata back: adata.X => n_cells x n_genes
    # We need genes x cells. We'll build a small data frame:
    
    expr_for_ssgsea = pd.DataFrame(
        adata.X.T,
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
    adata_df = pd.DataFrame(adata.X, 
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
