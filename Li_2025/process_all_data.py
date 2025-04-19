import os
import re
import pandas as pd
import numpy as np
import scanpy as sc
import gseapy as gp  # for ssGSEA
import seaborn as sns
import matplotlib.pyplot as plt
import tempfile
import shutil
import pickle
with open("signatures_CD8.pkl", "rb") as f:
    sigs_CD8 = pickle.load(f)

print({k: len(v) for k, v in sigs_CD8.items()})

with open("signatures_CD4.pkl", "rb") as f:
    sigs_CD4 = pickle.load(f)

print({k: len(v) for k, v in sigs_CD4.items()})

data_dir = "GSE288199_raw"

sample_dirs = [d for d in os.listdir(data_dir)
               if os.path.isdir(os.path.join(data_dir, d))]

def read_10x_custom(sample_dir, var_names='gene_symbols'):
    """
    Use sc.read_10x_mtx() on a folder that has nonstandard 10x filenames by
    copying them into a temporary folder (with the standard 10x names).
    """
    
    # Find the files in sample_dir that end with _matrix.mtx.gz, etc.
    filenames = os.listdir(sample_dir)
    mtx_file  = [f for f in filenames if f.endswith("_matrix.mtx.gz")]
    bc_file   = [f for f in filenames if f.endswith("_barcodes.tsv.gz")]
    feat_file = [f for f in filenames if f.endswith("_features.tsv.gz")]

    if not (mtx_file and bc_file and feat_file):
        raise FileNotFoundError(
            f"Could not find matrix/barcodes/features in {sample_dir} with the suffixes:"
            " _matrix.mtx.gz, _barcodes.tsv.gz, _features.tsv.gz"
        )

    mtx_path  = os.path.join(sample_dir, mtx_file[0])
    bc_path   = os.path.join(sample_dir, bc_file[0])
    feat_path = os.path.join(sample_dir, feat_file[0])

    # Create a temp folder
    tmpdir = tempfile.mkdtemp(prefix="10x_temp_")
    try:
        # Copy the three files with standard 10x names
        shutil.copy(mtx_path,  os.path.join(tmpdir, "matrix.mtx.gz"))
        shutil.copy(bc_path,   os.path.join(tmpdir, "barcodes.tsv.gz"))
        shutil.copy(feat_path, os.path.join(tmpdir, "features.tsv.gz"))

        # Now read_10x_mtx will see them as the standard 10x file names
        adata = sc.read_10x_mtx(
            tmpdir, 
            var_names=var_names, 
            make_unique=True
        )
    finally:
        # Clean up the temp directory
        shutil.rmtree(tmpdir)

    return adata
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
            if gene in se2.index:
                matched_genes.append(gene)
        
        # Print match info
        print(f"{s1}: {len(matched_genes)}/{len(sigs_CD8[s1])} genes are found.")
        
        # Store in gs
        gs[s1] = matched_genes
    return

final_df = pd.DataFrame()
for s1 in sample_dirs:
    print(s1)
    adata = read_10x_custom(f'GSE288199_raw/{s1}')

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)  # log transform
    new_index = [bc.rsplit('-', 1)[0] + f"-{s1}" for bc in adata.obs.index]
    adata.obs.index = new_index

    match_genes_in_sets(sigs_CD8, adata.var)
    match_genes_in_sets(sigs_CD4, adata.var)
        
    sets_ave = ["Hanada_pos_27g", "Oliveira_virus_26g", "Hanada_neg_5g"]
    sets_ave_CD4 = ["Hanada_pos_9g", "Hanada_neg_4g"]
    # Identify the rest for ssGSEA
    ssgsea_sets = {k: v for k, v in sigs_CD8.items() if k not in sets_ave}
    ssgsea_sets_CD4 = {k: v for k, v in sigs_CD4.items() if k not in sets_ave_CD4}
    # gseapy.ssgsea expects a DataFrame where rows = genes, columns = samples.
    # So let's invert adata back: adata.X => n_cells x n_genes
    # We need genes x cells. We'll build a small data frame:
    
    expr_for_ssgsea = pd.DataFrame(
        adata.X.transpose().toarray(),
        index=adata.var.index,   # Genes
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
                            columns=adata.var_names)
    
    
    
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
    signature_df_CD4["CD4_Hanada_pos_9g"] = ave_Hanada_pos_9g.values
    signature_df_CD4["CD4_Hanada_neg_4g"] = ave_Hanada_neg_4g.values

    df_combined_CD4 = df_wide_CD4.join(signature_df_CD4, how="inner")
    df_CD4_CD8 = df_combined.join(df_combined_CD4, how="inner")
    final_df = pd.concat([final_df, df_CD4_CD8], ignore_index=False)

final_df.to_csv('output.csv',index = True)