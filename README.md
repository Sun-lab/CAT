# CAT

Cancer Associated T cells

<!--
# Output files of the pipelines:

**cell_meta_data_CD4.csv, cell_meta_data_CD8.csv**: contains scores for each cell, basic cell information, whether a cell is cancer reactive or not, TCR (could have multiple) associated with the cell. These files can be multi-indexed if TCR info is available for each cell.

**cell_meta_data_CD4_cleaned.csv, cell_meta_data_CD4_cleaned.csv**: contains scores for each cell, basic cell information, whether a cell is cancer reactive or not, and one alpha chain and one beta chain from the TCR data (not multi-indexed). If more than 1 alpha/beta chain available for a single cell, we only take the one with highest UMI.

**tcr_meta_data_CD4_cells.csv, tcr_meta_data_CD4_cells.csv**: TCR data for alpha and beta chains including v gene, j gene, cdr3, and cell level information (average score, num_cells have this TCR, num_cancer_reactive_cells have this TCR).

## Run it using local conda environment


```
conda create --name CAT python=3.10.8
conda activate CAT
conda install -c conda-forge pandas numpy scanpy anndata seaborn matplotlib
conda install -c conda-forge python-igraph leidenalg
```

First run ```zz0_get_common_gene_set.ipynb``` to get the common set of genes to be used, and then run the codes within each study-specific folder. Typically to run ```process_all_data.py``` followed by ```process_output.ipynb```. Finally run ```zz1_prepare_data.Rmd``` and ```zz2_construct_training_data.Rmd``` in this folder. 
-->

## NN models for whether a T cell is cancer-associated or not 

- **classify_cancer_reactive_T**
  
NN models were trained for predicting whether a T cell is cancer-associated or not based on gene experssion. For each cell, the output is a numerical score between 0 and 1. One model was trained for CD4 and one model was for CD8.

### Build the conda environment

```
conda env create -f environment.yml
conda activate CAT
```

### Make prediction on test data

The test data should be T cells with gene expression count. The NN model for the corresponding cell type should be used. For example, for CD8 T cells, the NN model trained for CD8 should be used. 

To prepare the test data for making predictions on, these steps are needed:

1. restrict test genes to a set of genes, which are the genes on which the model development data library-size normalize was run. For the genes that do not exist in the test data, sert zero columns.
2. compute per-cell total_counts on the set of genes in point 1. 
3. library-size normalize to target_sum = 3000, then log1p
4. subset to the top-1000 DE genes used as model input, arranged in the same order as in model input

and these steps are all included in an example tutorial for using the CD8 NN model to make predictions on independent test data:

[step5_predict_on_Caushi_data_CD8.py](https://github.com/Sun-lab/CAT/blob/main/classify_cancer_reactive_T/step5_predict_on_Caushi_data_CD8.py)

  

  
