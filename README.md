# CAT

Cancer Associated TCR

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


## To run on the Fred Hutch Cluster, first
```
grabnode
```
and select hardwares (some larger files it is better to use 32 CPUs and 512GB memory). Then enter the working directory:
```
cd /fh/working/sun_w/CAT
```
Load the modules:
```
ml purge 
ml SciPy-bundle/2023.02-gfbf-2022b
ml scenicplus/1.0.0-foss-2022b
ml JupyterLab/4.0.3-GCCcore-12.2.0
jupyter lab --ip=$(hostname) --port=$(fhfreeport) --no-browser
```
Copy the links, open in browser and work from there.

If any additional library is needed (for example scanpy), one can use this command to install into the user directory:
```
pip3 install --user scanpy
```
or search Fred Hutch python modules in terminal:
```
module avail xxx
```


