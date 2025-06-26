This folder contains code to process Upadhye_2023 study data.

Download the data from:
```
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221776/suppl"

wget -c "$BASE/GSE221776_sc_NSCLC_PBT_annotation.txt.gz"
wget -c "$BASE/GSE221776_sc_NSCLC_PBT_umi.txt.gz"
wget -c "$BASE/GSE221776_sc_PBT_CD4_annotation.txt.gz"
wget -c "$BASE/GSE221776_sc_PBT_CD4_umi.txt.gz"
wget -c "$BASE/GSE221776_sc_PBT_CD8_annotation.txt.gz"
wget -c "$BASE/GSE221776_sc_PBT_CD8_umi.txt.gz"
```

### The run.sh will run the python file ``process_all_data.py`` for ssGSEA enrichment on each sample and cell.
On Fred Hutch cluster, run
```
sbatch run.sh
```
The output will be written in score csv files for each cell.

### Then, the ``process_output.ipynb`` will perform analysis, clustering, select cancer reactive cells and visualize the data. 

