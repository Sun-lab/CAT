This folder contains code to process Liu_2022 study data.

The data can be downloaded using wget:
```
wget -P . https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_RAW.tar
wget -P . https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_PBMC.bulkTCR.tsv.gz
wget -P . https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_Tcell.metadata.tsv.gz
wget -P . https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_all.Tcell.rawCounts.rds.gz
wget -P . https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_all.scTCR.tsv.gz
```

### The convert_to_csv.R file will write the rds to csv file that we will later process using python.
On Fred Hutch cluster, run this from the terminal:
```
Rscript convert_to_csv.R
```
### The run.sh will run the python file ``process_all_data.py`` for ssGSEA enrichment on each sample and cell.
On Fred Hutch cluster, run
```
sbatch run.sh
```
The output will be written in a score csv file for each cell.

### Then, the ``process_output.ipynb`` will perform analysis, write TCR meta data csv, clustering, select cancer reactive cells and visualize the data. 
