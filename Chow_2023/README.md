This folder contains code to process Zheng_2021 study data at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156728

### The run.sh will run the python file ``process_all_data.py`` for ssGSEA enrichment on each sample and cell.
On Fred Hutch cluster, run
```
sbatch run.sh
```
The output will be written in a csv file with each cell and TCR info (if any) for each cell.

### Then, the ``process_output.ipynb`` will perform analysis, write TCR meta data csv, clustering, select cancer reactive cells. 
