This folder contains code to process Chen_2024 study data.

### The run.sh will run the python file ``process_all_data.py`` for ssGSEA enrichment on each sample and cell.
On Fred Hutch cluster, run
```
sbatch run.sh
```
The output will be written in a score csv file for each cell.

### Then, the ``process_output.ipynb`` will perform analysis, write TCR meta data csv, clustering, select cancer reactive cells and visualize the data. 
