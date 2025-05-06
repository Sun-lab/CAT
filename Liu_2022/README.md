This folder contains code to process Chin_2025 study data.

### First, ``prepare_and_find_matches.ipynb`` will download the data.
We find matching tcr files with the corresponding LibraryName in the h5ad scRNA data file.

### Then the run.sh will run the python file ``process_all_data.py`` for ssGSEA enrichment on each sample and cell.
On Fred Hutch cluster, run
```
sbatch run.sh
```
The output will be written in a csv file with each cell and TCR info (if any) for each cell.

### Finally, the ``process_output.ipynb`` will perform analysis, write TCR meta data csv, clustering, select cancer reactive cells and visualize the data. 
