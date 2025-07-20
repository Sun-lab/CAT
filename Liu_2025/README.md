This folder contains code to process Chen_2024 study data.

### The python file ``process_all_data.py`` calculates cancer reactive score (e.g., mean value or ssGSEA enrichment) for each cell. The output will be written in output_CD4.csv and output_CD8.csv.

### Then, the ``process_output.ipynb`` will perform analysis, write TCR meta data csv, clustering, select cancer reactive cells and visualize the data, and the resutls are saved in ``cell_meta_data_CD4/CD8_cleaned.csv`` and ``tcr_meta_data_CD4_cells.csv``.

