This folder contains code to process Li_2025 study data.

### First, we download files: 
```
wget -r -np -nH --cut-dirs=6 \
  ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288199/suppl/
```

### Then put files to folders by each sample:
```
# 1) cd into the directory that has all these .gz files
cd /path/to/GSE288199_raw/

# 2) For each file, extract “HN01_post” or “HN01_pre” from the filename 
#    and use it to create (or reuse) a subfolder, then move the file there.

for f in GSE288199_*_T_*; do
    # “cut” by underscore, and re-join fields #2 and #3 => e.g. “HN01_post”
    sample=$(echo "$f" | cut -d'_' -f2,3)
    
    # Make a subfolder named by sample if it does not exist
    mkdir -p "$sample"
    
    # Move the file into that folder
    mv "$f" "$sample/"
done
```

### (optional, files already available in this repo) The pickle files are converted from RDS files:
```
################### Take RDS to pickle ###############################
# In R
library(jsonlite)

# Read the original RDS file
sigs_CD8 <- readRDS("sigs_CD4.rds")

# Write out as a JSON file
write_json(sigs_CD8, "signatures_CD4.json", pretty = TRUE)

# In Python
import json
import pickle
# Read the JSON file
with open("signatures_CD4.json", "r") as f:
   sigs_CD8 = json.load(f)

# Write the same data as a pickle
with open("signatures_CD4.pkl", "wb") as f:
   pickle.dump(sigs_CD8, f)
```

### First step, we run ``predict_cell_types.ipynb`` notebook, this will cluster cells and predict CD4/CD8 cell for each cell.
The output is written to ``cell_type.csv``.

### Then the run.sh will run the python file ``process_all_data.py`` for ssGSEA enrichment on each sample and cell.
On Fred Hutch cluster, run
```
sbatch run.sh
```
The output will be written in output_CD4/CD8.csv.

### Based on the cell types, we will run ``clustering.ipynb``, this will generate clusters for CD4/CD8 using the gene signature sets. 
The cluster result will be written to ``CD4_CD8_clusters.csv``.

### Finally, the ``select_cancer_reactive_cells_CD4_CD8.ipynb``will select cancer reactive cells and visualize the data. 
The final output is saved in ``cell_meta_data_CD4/CD8.csv``. CD8 data is usually more reliable.
