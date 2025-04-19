# #Dowload files: 

# wget -r -np -nH --cut-dirs=6 \
#   ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288199/suppl/

# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288199/matrix/GSE288199_series_matrix.txt.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288199/soft/GSE288199_family.soft.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288199/miniml/GSE288199_family.xml.tgz


# # Put files to folders each sample
# # 1) cd into the directory that has all these .gz files
# cd /path/to/GSE288199_raw/

# # 2) For each file, extract “HN01_post” or “HN01_pre” from the filename 
# #    and use it to create (or reuse) a subfolder, then move the file there.

# for f in GSE288199_*_T_*; do
#     # “cut” by underscore, and re-join fields #2 and #3 => e.g. “HN01_post”
#     sample=$(echo "$f" | cut -d'_' -f2,3)
    
#     # Make a subfolder named by sample if it does not exist
#     mkdir -p "$sample"
    
#     # Move the file into that folder
#     mv "$f" "$sample/"
# done

################### Take RDS to pickle ###############################
## # In R
## library(jsonlite)

## # Read the original RDS file
## sigs_CD8 <- readRDS("sigs_CD4.rds")

## # Write out as a JSON file
## write_json(sigs_CD8, "signatures_CD4.json", pretty = TRUE)

## # In Python
# import json
# import pickle
# # Read the JSON file
# with open("signatures_CD4.json", "r") as f:
#     sigs_CD8 = json.load(f)

# # Write the same data as a pickle
# with open("signatures_CD4.pkl", "wb") as f:
#     pickle.dump(sigs_CD8, f)