library(Matrix)

sparse_mat <- readRDS(gzfile("GSE179994_all.Tcell.rawCounts.rds.gz"))

# Transpose: now cells = rows, genes = cols
sparse_mat_t <- t(sparse_mat)

# Save matrix
writeMM(sparse_mat_t, "GSE179994_all.Tcell.rawCounts.mtx")

# Save corresponding names
writeLines(colnames(sparse_mat), "features.txt")  # genes → var_names
writeLines(rownames(sparse_mat), "barcodes.txt")  # cells → obs_names