## R packages for experiment1_tma (R >= 4.3). Run once on the cluster.
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran <- c(
  "Seurat", "SeuratObject",      # >= 5.0 (v5 assays, FOVs)
  "data.table", "R.utils",       # CosMx flat files (R.utils: reading .csv.gz with fread)
  "dplyr", "tidyr", "ggplot2", "patchwork", "ggrepel", "jsonlite",
  "dbscan", "RANN", "sp",        # dearraying, neighbours / registration, point-in-polygon
  "harmony",                     # slide batch correction
  "pheatmap", "Matrix"
)
optional <- c(
  "presto",                      # much faster FindAllMarkers on 19,000 genes (remotes::install_github("immunogenomics/presto"))
  "arrow", "hdf5r"               # only for a second segmentation delivered as cells.parquet / cell_feature_matrix.h5
)
bioc <- c(
  "UCell",                       # rank-based signature scores
  "edgeR", "limma",              # pseudobulk
  "variancePartition",           # Aim 2: variance explained by organ/patient/slide
  "speckle"                      # propeller: differential composition
)
install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE)
message("Optional, install if needed: ", paste(setdiff(optional, rownames(installed.packages())), collapse = ", "))
