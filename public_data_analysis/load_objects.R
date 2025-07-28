#!/usr/bin/env Rscript

# ======================== #
#   Load Required Libraries
# ======================== #
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ======================== #
#   Load In-House Dataset
# ======================== #
inhouse_dirs <- setNames(
  c("inhouse/CLL1", "inhouse/CLL2", "inhouse/CLL3"),
  c("CLL1", "CLL2", "CLL3")
)

inhouse_seurats <- lapply(names(inhouse_dirs), function(sample_name) {
  h5_path <- file.path(inhouse_dirs[[sample_name]], "filtered_feature_bc_matrix.h5")
  data <- Read10X_h5(h5_path)
  
  obj <- CreateSeuratObject(counts = data$`Gene Expression`, assay = "RNA")
  obj$Dataset <- "inhouse"
  obj$Sample <- sample_name
  obj$barcode <- colnames(obj)
  obj
})

# ======================== #
#   Load GSE111014 Dataset
# ======================== #
gse111014_data <- Read10X(data.dir = "GSE111014/")
gse111014_obj <- CreateSeuratObject(counts = gse111014_data, project = "GSE111014")
gse111014_obj$Dataset <- "GSE111014"

# Infer sample name from barcode (if possible)
gse111014_obj$Sample <- sub(".*-(scRNA-seq_.*)", "\\1", colnames(gse111014_obj))
if (all(gse111014_obj$Sample == colnames(gse111014_obj))) {
  gse111014_obj$Sample <- "GSE111014"
}
gse111014_obj$barcode <- colnames(gse111014_obj)

# =============================== #
#   Load GSE165087 RNA-Only Data
# =============================== #
gse165087_dirs <- list.dirs("GSE165087_RAW", full.names = TRUE, recursive = FALSE)
names(gse165087_dirs) <- basename(gse165087_dirs)

gse165087_seurats <- lapply(names(gse165087_dirs), function(sample_name) {
  message("Loading GSE165087 sample: ", sample_name)
  sample_dir <- gse165087_dirs[[sample_name]]
  
  data <- tryCatch({
    Read10X(data.dir = sample_dir)
  }, error = function(e) {
    warning(paste("Failed to read", sample_name, ":", e$message))
    return(NULL)
  })
  
  if (is.null(data)) return(NULL)
  
  # If multimodal, extract only the RNA assay
  if (is.list(data)) {
    if (!"Gene Expression" %in% names(data)) {
      warning(paste("Sample", sample_name, "has no 'Gene Expression' matrix. Skipping."))
      return(NULL)
    }
    data <- data[["Gene Expression"]]
  }
  
  if (!inherits(data, "dgCMatrix") || is.null(colnames(data))) {
    warning(paste("Sample", sample_name, "contains invalid RNA matrix. Skipping."))
    return(NULL)
  }
  
  obj <- CreateSeuratObject(counts = data, assay = "RNA")
  obj$Dataset <- "GSE165087"
  obj$Sample <- sample_name
  obj$barcode <- colnames(obj)
  obj
})

# Remove failed samples
gse165087_seurats <- Filter(Negate(is.null), gse165087_seurats)

# ======================== #
#   Merge All Datasets
# ======================== #
all_seurats <- c(inhouse_seurats, list(gse111014_obj), gse165087_seurats)

merged_obj <- Reduce(function(x, y) {
  merge(x, y, add.cell.ids = c(x$Sample[1], y$Sample[1]))
}, all_seurats)

# ======================== #
#   Summary and Save Output
# ======================== #
message("Dataset distribution:")
print(table(merged_obj$Dataset))

message("Sample distribution:")
print(table(merged_obj$Sample))

saveRDS(merged_obj, file = "merged_seurat_all_samples.rds")
message("Saved merged Seurat object to 'merged_seurat_all_samples.rds'")

