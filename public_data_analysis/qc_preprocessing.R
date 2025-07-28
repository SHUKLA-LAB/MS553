#!/usr/bin/env Rscript

# ======================== #
#   Load Required Libraries
# ======================== #
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ======================== #
#   Define Command-Line Options
# ======================== #
option_list <- list(
  make_option(
    c("-i", "--input"), type = "character", default = "merged_seurat_all_samples.rds",
    help = "Path to input merged Seurat RDS file [default: %default]"
  ),
  make_option(
    c("-o", "--output"), type = "character", default = "processed_seurat_qc.rds",
    help = "Output RDS file name for processed Seurat object [default: %default]"
  ),
  make_option(
    c("--min_umi"), type = "integer", default = 300,
    help = "Minimum UMI count threshold [default: %default]"
  ),
  make_option(
    c("--max_umi"), type = "integer", default = 6000,
    help = "Maximum UMI count threshold [default: %default]"
  ),
  make_option(
    c("--min_log10_genes_per_umi"), type = "double", default = 0.80,
    help = "Minimum log10(genes/UMI) threshold [default: %default]"
  ),
  make_option(
    c("--max_mito_ratio"), type = "double", default = 0.25,
    help = "Maximum mitochondrial gene ratio [default: %default]"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ======================== #
#   Load Merged Seurat Object
# ======================== #
message("Reading input Seurat object: ", opt$input)
seurat_obj <- readRDS(opt$input)

# ======================== #
#   Calculate QC Metrics
# ======================== #
message("Calculating QC metrics...")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj$mitoRatio <- seurat_obj$percent.mt / 100
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

# ======================== #
#   Filter Low-Quality Cells
# ======================== #
message("Filtering cells using thresholds:")
message(" - UMI count: [", opt$min_umi, ", ", opt$max_umi, "]")
message(" - log10(genes/UMI) > ", opt$min_log10_genes_per_umi)
message(" - Mitochondrial ratio < ", opt$max_mito_ratio)

seurat_obj <- subset(
  x = seurat_obj,
  subset = 
    nCount_RNA > opt$min_umi &
    nCount_RNA < opt$max_umi &
    log10GenesPerUMI > opt$min_log10_genes_per_umi &
    mitoRatio < opt$max_mito_ratio
)

# ======================== #
#   Normalize & Cluster
# ======================== #
message("Running SCTransform and clustering...")
seurat_obj <- SCTransform(seurat_obj, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

# ======================== #
#   Save Output
# ======================== #
message("Saving processed object to: ", opt$output)
saveRDS(seurat_obj, file = opt$output)

message("Seurat QC and normalization complete.")

