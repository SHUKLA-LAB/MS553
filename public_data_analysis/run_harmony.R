#!/usr/bin/env Rscript

# ======================== #
#   Load Required Libraries
# ======================== #
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

# ======================== #
#   Define CLI Options
# ======================== #
option_list <- list(
  make_option(
    c("-i", "--input"), type = "character", default = "processed_seurat_qc.rds",
    help = "Path to input Seurat object (post-QC) [default: %default]"
  ),
  make_option(
    c("-o", "--output"), type = "character", default = "processed_seurat_qc_harmony.rds",
    help = "Output file name for Harmony-integrated Seurat object [default: %default]"
  ),
  make_option(
    c("-p", "--pdf"), type = "character", default = "qc_clustering_umaps.pdf",
    help = "Output PDF file for UMAP plots [default: %default]"
  ),
  make_option(
    c("--group_by"), type = "character", default = "Sample",
    help = "Metadata column to integrate across with Harmony [default: %default]"
  ),
  make_option(
    c("--dims"), type = "integer", default = 30,
    help = "Number of Harmony/UMAP dimensions [default: %default]"
  ),
  make_option(
    c("--resolution"), type = "double", default = 0.05,
    help = "Resolution for clustering [default: %default]"
  ),
  make_option(
    c("--min_dist"), type = "double", default = 0.3,
    help = "UMAP minimum distance [default: %default]"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ======================== #
#   Load Seurat Object
# ======================== #
message("Loading input Seurat object: ", opt$input)
seurat_obj <- readRDS(opt$input)

# ======================== #
#   Run Harmony Integration
# ======================== #
message("Running Harmony on group: ", opt$group_by)
seurat_obj <- RunHarmony(
  object = seurat_obj,
  group.by.vars = opt$group_by,
  reduction = "pca",
  assay.use = "SCT",
  verbose = TRUE
)

# ======================== #
#   UMAP on Harmony Reduction
# ======================== #
message("Running UMAP with Harmony reduction...")
seurat_obj <- RunUMAP(
  seurat_obj,
  dims = 1:opt$dims,
  min.dist = opt$min_dist,
  reduction = "harmony",
  reduction.name = "umap_RNA"
)

# ======================== #
#   Clustering
# ======================== #
message("Finding neighbors and clustering (resolution = ", opt$resolution, ")...")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:opt$dims, reduction = "harmony")
seurat_obj <- FindClusters(seurat_obj, resolution = opt$resolution)
seurat_obj$seurat_clusters_RNA <- seurat_obj$seurat_clusters

# ======================== #
#   Plotting
# ======================== #
message("Generating UMAP plots...")
p1 <- DimPlot(seurat_obj, group.by = "seurat_clusters_RNA", label = TRUE, reduction = "umap_RNA", repel = TRUE) +
  NoLegend()
p2 <- DimPlot(seurat_obj, group.by = opt$group_by, label = TRUE, reduction = "umap_RNA")
p3 <- FeaturePlot(seurat_obj, features = "nCount_RNA", reduction = "umap_RNA", order = TRUE)

# ======================== #
#   Save UMAP PDF
# ======================== #
message("Saving UMAP plots to: ", opt$pdf)
pdf(opt$pdf, width = 15, height = 6)
p1 | p2 | p3
dev.off()

# ======================== #
#   Save Final Object
# ======================== #
message("Saving Harmony-processed Seurat object to: ", opt$output)
saveRDS(seurat_obj, file = opt$output)

message("Harmony integration and plotting complete.")

