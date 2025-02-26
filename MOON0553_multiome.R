# Clean environment
rm(list = ls(all.names = TRUE))
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation


###############################################################################
######################## Install and/or Load libraries ########################
###############################################################################

# Load required package for Bioconductor installations
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Function to install and load packages from CRAN, Bioconductor, or GitHub
load_and_install <- function(pkg, source = "CRAN") {
  if (!require(pkg, character.only = TRUE)) {
    if (source == "CRAN") {
      install.packages(pkg, repos = "http://cran.us.r-project.org")
    } else if (source == "Bioconductor") {
      BiocManager::install(pkg)
    } else if (source == "GitHub") {
      devtools::install_github(pkg)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

# List of CRAN packages
cran_pkgs <- c(
  "tidyverse", "Matrix", "RCurl", "scales", "cowplot", "ggplot2", "ggrepel",
  "dplyr", "patchwork", "RColorBrewer","scCustomize", "gridExtra"
)

# List of Bioconductor packages
bioc_pkgs <- c(
  "Seurat", "SingleCellExperiment", "Signac", "harmony", "SingleR",
  "AnnotationHub", "ensembldb", "EnsDb.Hsapiens.v86", 
  "BSgenome.Hsapiens.UCSC.hg38", "qlcMatrix", "Herper", "metap", "multtest", "glmGamPoi"
)

# Load or install CRAN packages
lapply(cran_pkgs, load_and_install)

# Load or install Bioconductor packages
lapply(bioc_pkgs, load_and_install, source = "Bioconductor")

###############################################################################
############################## Prepare Input Data #############################
###############################################################################

sample_dirs <- c(
  'CLL1/',
  'CLL2/',
  'CLL3/'
)

names(sample_dirs) <- c('CLL1', 'CLL2', 'CLL3')

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# Add 'chr' prefix to seqlevels
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# Map seqlevels to UCSC format and remove NA values
mapped_seqlevels <- mapSeqlevels(seqlevels(annotation), "UCSC")
mapped_seqlevels <- mapped_seqlevels[!is.na(mapped_seqlevels)]  # Remove NAs

# Subset the annotation object to include only valid seqlevels
annotation <- keepSeqlevels(annotation, names(mapped_seqlevels), pruning.mode = "coarse")

# Rename the seqlevels to the UCSC format
annotation <- renameSeqlevels(annotation, mapped_seqlevels)

# load counts matrices for each sample:
counts_list <- lapply(sample_dirs, function(x){Seurat::Read10X_h5(paste0(x, 'filtered_feature_bc_matrix.h5'))})

# load isoform counts matrices for each sample:
fragments_list <- lapply(sample_dirs, function(x){paste0(x, 'atac_fragments.tsv.gz')})

# create individual Seurat objects for each sample
seurat_list <- lapply(1:length(sample_dirs), function(i){
  cur <- Seurat::CreateSeuratObject(counts_list[[i]]$`Gene Expression`, assay = "RNA");
  
  # add a column for the original barcode
  cur$barcode <- colnames(cur)
  
  # add a column indicating the sample name
  cur$Sample <- names(sample_dirs)[[i]]
  
  # add the ATAC assay to the seurat object 
  cur[["ATAC"]] <- Signac::CreateChromatinAssay(counts = counts_list[[i]]$Peaks, 
                                                sep = c(":", "-"),
                                                fragments = fragments_list[[i]],
                                                annotation = annotation)
  cur
})

rm(counts_list, annotation)
# merge samples into one Seurat object
seurat_obj <- Reduce(merge, seurat_list)
rm(seurat_list)

###############################################################################
############################### Quality control ###############################
###############################################################################

# Set ATAC as the default assay for the Seurat object
DefaultAssay(seurat_obj) <- "ATAC"

# Calculate nucleosome signal per cell, a metric for chromatin accessibility
seurat_obj <- NucleosomeSignal(seurat_obj)

# Calculate TSS enrichment score per cell, a quality control metric for ATAC-seq data
seurat_obj <- TSSEnrichment(seurat_obj)

# Generate and save a violin plot to visualize QC metrics 
pdf("qc_metrics_violin.pdf", width = 15, height = 6)
VlnPlot(
  object = seurat_obj,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0,
  group.by = "Sample") & theme(axis.title.x = element_blank())
dev.off() 

# Add log10 transformed ratio of genes per UMI to the metadata (quality control metric for RNA-seq data)
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

# Calculate mitochondrial gene content ratio as a percentage and store it in the metadata
seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100  # Normalize mitoRatio for downstream filtering

# Filter out low-quality cells based on multiple criteria:
seurat_obj <- subset(
  x = seurat_obj,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 6000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 300 &
    nucleosome_signal < 2 &
    TSS.enrichment > 4 &
    log10GenesPerUMI > 0.80 & 
    mitoRatio < 0.25)

###############################################################################
################## Normalization, Clustering and Integration ##################
###############################################################################

# process ATAC
DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 5) %>%
  RunTFIDF() %>%
  RunSVD()

# process RNA
DefaultAssay(seurat_obj) <- 'RNA'

# run normalization, feature selection, scaling, and linear dimensional reduction
seurat_obj <- SCTransform(seurat_obj, vst.flavor = "v2") %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

# run harmony to integrate the samples
seurat_obj <- RunHarmony(
  seurat_obj, 
  group.by.vars = 'Sample',
  reduction = 'pca',
  reduction.save = 'harmony_RNA',
  assay.use = 'SCT',
  lambda=NULL)

# run UMAP
seurat_obj <- RunUMAP(
  seurat_obj,
  dims=1:30,
  min.dist=0.3,
  reduction='harmony_RNA',
  reduction.name = 'umap_RNA')

# run clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction='harmony_RNA')
seurat_obj <- FindClusters(seurat_obj, resolution = 0.05)
seurat_obj$seurat_clusters_RNA <- seurat_obj$seurat_clusters

# plot clusters and samples on the UMAP
p1 <- DimPlot(seurat_obj, group.by = 'seurat_clusters_RNA', label=TRUE, reduction='umap_RNA', repel = T)  + 
  NoLegend()
p2 <- DimPlot(seurat_obj, group.by = 'Sample', label=TRUE, reduction='umap_RNA') 
p3 <- FeaturePlot(seurat_obj, features='nCount_RNA', reduction='umap_RNA', order=TRUE) 

# show the plots
pdf("qc_clustering_umaps.pdf", width = 15, height = 6)
p1 | p2 | p3
dev.off()

# Gene markers
tiff("dotplot_immune_markers.tiff", width = 15, height = 7, units = "in", res = 300)

DotPlot(object = seurat_obj, assay = "SCT", features = c("MS4A1", "PAX5", "CD79A", "CD19", "ROR1", "CD40",
                                                         "CD5", "CD3D", "CD3E", "CD3G", "CD28", "IL7R", "CD8A", "CD7",
                                                         "CD247", "NCR1", "FCGR3A", "NCAM1",
                                                         "CD86", "ITGAM", "ITGAX", "CD4", "LYZ", "CD68", "CD14", "CD16")) +
  theme(axis.title.x = element_blank(),               
        axis.text.x = element_text(angle = 90,      
                                   vjust = 0.5, 
                                   hjust = 1)) + 
  labs(y = "Cluster")
dev.off()

 # CLL Clusters subset 
CLL_clusters <- c(0, 1)
DefaultAssay(seurat_obj) <- "SCT"
CLL_seurat <- seurat_obj[, seurat_obj$seurat_clusters %in% CLL_clusters]

###############################################################################
################################# Cell typing #################################
###############################################################################

# Annotate cells as either tumor (if they belong to CLL clusters) or normal
seurat_obj$stage <- ifelse(seurat_obj$seurat_clusters_RNA %in% CLL_clusters, "tumor", "normal")

# Extract normalized counts from the SCT (single-cell transformed) assay data
norm_counts <- LayerData(seurat_obj, assay="SCT", layer="data")

# Load reference data from the Human Primary Cell Atlas for cell type annotation
# Subset the reference to include only specific immune-related cell types of interest
ref <- celldex::HumanPrimaryCellAtlasData()
ref <- ref[,grepl('DC|B_cell|Neutrophils|T_cells|Monocyte|Erythroblast|Macrophage|NK_cell|Platelets|Myelocyte', ref$label.main)]

# Perform cell type annotation using SingleR, comparing normalized counts to the reference
ct_annotated <- SingleR(
  test = norm_counts, 
  ref = ref, 
  labels = ref$label.main, 
  de.method = 'wilcox'
)

# Add the pruned SingleR annotations to the Seurat object metadata
seurat_obj <- AddMetaData(seurat_obj, ct_annotated$pruned.labels, col.name = 'SingleR_HCA')

# Set cell identities (Idents) to the newly added SingleR annotations
Idents(seurat_obj) <- 'SingleR_HCA'

# Visualize the annotated clusters using UMAP and save the plot as a PDF
pdf("singleR_cluster_annotations.pdf", width = 8, height = 6)
DimPlot(seurat_obj, label = TRUE, reduction = "umap_RNA")  
dev.off()

# Rename cluster identities
new.cluster.ids <- c(
  "CLL Cluster 1", "CLL Cluster 2", "T cells","Non-Classical Monocytes", 
  "Mature B cells", "Classical Monocytes")
metadata <- seurat_obj@meta.data
metadata <- metadata %>%
  mutate(cell_label = new.cluster.ids[seurat_clusters])
seurat_obj@meta.data <- metadata

# Visualize the renamed clusters with UMAP and save the plot as a PDF
tiff("annotated_clusters_umap.pdf", width = 8, height = 6, units = "in", res = 300)
Idents(seurat_obj) <- "cell_label"
DimPlot(seurat_obj, reduction = "umap_RNA", label = TRUE, pt.size = 0.3) + NoLegend()
dev.off()

# Plot the clusters for each patient in a UMAP
tiff("clustering_umaps_split_by_sample.tiff", width = 15, height = 7, units = "in", res = 600)
seurat_obj$cell_label <- factor(seurat_obj$cell_label, levels = c("CLL Cluster 1", "CLL Cluster 2", "Mature B cells", "T cells", "Classical Monocytes", "Non-Classical Monocytes"))
Idents(seurat_obj) <- 'cell_label'
DimPlot(seurat_obj, 
        group.by = 'cell_label', 
        label = TRUE, 
        reduction = 'umap_RNA', 
        repel = TRUE, 
        split.by = 'Sample', 
        label.size = 5,  
        pt.size = 1.5    
) + 
  ggtitle("") + 
  xlab("UMAP 1") + 
  ylab("UMAP 2") + 
  NoLegend() + 
  theme(text = element_text(size = 18), strip.text.x = element_text(size = 18, face="bold")
        )  
dev.off()

###############################################################################
############################### BCR gene score ################################
###############################################################################

# List of BCR pathway genes
 bcr_genes <- c(
  "AKT3", "BCL10", "BLNK", "BTK", "CARD11", "CD19", "CD22", "CD72", "CD79A", "CD79B", "CD81", "CHP1", "CHP2", "CHUK", "CR2", "DAPP1", "FCGR2B", "FOS", "GRB2", 
  "GSK3B", "HRAS", "IFITM1", "IKBKB", "IKBKG", "JUN", "KRAS", "LILRB3", "LYN", "MALT1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "NFAT5", "NFATC1", "NFATC2", 
  "NFATC3", "NFATC4", "NFKB1", "NFKBIA", "NFKBIB", "NFKBIE", "NRAS", "PIK3AP1", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", 
  "PIK3R5", "PLCG2", "PPP3CA", "PPP3CB", "PPP3CC", "PPP3R1", "PPP3R2", "PRKCB", "PTPN6", "RAC1", "RAC2", "RAF1", "RASGRP3", "RELA", "SOS1", "SOS2", 
  "SYK", "VAV1", "VAV2", "VAV3", "PAX5")

plot_bcr_signaling <- function(seurat_obj, bcr_genes, sample) {
  # Set the default assay to SCT
  DefaultAssay(seurat_obj) <- "SCT"
  
  # Subset the Seurat object for the BCR genes if a sample is provided
  if (sample != "all") {
    seurat_sub <- subset(x = seurat_obj, subset = Sample == sample)
  } else {
    seurat_sub <- seurat_obj
  }
  
  # Subset by BCR genes
  seurat_sub <- subset(seurat_sub, features = bcr_genes)
  
  Idents(seurat_sub) <- "cell_label"
  
  # Compute the average expression for the BCR genes
  average_expression <- colMeans(seurat_sub@assays$SCT@data, na.rm = TRUE)
  
  # Add the BCR signaling pathway score to the metadata
  seurat_sub@meta.data$BCR.signaling.pathway <- average_expression
  
  # Generate the FeaturePlot
  plot <- FeaturePlot(object = seurat_sub, features = "BCR.signaling.pathway", reduction = "umap_RNA", label = T, repel = T, pt.size = 0.7, label.size = 5) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  return(plot)
}

plot_bcr_density <- function(seurat_obj, bcr_genes) {
  # Set the default assay to SCT
  DefaultAssay(seurat_obj) <- "SCT"
  
  # Subset the Seurat object for the BCR genes
  seurat_sub <- subset(seurat_obj, features = bcr_genes)
  
  # Compute the average expression for the BCR genes
  average_expression <- colMeans(seurat_sub@assays$SCT@data, na.rm = TRUE)
  
  # Add the BCR signaling pathway score to the metadata
  seurat_sub@meta.data$`BCR signaling pathway` <- average_expression
  
  # Generate the FeaturePlot
  plot <- Plot_Density_Custom(seurat_object = seurat_sub, features = "BCR signaling pathway", reduction = "umap_RNA")
  return(plot)
}

# Plot the BCR signaling pathway expression for each sample and stack them in a single plot
bcr_cll1 <- plot_bcr_signaling(seurat_obj, bcr_genes, "CLL1") + ggtitle("CLL1") + xlab("") + ylab("UMAP 2") + theme(legend.position = "right") + theme(legend.position = "none")
bcr_cll2 <- plot_bcr_signaling(seurat_obj, bcr_genes, "CLL2") + ggtitle("CLL2") + xlab("UMAP 1") + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line.y = element_blank()) + NoLegend() 
bcr_cll3 <- plot_bcr_signaling(seurat_obj, bcr_genes, "CLL3")  + ggtitle("CLL3") + xlab("") + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line.y = element_blank()) + theme(legend.position = "right") + labs(color = "BCR score") + theme(legend.title = element_text(face = "bold"))
tiff("BCR_signaling_umaps_split_by_sample.tiff", width = 16, height = 7, units = "in", res = 300)
bcr_umaps <- grid.arrange(bcr_cll1, bcr_cll2, bcr_cll3, ncol=3)
dev.off()

###############################################################################
########################## PRKC family - expression ###########################
###############################################################################

# List of PKC - family genes
prkc_genes <- c("PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCH", "PRKCG", 
                "PRKCI", "PRKD1", "PRKD3", "PRKCQ", "PRKCZ")

# Find DEGs between tumor and normal cells
Idents(seurat_obj) <- "stage"
seurat_obj <- PrepSCTFindMarkers(seurat_obj)
deg_results <- FindMarkers(seurat_obj, ident.1 = "tumor", ident.2 = "normal", logfc.threshold = 0.25, test.use = "wilcox")

# Check the results
head(deg_results)

# Prepare data for volcano plot
deg_results$gene <- rownames(deg_results)
min_non_zero_pval <- min(deg_results$p_val_adj[deg_results$p_val_adj > 0])
deg_results$p_val_adj <- ifelse(deg_results$p_val_adj == 0, min_non_zero_pval, deg_results$p_val_adj)
deg_results$NegLogPValue <- -log10(deg_results$p_val_adj)
deg_results$significant <- ifelse(deg_results$avg_log2FC <= -0.5 & deg_results$p_val_adj < 0.05, "Downregulated",
                        ifelse(deg_results$avg_log2FC >= 0.5 & deg_results$p_val_adj < 0.05, "Upregulated", "Not significant"))
write.table(deg_results, file = "DEG-table_CLL-vs-all_PKC-family.csv", row.names = T, sep = ",", quote = F)

# Volcano plot
create_volcano_plot <- function(deg_results, genes, title) {
  # Create base volcano plot
  max_value <- max(abs(deg_results$avg_log2FC), na.rm = TRUE)
  volcano <- ggplot(deg_results, aes(x = avg_log2FC, y = NegLogPValue)) +
    geom_point(aes(color = significant), alpha = 0.6) +
    scale_color_manual(values = c("#4393C3", "grey44", "#D6604D")) +
    theme_minimal() +
    labs(
      title = title,
      x = expression(paste("log"[2], " fold change")),
      y = expression(paste("-log"[10], "(adjusted p-value)"))
    ) + 
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    guides(color = guide_legend(title = NULL)) +
    scale_x_continuous(limits = c(-max_value, max_value),
                       expand = expansion(mult = c(0.05, 0.05))) +
    theme_classic()
  
  # Add labels for the specified genes
  volcano_labels <- volcano +
    geom_text_repel(
      data = subset(deg_results, gene %in% genes),
      aes(label = paste0("italic(",gene, ")")),
      size = 4,
      parse = TRUE,
      color = "black",
      box.padding = 0.5,  
      point.padding = 0.6,  
      segment.color = "black",  
      segment.size = 0.8,  
      max.overlaps = Inf,  
      min.segment.length = 0,
      position = position_nudge(x = 0.05, y = 0.4) 
    )
  
  # Return the final plot with adjusted theme
  print(volcano_labels + theme(text = element_text(size = 16)))
}

# Create volcano plots for differentially expressed genes (DEG) results, 
pdf("volcano_CLL-vs-all_PKC-family-1.pdf", width = 9, height = 6)
create_volcano_plot(deg_results, prkc_genes, "CLL Clusters 1 & 2 vs All Other Cells")
dev.off()

# Generate a volcano plot specifically for a predefined list of PRKC genes
tiff("volcano_CLL-vs-all_PKC-family-2.tiff", width = 9, height = 6, units = "in", res = 300)
create_volcano_plot(deg_results, c("PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCH", "PRKCG", "PRKCQ", "TCF3", "TCF4", "CDK14"), "CLL Clusters 1 & 2 vs All Other Cells")
dev.off()

# Generate a dot plot for the expression of PRKC genes across cell types
tiff("dotplot_expression_PKC-family.tiff", width = 9, height = 6, units = "in", res = 300)
seurat_obj$cell_label <- factor(seurat_obj$cell_label, levels = rev(c("CLL Cluster 1", "CLL Cluster 2", "Mature B cells", "T cells", "Classical Monocytes", "Non-Classical Monocytes")))
Idents(seurat_obj) <- "cell_label"
DotPlot(seurat_obj, features = prkc_genes, dot.scale = 8, cols = c("#4393C3", "#D6604D")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
dev.off()

# Loop through each gene in the PRKC family and create a UMAP plot for expression per sample
pdf("PRKC-family_expression_per_sample-umap.pdf", height = 5, width = 15)
for (gene in prkc_genes) {
  plot <- FeaturePlot(seurat_obj,features = gene, split.by = "Sample", reduction = 'umap_RNA', 
                      cols = c("azure3", "red1"), max.cutoff = 2, pt.size = 1) + theme(legend.position = "right")
  print(plot) 
}
dev.off() 

# Plot the PRKCB expression for each sample and stack them in a single plot
Idents(seurat_obj) <- "cell_label"

adjusted_colors <- c("#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")

prkcb1 <- FeaturePlot(object = subset(x = seurat_obj, subset = Sample == "CLL1"), features = "PRKCB", reduction = "umap_RNA", label = T, repel = T, pt.size=0.7, label.size = 5) + ggtitle("CLL1") + xlab("") + ylab("UMAP 2") + theme(legend.position = "right")  + guides(color = "none") + scale_colour_gradientn(colours = adjusted_colors)
prkcb2 <- FeaturePlot(object = subset(x = seurat_obj, subset = Sample == "CLL2"), features = "PRKCB", reduction = "umap_RNA", label = T, repel = T, pt.size=0.7, label.size = 5 ) + ggtitle("CLL2") + xlab("UMAP 1") + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line.y = element_blank())  + NoLegend() + scale_colour_gradientn(colours = adjusted_colors)
prkcb3 <- FeaturePlot(object = subset(x = seurat_obj, subset = Sample == "CLL3"), features = "PRKCB", reduction = "umap_RNA", label = T, repel = T, pt.size=0.7, label.size = 5) + ggtitle("CLL3") + xlab("") + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line.y = element_blank()) + theme(legend.position = "right") + scale_colour_gradientn(colours = adjusted_colors) + labs(color = "PRKCB expression") + theme(legend.title = element_text(face = "bold"))
tiff("PRKCB_umaps_split_by_sample.tiff", width = 16, height = 7, units = "in", res = 300)
prkcb_exp <- grid.arrange(prkcb1, prkcb2, prkcb3, ncol=3)
dev.off()

## Heatmaps for BCR pathway genes expression across samples
bcr_genes_sorted <- sort(bcr_genes)

plot_heatmap <- function(seurat_obj, bcr_genes_sorted, sample) {
  sub_df <- subset(x = seurat_obj, subset = Sample == sample)
  Idents(sub_df) <- "cell_label"
  p <- DoHeatmap(sub_df, features = bcr_genes, label = F) + 
    scale_fill_gradientn(colors = c("#4393C3", "black", "#D6604D"))
  p$data$Feature <- factor(p$data$Feature, levels = rev(bcr_genes_sorted))
  tiff(filename=paste0(sample, "_heatmap.tiff"), width = 9, height = 7, units = "in", res = 600)
  print(p)
  dev.off()
}

for (sample in c("CLL1", "CLL2", "CLL3")) {
  plot_heatmap(seurat_obj, bcr_genes_sorted, sample)
}

###############################################################################
######################## BCR pathway - PRKCB density ##########################
###############################################################################

# PRKCB expression umap 
prkcb_exp_umap <- FeaturePlot_scCustom(seurat_obj,
                              features = "PRKCB",
                              reduction = "umap_RNA",
                              label = TRUE, repel = TRUE, label.size = 5) +  
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank())

# Customized BCR signaling UMAP
DefaultAssay(seurat_obj) <- "SCT"
seurat_sub <- subset(seurat_obj, features = bcr_genes)
Idents(seurat_sub) <- "cell_label"
average_expression <- colMeans(seurat_sub@assays$SCT@data, na.rm = TRUE)
seurat_sub@meta.data$`BCR signaling pathway` <- average_expression

bcr_score_umap <- FeaturePlot_scCustom(seurat_object = seurat_sub, 
                                       features = "BCR signaling pathway", 
                                       reduction = "umap_RNA", 
                                       label = TRUE, 
                                       repel = TRUE, 
                                       label.size = 5) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + 
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()) 

# Plot density for PRKCB and BCR signaling pathway
prkcb_density <- Plot_Density_Custom(seurat_object = seurat_obj, features = "PRKCB", reduction = "umap_RNA")
bcr_density <- plot_bcr_density(seurat_obj, bcr_genes)

# Define the y-axis breaks
y_breaks <- c(5, 0, -5, -10, -15)

# Re-create the plots using common fill limits and adjusting the y-axis breaks
common_limits <- c(0, 0.03)

# Modify the prkcb_density plot
prkcb_density <- prkcb_density + 
  scale_fill_gradient(limits = common_limits) + 
  scale_y_continuous(breaks = y_breaks, expand = expansion(mult = c(0.1, 0.1))) + 
  theme(axis.title = element_blank(), 
        axis.text = element_text(size = 12, face = "plain"),  
        axis.line = element_line(size = 0.5),  
        axis.ticks = element_line(size = 0.5),  
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        plot.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) 

# Modify the bcr_density plot
bcr_density <- bcr_density + 
  scale_fill_gradient(limits = common_limits) + 
  scale_y_continuous(breaks = y_breaks, expand = expansion(mult = c(0.1, 0.1))) + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12, face = "plain"),  
        axis.line = element_line(size = 0.5), 
        axis.ticks = element_line(size = 0.5),  
        plot.title = element_blank()) + 
  NoLegend()

# Embed all the plots in a single visualization
bcr_prkcb_density_umaps <- grid.arrange(
  bcr_score_umap, 
  prkcb_exp_umap, 
  bcr_density, 
  prkcb_density, 
  ncol = 2,
  bottom = textGrob("UMAP 1", gp = gpar(fontsize = 14)),
  left = textGrob("UMAP 2", gp = gpar(fontsize = 14), just = "center", rot = 90)
)

###############################################################################
############################# Link peaks to genes #############################
###############################################################################

# Modify the plot_gene_coverage function to use tryCatch
plot_gene_coverage <- function(seurat_obj, gene){
  # Set the default assay to ATAC
  DefaultAssay(seurat_obj) <- "ATAC"
  Idents(seurat_obj) <- "cell_label"
  
  # Compute the GC content for each peak
  seurat_obj <- RegionStats(seurat_obj, genome = BSgenome.Hsapiens.UCSC.hg38)
  
  # Try to link peaks to the specified gene and create the plot
  tryCatch({
    # Link peaks to the specified gene
    seurat_obj <- LinkPeaks(
      object = seurat_obj,
      peak.assay = "ATAC",
      expression.assay = "SCT",
      genes.use = c(gene)
    )
    
    # Generate the coverage plot for the specified gene
    plot <- CoveragePlot(
      object = seurat_obj,
      region = gene,
      features = gene,
      expression.assay = "SCT",
      extend.upstream = 500,
      extend.downstream = 10000
    )
    
    return(plot) # Return the plot if successful
    
  }, error = function(e) {
    message(paste("Skipping gene", gene, "due to error:", e$message))
    return(NULL) # Return NULL if an error occurs
  })
}

# Loop through each gene in the PRKC family and create a coverage plot per sample
tiff("PRKC-family_chromatin_accessibility.tiff", height = 6, width = 11, units = "in", res = 300)
for (gene in c("PRKCB")) {
  plot <- plot_gene_coverage(seurat_obj, gene)
  
  # Only print the plot if it was successfully generated
  if (!is.null(plot)) {
    print(plot)
  }
}
dev.off()

###############################################################################
# Save  object
saveRDS(seurat_obj, file = "CLL1_2_3_multiome_integrated.rds")
