#!/usr/bin/env Rscript

# --- Libraries ---
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
})

# --- Argument parser --- #
option_list <- list(
  make_option(c("--input"), type = "character", help = "Path to preprocessed RDS object")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ===== # ===== Input Setup ===== # ===== #

seurat_obj <- readRDS(opt$input) # "processed_seurat_qc_harmony.rds"

comparisons <- list(
  c("Relapsed: In-house", "Naïve: GSE111014"),
  c("Relapsed: In-house", "Naïve: GSE165087")
)

# ===== # ===== Immune markers ===== # ===== #

tiff("dotplot_immune_markers.tiff", width = 15, height = 7, units = "in", res = 300)
DotPlot(object = seurat_obj, assay = "SCT", features = c("MS4A1", "PAX5", "CD79A", "CD19", "ROR1", "CD40","CD5", "CD3D", "CD3E", "CD3G", "CD28", "IL7R", "CD8A", "CD7","CD247", "NCR1", "FCGR3A", "NCAM1","CD86", "ITGAM", "ITGAX", "CD4", "LYZ", "CD68", "CD14", "CD16")) + theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + labs(y = "Cluster")
dev.off()

# ===== # ===== Score Panels ===== # ===== #

bcr_genes <- c(
  "AKT3", "BCL10", "BLNK", "BTK", "CARD11", "CD19", "CD22", "CD72", "CD79A", "CD79B", "CD81", "CHP1", "CHP2", "CHUK", "CR2", "DAPP1", "FCGR2B", "FOS", "GRB2", 
  "GSK3B", "HRAS", "IFITM1", "IKBKB", "IKBKG", "JUN", "KRAS", "LILRB3", "LYN", "MALT1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "NFAT5", "NFATC1", "NFATC2", 
  "NFATC3", "NFATC4", "NFKB1", "NFKBIA", "NFKBIB", "NFKBIE", "NRAS", "PIK3AP1", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", 
  "PIK3R5", "PLCG2", "PPP3CA", "PPP3CB", "PPP3CC", "PPP3R1", "PPP3R2", "PRKCB", "PTPN6", "RAC1", "RAC2", "RAF1", "RASGRP3", "RELA", "SOS1", "SOS2", 
  "SYK", "VAV1", "VAV2", "VAV3", "PAX5")

Wnt_Catenin <- c(
  "ACVR1", "ACVR1B", "ACVR1C", "ACVR2A", "ACVR2B", "AKT1", "AKT2", "AKT3", "APC", "APC2",
  "APPL1", "APPL2", "AXIN1", "AXIN2", "BCL9", "BMPR2", "BTRC", "CCND1", "CD44", "CDH1",
  "CDH12", "CDH2", "CDH3", "CDH5", "CDKN2A", "CREBBP", "CSNK1A1", "CSNK1D", "CSNK1E", "CSNK1G1",
  "CSNK1G2", "CSNK1G3", "CSNK2A1", "CSNK2A2", "CSNK2B", "CTNNB1", "DKK1", "DKK2", "DKK3", "DKK4",
  "DKKL1", "DVL1", "DVL2", "DVL3", "EP300", "FRAT1", "FRZB", "FZD1", "FZD10", "FZD2",
  "FZD3", "FZD4", "FZD5", "FZD6", "FZD7", "FZD8", "FZD9", "GJA1", "GNAO1", "GNAQ",
  "GSK3A", "GSK3B", "H2BW2", "HDAC1", "HNF1A", "ILK", "JUN", "KREMEN1", "KREMEN2", "LEF1",
  "LRP1", "LRP5", "LRP6", "MAP3K7", "MAP4K1", "MARK2", "MDM2", "MMP7", "MYC", "NLK",
  "NR5A2", "PIN1", "POU5F1", "PPARD", "PPM1J", "PPM1L", "PPP2CA", "PPP2CB", "PPP2R1A", "PPP2R1B",
  "PPP2R2A", "PPP2R2B", "PPP2R2C", "PPP2R2D", "PPP2R3A", "PPP2R3B", "PPP2R5A", "PPP2R5B", "PPP2R5C", "PPP2R5D",
  "PPP2R5E", "PTPA", "RARA", "RARB", "RARG", "RPS27A", "RUVBL2", "SFRP1", "SFRP2", "SFRP4",
  "SFRP5", "SMO", "SOX1", "SOX10", "SOX11", "SOX12", "SOX13", "SOX14", "SOX15", "SOX17",
  "SOX18", "SOX2", "SOX21", "SOX3", "SOX4", "SOX5", "SOX6", "SOX7", "SOX8", "SOX9",
  "SRC", "TAB1", "TCF3", "TCF4", "TCF7", "TCF7L1", "TCF7L2", "TGFB1", "TGFB2", "TGFB3",
  "TGFBR1", "TGFBR2", "TGFBR3", "TLE1", "TLE3", "TLE4", "TP53", "UBA52", "UBB", "UBC",
  "UBD", "WIF1", "WNT1", "WNT10A", "WNT10B", "WNT11", "WNT16", "WNT2", "WNT2B", "WNT3",
  "WNT3A", "WNT4", "WNT5A", "WNT5B", "WNT6", "WNT7A"
)


# ===== # ===== Subset B/CLL Clusters and UMAPs ===== # ===== #

DefaultAssay(seurat_obj) <- "SCT" 

clusters_umap <- DimPlot(seurat_obj, group.by = 'seurat_clusters_RNA', label = F, reduction = "umap_RNA") & ggtitle("Clusters UMAP")

Idents(seurat_obj) <- "seurat_clusters_RNA"

# <- Based on the manual inspection of the immune markers identified above -> #
seurat_obj <- RenameIdents(seurat_obj,
                           "1" = "T/NK",
                           "7" = "T/NK",
                           "10" = "T/NK",
                           "4" = "Monocytes/Macrophages",
                           "11" = "B/CLL + T/NK",
                           "0" = "B/CLL Cells",
                           "2" = "B/CLL Cells",
                           "3" = "B/CLL Cells",
                           "5" = "B/CLL Cells",
                           "6" = "B/CLL Cells",
                           "8" = "B/CLL Cells",
                           "9" = "B/CLL Cells",
                           "12" = "B/CLL Cells"
)

annotated_clusters_umap <- DimPlot(seurat_obj, label = T, reduction = "umap_RNA") & ggtitle("Annotated Clusters UMAP")

tiff("clusters_umap.tiff", width = 11, height = 5, units = "in", res = 300)
clusters_umap | annotated_clusters_umap
dev.off()


CLL_clusters <- c(0, 2, 3, 5, 6, 8, 9, 12) 

# Add 'Cell_Source' column to metadata and UMAP
seurat_obj$Cell_Source <- ifelse(
  seurat_obj$seurat_clusters_RNA %in% CLL_clusters,
  "B/CLL Cells",
  "Others"
)

cell_source_umap <- DimPlot(seurat_obj, group.by = "Cell_Source", label = F, reduction = "umap_RNA") & ggtitle("Cell Source")

# Remove non-B/CLL cells from the Seurat object
seurat_obj <- seurat_obj[, seurat_obj$seurat_clusters %in% CLL_clusters]

# UMAP for B/CLL cells grouped by Dataset
cll_dataset_umaps <- DimPlot(seurat_obj, group.by = "Dataset", label = F, reduction = "umap_RNA")

tiff("public_datasets_umap.tiff", width = 11, height = 5, units = "in", res = 300)
cell_source_umap | cll_dataset_umaps
dev.off()

# Violins for B/CLL UMI counts in each dataset
counts_violin <- ggplot(seurat_obj@meta.data, aes(x = Dataset, y = nCount_RNA, fill=Dataset)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_bw() +
  labs(title = "Total UMI counts per cell by Dataset", y = "UMI counts", x = "")

tiff("violin_counts_per_dataset.tiff", width = 7, height = 4, units = "in", res = 300)
counts_violin
dev.off()

# ---- UMAP PLOT + SCORE CALC ---- #
compute_pathway_score_plot_umap <- function(seurat_obj, gene_set, score_name, umap_file) {
  DefaultAssay(seurat_obj) <- "SCT"
  subset_obj <- subset(seurat_obj, features = gene_set)
  avg_score <- colMeans(subset_obj@assays$SCT@data, na.rm = TRUE)
  subset_obj@meta.data[[score_name]] <- avg_score
  
  # Set color scale limits
  all_scores <- unlist(subset_obj@meta.data[subset_obj$Dataset %in% c("GSE111014", "GSE165087", "inhouse"), score_name])
  score_min <- min(all_scores, na.rm = TRUE)
  score_max <- max(all_scores, na.rm = TRUE)
  
  plot_list <- lapply(c("GSE111014", "GSE165087", "inhouse"), function(ds) {
    title <- switch(ds,
                    "GSE111014" = "Naïve: GSE111014",
                    "GSE165087" = "Naïve: GSE165087",
                    "inhouse" = "Relapsed: In-house"
    )
    
    FeaturePlot(
      object = subset(subset_obj, subset = Dataset == ds),
      reduction = "umap_RNA",
      features = score_name,
      label = FALSE, repel = TRUE,
      pt.size = 0.7, label.size = 5
    ) & scale_colour_gradientn(
      colours = rev(brewer.pal(n = 11, name = "RdBu")),
      limits = c(score_min, score_max)
    ) & ggtitle(title) & xlab("UMAP 1") & ylab("UMAP 2")
  })
  
  # Save side-by-side UMAPs
  tiff(umap_file, width = 15, height = 5, units = "in", res = 300)
  print(patchwork::wrap_plots(plot_list, nrow = 1))
  dev.off()
  
  return(subset_obj)
}


# === Function to plot FeaturePlot UMAP === #
plot_umap <- function(seurat_obj, gene, output_file) {
  # Compute global range of expression
  expr_vals <- FetchData(seurat_obj, vars = gene)[, 1]
  expr_range <- range(expr_vals, na.rm = TRUE)
  
  # Define datasets (adjust based on your object)
  datasets <- unique(seurat_obj$Dataset)
  
  # Generate plots manually per dataset
  plot_list <- lapply(datasets, function(ds) {
    FeaturePlot(
      object = subset(seurat_obj, subset = Dataset == ds),
      reduction = "umap_RNA",
      features = gene,
      pt.size = 0.7, label = FALSE, repel = TRUE
    ) & scale_colour_gradientn(
      colours = rev(brewer.pal(n = 11, name = "RdBu")),
      limits = expr_range
    ) & ggtitle(ds) & xlab("UMAP 1") & ylab("UMAP 2")
  })
  
  # Save output
  tiff(output_file, width = 15, height = 5, units = "in", res = 300)
  print(patchwork::wrap_plots(plot_list, nrow = 1))
  dev.off()
}

# === Function to plot violin with brackets and Wilcoxon tests === #
plot_violin_with_stats <- function(seurat_obj, gene, output_file, y_label) {
  
  df <- FetchData(seurat_obj, vars = c(gene, "Dataset"), slot = "data") %>%
    as_tibble() %>%
    mutate(Dataset = case_when(
      Dataset == "GSE111014" ~ "Naïve: GSE111014",
      Dataset == "GSE165087" ~ "Naïve: GSE165087",
      Dataset == "inhouse" ~ "Relapsed: In-house",
      TRUE ~ Dataset
    )) %>%
    mutate(Dataset = factor(Dataset, levels = c("Naïve: GSE111014", "Naïve: GSE165087", "Relapsed: In-house")))
  
  # Wilcoxon tests
  wilcox_1 <- wilcox.test(as.formula(paste0("`", gene, "` ~ Dataset")), data = df %>% filter(Dataset %in% c("Naïve: GSE111014", "Relapsed: In-house")))
  wilcox_2 <- wilcox.test(as.formula(paste0("`", gene, "` ~ Dataset")), data = df %>% filter(Dataset %in% c("Naïve: GSE165087", "Relapsed: In-house")))
  
  # Stars
  get_stars <- function(pval) {
    if (pval < 0.01) return("***")
    else if (pval < 0.05) return("*")
    else return("ns")
  }
  
  star1 <- get_stars(wilcox_1$p.value)
  star2 <- get_stars(wilcox_2$p.value)
  y_max <- max(df[[gene]], na.rm = TRUE)
  
  # Plot
  p <- ggplot(df, aes(x = Dataset, y = .data[[gene]], fill = Dataset)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    theme_classic() +
    theme(legend.position = "none") +
    ylab(y_label) +
    xlab("") +
    
    # Bracket 1
    geom_segment(aes(x = 1, xend = 3, y = y_max * 1.05, yend = y_max * 1.05)) +
    geom_segment(aes(x = 1, xend = 1, y = y_max * 1.05, yend = y_max * 1.03)) +
    geom_segment(aes(x = 3, xend = 3, y = y_max * 1.05, yend = y_max * 1.03)) +
    annotate("text", x = 2, y = y_max * 1.06, label = star1, size = 5) +
    
    # Bracket 2
    geom_segment(aes(x = 2, xend = 3, y = y_max * 1.13, yend = y_max * 1.13)) +
    geom_segment(aes(x = 2, xend = 2, y = y_max * 1.13, yend = y_max * 1.11)) +
    geom_segment(aes(x = 3, xend = 3, y = y_max * 1.13, yend = y_max * 1.11)) +
    annotate("text", x = 2.5, y = y_max * 1.14, label = star2, size = 5)
  
  tiff(output_file, width = 7, height = 4, units = "in", res = 300)
  print(p)
  dev.off()
}

#########---*-----*-----#########
######## GENERATING PLOTS #######
#########---*-----*-----#########

# --- BCR Pathway --- #
bcr_obj <- compute_pathway_score_plot_umap(
  seurat_obj,
  gene_set = bcr_genes,
  score_name = "BCR.signaling.pathway",
  umap_file = "BCR_signaling_umaps_split_by_sample.tiff"
)

plot_violin_with_stats(
  bcr_obj,
  "BCR.signaling.pathway",
  "BCR_signaling_violins_brackets.tiff",
  "BCR signaling pathway score"
)

# --- Wnt-Catenin Pathway --- #
wnt_obj <- compute_pathway_score_plot_umap(
  seurat_obj,
  gene_set = Wnt_Catenin,
  score_name = "Wnt-Catenin.pathway",
  umap_file = "Wnt-Catenin_umaps_split_by_sample.tiff"
)

plot_violin_with_stats(
  wnt_obj,
  "Wnt-Catenin.pathway",
  "Wnt-Catenin_violins_brackets.tiff",
  "Wnt-Catenin pathway score"
)

# --- PRKCB --- #

plot_umap(seurat_obj, "PRKCB", "PRKCB_split_by_sample.tiff")
plot_violin_with_stats(seurat_obj, "PRKCB", "PRKCB_violins_brackets.tiff", "PRKCB")

# --- housekeeping genes --- #
housekeeping_genes <- c("ACTB", "GAPDH", "RPL13A")


for (gene in housekeeping_genes) {
  plot_umap(seurat_obj, gene, paste0(gene, "_split_by_sample.tiff"))
  plot_violin_with_stats(seurat_obj, gene, paste0(gene, "_violins_brackets.tiff"), gene)
}





