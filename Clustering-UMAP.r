# Load necessary packages
library(Seurat)
library(dplyr)
library(ggplot2)

# Set the working directory
# setwd(file.path("/experiment/scRNC-seq/summary/QCsample"))

# Read all datasets
data_CTRL <- readRDS("CTRL.rds")
# data_MPS1ip1 <- readRDS("MPS1ip1.rds")
# data_MPS1ip6 <- readRDS("MPS1ip6.rds")
# data_MPS1ip12 <- readRDS("MPS1ip12.rds")
data_WGDp6 <- readRDS("WGDp6.rds")
data_WGDp12 <- readRDS("WGDp12.rds")

# Extract specific layers from the data
counts_CTRL <- data_CTRL[["RNA"]]$`counts.Gene Expression.CTRL`
# counts_MPS1ip1 <- data_MPS1ip1[["RNA"]]$`counts.Gene Expression.MPS1ip1`
# counts_MPS1ip6 <- data_MPS1ip6[["RNA"]]$`counts.Gene Expression.MPS1ip6`
# counts_MPS1ip12 <- data_MPS1ip12[["RNA"]]$`counts.Gene Expression.MPS1ip12`
counts_WGDp6 <- data_WGDp6[["RNA"]]$`counts.Gene Expression.WGDp6`  # Corrected key name
counts_WGDp12 <- data_WGDp12[["RNA"]]$`counts.Gene Expression.WGDp12`  # Corrected key name and quotes

# Create Seurat objects
seurat_CTRL <- CreateSeuratObject(counts = counts_CTRL)
seurat_CTRL$group <- "CTRL"

# seurat_MPS1ip1 <- CreateSeuratObject(counts = counts_MPS1ip1)
# seurat_MPS1ip1$group <- "MPS1ip1"

# seurat_MPS1ip6 <- CreateSeuratObject(counts = counts_MPS1ip6)
# seurat_MPS1ip6$group <- "MPS1ip6"

# seurat_MPS1ip12 <- CreateSeuratObject(counts = counts_MPS1ip12)
# seurat_MPS1ip12$group <- "MPS1ip12"

seurat_WGDp6 <- CreateSeuratObject(counts = counts_WGDp6)
seurat_WGDp6$group <- "WGDp6"

seurat_WGDp12 <- CreateSeuratObject(counts = counts_WGDp12)
seurat_WGDp12$group <- "WGDp12"

# Merge all datasets
seurat_combined <- merge(seurat_CTRL, y = list(seurat_WGDp6, seurat_WGDp12))

# Normalize data, select variable features, and perform PCA analysis on the combined data
seurat_combined <- NormalizeData(seurat_combined)
seurat_combined <- FindVariableFeatures(seurat_combined)
seurat_combined <- ScaleData(seurat_combined)
seurat_combined <- RunPCA(seurat_combined)

# Perform clustering analysis before running UMAP
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:10)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5)

# Run UMAP for dimensionality reduction
seurat_combined <- RunUMAP(seurat_combined, dims = 1:10)

# Plot UMAP and color by group information
DimPlot(seurat_combined, reduction = "umap", label = TRUE, repel = TRUE, group.by = "group") + 
  ggtitle("UMAP of All Groups")

# If you want to color by cluster results instead of group, run the following code:
DimPlot(seurat_combined, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("UMAP of All Groups")
