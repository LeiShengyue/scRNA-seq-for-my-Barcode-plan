library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Set working directory
setwd(file.path("/path/experiment/scRNA-seq/summary/QCsample"))

# Load pre-processed objects
data1 <- readRDS("CTRL.rds")
data2 <- readRDS("WGDp12.rds")

# Extract specific layer data
counts_CTRL <- data1[["RNA"]]$`counts.Gene Expression.CTRL`
counts_WGDp12 <- data2[["RNA"]]$`counts.Gene Expression.WGDp12`

# Create Seurat objects
seurat_obj1 <- CreateSeuratObject(counts = counts_CTRL)
seurat_obj1$group <- "CTRL"  # Add group column labeled as CTRL
seurat_obj2 <- CreateSeuratObject(counts = counts_WGDp12)
seurat_obj2$group <- "WGDp12"  # Add group column labeled as WGDp12

# Merge and normalize data
seurat_obj <- merge(seurat_obj1, y = seurat_obj2)
seurat_obj <- NormalizeData(seurat_obj)

# Select variable features
seurat_obj <- FindVariableFeatures(seurat_obj)

# Scale data and run PCA analysis
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

# Join data layers
seurat_obj <- JoinLayers(seurat_obj)

# Confirm cluster assignment
table(Idents(seurat_obj))

# Perform differential gene expression analysis
Idents(seurat_obj) <- seurat_obj$group
group_markers <- FindMarkers(seurat_obj, ident.1 = "WGDp12", ident.2 = "CTRL", logfc.threshold = 0)

head(group_markers)

# Save differential expression results
write.table(group_markers, quote = FALSE, paste0("DEGs_WGDp12_vs_CTRL.txt"), sep = "\t")
