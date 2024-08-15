library(Seurat)
library(dplyr)
library(patchwork)

# Load the data for all samples
MPS1ip6 <- readRDS("/path/to/MPS1ip6.rds")
MPS1ip1 <- readRDS("/path/to/MPS1ip1.rds")
MPS1ip12 <- readRDS("/path/to/MPS1ip12.rds")
WGDp1 <- readRDS("/path/to/WGDp1.rds")
WGDp7 <- readRDS("/path/to/WGDp7.rds")
CTRL <- readRDS("/path/to/CTRL.rds")

# Create Seurat objects for each sample
sdata_MPS1ip6 <- CreateSeuratObject(MPS1ip6, project = "MPS1ip6")
sdata_MPS1ip1 <- CreateSeuratObject(MPS1ip1, project = "MPS1ip1")
sdata_MPS1ip12 <- CreateSeuratObject(MPS1ip12, project = "MPS1ip12")
sdata_WGDp1 <- CreateSeuratObject(WGDp1, project = "WGDp1")
sdata_WGDp7 <- CreateSeuratObject(WGDp7, project = "WGDp7")
sdata_CTRL <- CreateSeuratObject(CTRL, project = "CTRL")

# Add sample information as metadata
sdata_MPS1ip6$sample <- "MPS1ip6"
sdata_MPS1ip1$sample <- "MPS1ip1"
sdata_MPS1ip12$sample <- "MPS1ip12"
sdata_WGDp1$sample <- "WGDp1"
sdata_WGDp7$sample <- "WGDp7"
sdata_CTRL$sample <- "CTRL"

# Merge all samples into one Seurat object
combined_data <- merge(sdata_CTRL, y = c(sdata_MPS1ip6, sdata_MPS1ip1, sdata_MPS1ip12, sdata_WGDp1, sdata_WGDp7))

# Normalize the data
combined_data <- NormalizeData(combined_data)

# Find highly variable features
combined_data <- FindVariableFeatures(combined_data, selection.method = "vst", nfeatures = 2000)

# Scale the data
combined_data <- ScaleData(combined_data)

# Perform PCA
combined_data <- RunPCA(combined_data, features = VariableFeatures(object = combined_data))

# Determine the number of PCs to use by examining the ElbowPlot
ElbowPlot(combined_data)

# Cluster the cells
combined_data <- FindNeighbors(combined_data, dims = 1:20)  # Adjust the number of PCs based on the ElbowPlot
combined_data <- FindClusters(combined_data, resolution = 0.5)

# Run UMAP for visualization
combined_data <- RunUMAP(combined_data, dims = 1:20)

# Plot the UMAP colored by cluster
DimPlot(combined_data, reduction = "umap", group.by = "seurat_clusters") + ggtitle("UMAP by Clusters")

# Plot the UMAP colored by sample
DimPlot(combined_data, reduction = "umap", group.by = "sample") + ggtitle("UMAP by Sample")

# Identify cluster markers
cluster_markers <- FindAllMarkers(combined_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save the markers to a file
write.csv(cluster_markers, "/path/to/cluster_markers.csv")

# Save the Seurat object
saveRDS(combined_data, file = "/path/to/combined_data_seurat.rds")
