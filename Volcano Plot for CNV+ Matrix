# Volcano Plot for CNV+ Matrix
# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Load processed objects
data1 <- readRDS("/path/to/CTRL.rds")
data2 <- readRDS("/path/to/s23-MPS1ip6.rds")
#data_s23 <- readRDS("/path/to/s23-MPS1ip6.rds")
#data_s5 <- readRDS("/path/to/s5-MPS1ip6.rds")
#data_s10 <- readRDS("/path/to/s10-MPS1ip1.rds")
#data_s20 <- readRDS("/path/to/s2-MPS1ip1.rds")

# Extract specific layer data
counts_CTRL <- data1[["RNA"]]$`counts.Gene Expression.CTRL`
#counts_MPS1ip1 <- data2[["RNA"]]$`counts.Gene Expression.MPS1ip1`

# Create Seurat object
seurat_obj1 <- CreateSeuratObject(counts = counts_CTRL)
seurat_obj1$group <- "CTRL"  # Add a group column labeled as CTRL
seurat_obj2 <- CreateSeuratObject(data2)
seurat_obj2$group <- "s23_MPS1ip6"  # Add a group column labeled as MPS1i

# Data integration and normalization
seurat_obj <- merge(seurat_obj1, y = seurat_obj2)
seurat_obj <- NormalizeData(seurat_obj)

# Select variable features
seurat_obj <- FindVariableFeatures(seurat_obj)

# Scaling and PCA analysis
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

# Join data layers
seurat_obj <- JoinLayers(seurat_obj)

# Check cluster assignments
table(Idents(seurat_obj))

# Differential gene expression analysis
Idents(seurat_obj) <- seurat_obj$group
group_markers <- FindMarkers(seurat_obj, ident.1 = "MPS1ip1", ident.2 = "CTRL", min.pct = 0.2)

# Join data layers
seurat_obj <- JoinLayers(seurat_obj)

# Mark significant differential genes
group_markers$logp <- -log10(group_markers$p_val_adj)
group_markers$gene <- rownames(group_markers)
group_markers$color <- ifelse(group_markers$p_val_adj < 0.05 & abs(group_markers$avg_log2FC) > 0.5, 
                              ifelse(group_markers$avg_log2FC < -0.5 & group_markers$p_val_adj < 0.05, "lightblue", "pink"), "lightgrey")

# Select significant genes
data111 <- subset(group_markers, group_markers$p_val_adj < 0.05 & abs(group_markers$avg_log2FC) > 1)

# Plot the volcano plot
ggplot(group_markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = color)) +
    scale_color_manual(values = c("lightblue" = "lightblue", "pink" = "pink", "lightgrey" = "lightgrey")) +
    geom_text_repel(data = data111, aes(label = gene), size = 3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
    xlab("log2 Fold Change") +
    ylab("-log10 P-Value Adj") +
    ggtitle("CTRL vs. MPS1ip1") +
    theme_minimal() +
    coord_cartesian(xlim = c(-4, 4)) +
    theme(plot.title = element_text(size = 25, face = 'bold', vjust = 0.5, hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 10, face = 'bold'),
          legend.position = 'right',
          legend.key.size = unit(0.2, 'cm'),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 10, face = 'bold', vjust = 0.5, hjust = 0.5),
          axis.text.y = element_text(size = 10, face = 'bold', vjust = 0.5, hjust = 0.5),
          axis.title.x = element_text(size = 10, face = 'bold', vjust = 0.5, hjust = 0.5),
          axis.title.y = element_text(size = 10, face = 'bold', vjust = 0.5, hjust = 0.5),
          panel.background = element_rect(fill = 'transparent', colour = 'black'),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.background = element_rect(fill = 'transparent', colour = 'black'))

# Filter genes with increased expression
up_genes <- group_markers %>%
  filter(avg_log2FC > 0.5 & -log10(p_val_adj) > 10)

down_genes <- group_markers %>%
  filter(avg_log2FC < -0.5 & -log10(p_val_adj) > 10)

# Export to text files
#write.table(up_genes, "/path/to/Up_MPS1ip1.txt", row.names = FALSE, quote = FALSE, sep = "\t")
#write.table(down_genes, "/path/to/Down_MPS1ip1.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(up_genes, "/path/to/Up_s23_MPS1ip6.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(down_genes, "/path/to/Down_s23_MPS1ip6.txt", row.names = FALSE, quote = FALSE, sep = "\t")
