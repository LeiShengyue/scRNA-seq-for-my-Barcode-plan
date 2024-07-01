#CNV+矩阵的火山图制作
#QC data
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

#读取已处理的对象
data1 <- readRDS("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/CTRL.rds")
data2 <- readRDS("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ matrix data/CNV+ subcluster matrix/s23-MPS1ip6.rds")
#data_s23 <- readRDS("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ matrix data/CNV+ subcluster matrix/s23-MPS1ip6.rds")
#data_s5 <- readRDS("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ matrix data/CNV+ subcluster matrix/s5-MPS1ip6.rds")
#data_s10 <- readRDS("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ matrix data/CNV+ subcluster matrix/s10-MPS1ip1.rds")
#data_s20 <- readRDS("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ matrix data/CNV+ subcluster matrix/s2-MPS1ip1.rds")


# 提取特定层的数据
counts_CTRL <- data1[["RNA"]]$`counts.Gene Expression.CTRL`
#counts_MPS1ip1 <- data2[["RNA"]]$`counts.Gene Expression.MPS1ip1`

# 创建Seurat对象
seurat_obj1 <- CreateSeuratObject(counts = counts_CTRL)
seurat_obj1$group <- "CTRL"  # 添加group列并标记为CTRL
seurat_obj2 <- CreateSeuratObject(data2)
seurat_obj2$group <- "s23_MPS1ip6"  # 添加group列并标记为MPS1i

# 数据整合与标准化
seurat_obj <- merge(seurat_obj1, y = seurat_obj2)
seurat_obj <- NormalizeData(seurat_obj)

# 选择变量特征
seurat_obj <- FindVariableFeatures(seurat_obj)

# 数据缩放和PCA分析
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

# 连接数据层
seurat_obj <- JoinLayers(seurat_obj)

# 确认群集分配
table(Idents(seurat_obj))

# 差异基因表达分析
Idents(seurat_obj) <- seurat_obj$group
group_markers <- FindMarkers(seurat_obj, ident.1 = "MPS1ip1", ident.2 = "CTRL", min.pct = 0.2)

# 连接数据层
seurat_obj <- JoinLayers(seurat_obj)

# 标记显著差异基因
group_markers$logp <- -log10(group_markers$p_val_adj)
group_markers$gene <- rownames(group_markers)
group_markers$color <- ifelse(group_markers$p_val_adj < 0.05 & abs(group_markers$avg_log2FC) > 0.5, 
                              ifelse(group_markers$avg_log2FC < -0.5 & group_markers$p_val_adj < 0.05, "lightblue", "pink"), "lightgrey")

# 选择显著基因
data111 <- subset(group_markers, group_markers$p_val_adj < 0.05 & abs(group_markers$avg_log2FC) > 1)

# 绘制火山图
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

# 筛选表达上升的基因
up_genes <- group_markers %>%
  filter(avg_log2FC > 0.5 & -log10(p_val_adj) > 10)

down_genes <- group_markers %>%
  filter(avg_log2FC < -0.5 & -log10(p_val_adj) > 10)

# 导出为文本文件
#write.table(up_genes, "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ DEG/genes/Up_MPS1ip1.txt", row.names = FALSE, quote = FALSE, sep = "\t")
#write.table(down_genes, "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ DEG/genes/Down_MPS1ip1.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(up_genes, "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ subcluster-MPS1i/genes/Up_s23_MPS1ip6.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(down_genes, "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ subcluster-MPS1i/genes/Down_s23_MPS1ip6.txt", row.names = FALSE, quote = FALSE, sep = "\t")


