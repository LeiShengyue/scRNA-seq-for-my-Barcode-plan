#UMAP-combined
# 加载必要的包
library(Seurat)
library(dplyr)
library(ggplot2)
library(infercnv)

# 设置工作目录
setwd(file.path("/Users/scRNC-seq/summary/QCsample"))

# 读取所有数据集
data_CTRL <- readRDS("CTRL.rds")
data_MPS1ip1 <- readRDS("MPS1ip1.rds")
data_MPS1ip6 <- readRDS("MPS1ip6.rds")
data_MPS1ip12 <- readRDS("MPS1ip12.rds")
data_MPS1ip29 <- readRDS("MPS1ip29.rds")
data_MPS1ip52 <- readRDS("MPS1ip52.rds")

data_WGDp6 <- readRDS("WGDp6.rds")
data_WGDp12 <- readRDS("WGDp12.rds")
data_WGDp29 <- readRDS("WGDp29.rds")
data_WGDp52 <- readRDS("WGDp52.rds")

# 提取特定层的数据
counts_CTRL <- data_CTRL[["RNA"]]$`counts.Gene Expression.CTRL`
counts_MPS1ip1 <- data_MPS1ip1[["RNA"]]$`counts.Gene Expression.MPS1ip1`
counts_MPS1ip6 <- data_MPS1ip6[["RNA"]]$`counts.Gene Expression.MPS1ip6`
counts_MPS1ip12 <- data_MPS1ip12[["RNA"]]$`counts.Gene Expression.MPS1ip12`
counts_MPS1ip29 <- data_MPS1ip29[["RNA"]]$`counts.Gene Expression.MPS1ip29`
counts_MPS1ip52 <- data_MPS1ip52[["RNA"]]$`counts.Gene Expression.MPS1ip52`


counts_WGDp6 <- data_WGDp6[["RNA"]]$`counts.Gene Expression.WGDp1`  # 修正键名
counts_WGDp12 <- data_WGDp12[["RNA"]]$`counts.Gene Expression.WGDp7`  # 修正键名和引号
counts_WGDp29 <- data_WGDp29[["RNA"]]$`counts.Gene Expression.WGDp29`  # 修正键名
counts_WGDp52 <- data_WGDp52[["RNA"]]$`counts.Gene Expression.WGDp52`  

# 合并数据矩阵
combined_matrix <- cbind(counts_CTRL, counts_MPS1ip1, counts_MPS1ip6, counts_MPS1ip12, counts_MPS1ip29, counts_MPS1ip52)

# 合并数据矩阵
combined_matrix <- cbind(counts_CTRL, counts_WGDp6, counts_WGDp12, counts_WGDp29, counts_WGDp52)

# 提取前缀
prefixes <- sub("_.+$", "", colnames(combined_matrix))

# 统计前缀数量
prefix_counts <- table(prefixes)

# 打印前缀数量
print(prefix_counts)

# 删除前缀 "CTRL_" 和 "MPS1i_"
colnames(combined_matrix) <- sub("^(CTRL_|MPS1ip1_|MPS1ip6_|MPS1ip12_|MPS1ip29_|MPS1ip52_)", "", colnames(combined_matrix))

colnames(combined_matrix) <- sub("^(CTRL_|WGDp1_|WGDp7_|WGDp29_|WGDp52_)", "", colnames(combined_matrix))

# 查看列名
column_names <- colnames(combined_matrix)
print(column_names)

# 查找重复的列名
duplicated_columns <- duplicated(column_names)

# 打印重复的列名
duplicated_column_names <- column_names[duplicated_columns]
print(duplicated_column_names)

# 删除重复的列
unique_combined_matrix <- combined_matrix[, !duplicated_columns]

#创建Seurat对象
seurat_object <- CreateSeuratObject(counts = unique_combined_matrix, project = "10X_Project")


seurat_object = infercnv::add_to_seurat(infercnv_output_path="/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/20241028-Passage29+52/MPS1i-CNVoutput",
                                     seurat_obj=seurat_object, # optional
                                     top_n=10)

seurat_object = infercnv::add_to_seurat(infercnv_output_path="/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/20241028-Passage29+52/WGD-CNVoutput",
                                     seurat_obj=seurat_object, # optional
                                     top_n=10)                                   

# saveRDS(seurat_object,file="inferCNV_seurat_obj_MPS1ip1-52.rds")
# saveRDS(seurat_object,file="inferCNV_seurat_obj_WGDp6-52.rds")

# library(harmony)

# head(seurat_object@meta.data)

# 初始化 group 列为 NA
seurat_object$group <- NA

# 按关键词进行分组
seurat_object$group[grepl("CTRL", seurat_object$infercnv_subcluster)] <- "CTRL"
seurat_object$group[grepl("MPS1ip1", seurat_object$infercnv_subcluster)] <- "MPS1ip1"
seurat_object$group[grepl("MPS1ip6", seurat_object$infercnv_subcluster)] <- "MPS1ip6"
seurat_object$group[grepl("MPS1ip12", seurat_object$infercnv_subcluster)] <- "MPS1ip12"
seurat_object$group[grepl("MPS1ip29", seurat_object$infercnv_subcluster)] <- "MPS1ip29"
seurat_object$group[grepl("MPS1ip52", seurat_object$infercnv_subcluster)] <- "MPS1ip52"

seurat_object$group[grepl("CTRL", seurat_object$infercnv_subcluster)] <- "CTRL"
seurat_object$group[grepl("WGDp6", seurat_object$infercnv_subcluster)] <- "WGDp6"
seurat_object$group[grepl("WGDp12", seurat_object$infercnv_subcluster)] <- "WGDp12"
seurat_object$group[grepl("WGDp29", seurat_object$infercnv_subcluster)] <- "WGDp29"
seurat_object$group[grepl("WGDp52", seurat_object$infercnv_subcluster)] <- "WGDp52"

# 检查分组结果
table(seurat_object$group)

table(is.na(seurat_object$group))

# 对合并数据进行标准化、选择变量特征和PCA分析
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)

seurat_object <- RunPCA(seurat_object, npcs = 30)


# 运行UMAP前进行聚类分析
seurat_object <- FindNeighbors(seurat_object, dims = 1:20)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

# 运行UMAP进行降维分析
seurat_object <- RunUMAP(seurat_object, dims = 1:20)

DimPlot(seurat_object, reduction = "umap",
        group.by = 'group', label=TRUE, repel = TRUE, 
        cols = c('CTRL'= '#474747', 'MPS1ip1'='#ff9aa6',
                 'MPS1ip6'='#609dff', 'MPS1ip12'='#ffed86', 'MPS1ip29'='#31C53F', 'MPS1ip52'='#1185c8'))

DimPlot(seurat_object, reduction = "umap",
        group.by = 'group', label=TRUE, repel = TRUE, 
        cols = c('CTRL'= '#474747', 
                 'WGDp6'='#609dff', 'WGDp12'='#ffed86', 'WGDp29'='#31C53F', 'WGDp52'='#1185c8'))


chr1 <- FeaturePlot(seurat_object, reduction="umap", features="has_cnv_chr1") + ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))


chr1_pro <- FeaturePlot(seurat_object, reduction="umap", features="proportion_cnv_chr1") + ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))
chr1_pro

top1 <- DimPlot(seurat_object, reduction="umap", group.by="top_loss_2", pt.size=0.5 )

# 初始化一个新列为 "No CNV"
seurat_object$cnv_status <- "No CNV"

# 查找所有包含 "has_cnv_chr" 的列
cnv_columns <- grep("^has_cnv_chr", colnames(seurat_object@meta.data), value = TRUE)

# 标记有 CNV 的细胞
seurat_object$cnv_status[apply(seurat_object@meta.data[, cnv_columns], 1, any)] <- "CNV"

# 检查标记结果
table(seurat_object$cnv_status)

cnv <- DimPlot(seurat_object, reduction = "umap", group.by = "cnv_status", cols = c('CNV'= '#ffa9b2', 'No CNV'='#a09f9f'))
cnv <- DimPlot(seurat_object, reduction = "umap", group.by = "cnv_status", cols = c('CNV'= '#a9c0ff', 'No CNV'='#a09f9f'))
cnv


library(harmony)

# 使用 Harmony 进行批次校正
seurat_harmony <- RunHarmony(seurat_object, group.by.vars = "group")

# 使用批次校正后的数据进行邻近分析和聚类，选择前 20 个主成分
seurat_harmony <- FindNeighbors(seurat_harmony, reduction = "harmony", dims = 1:20)
seurat_harmony <- FindClusters(seurat_harmony, resolution = 0.5)

# 使用 Harmony 校正后的数据进行 UMAP 降维
seurat_harmony <- RunUMAP(seurat_harmony, reduction = "harmony", dims = 1:20)

DimPlot(seurat_harmony, reduction = "umap",
        group.by = 'group', label=TRUE, repel = TRUE, 
        cols = c('CTRL'= '#474747', 'MPS1ip1'='#ff9aa6',
                 'MPS1ip6'='#609dff', 'MPS1ip12'='#ffed86', 'MPS1ip29'='#31C53F', 'MPS1ip52'='#1185c8'))

DimPlot(seurat_harmony, reduction = "umap",
        group.by = 'group', label=TRUE, repel = TRUE, 
        cols = c('CTRL'= '#474747', 
                 'WGDp6'='#609dff', 'WGDp12'='#ffed86', 'WGDp29'='#31C53F', 'WGDp52'='#1185c8'))




#如果有inferCNV的result和seurat_object不匹配的情况
#提取细胞名称
infercnv_cells <- colnames(infercnv_result@expr.data)
seurat_cells <- colnames(seurat_object)

# 再次检查交集
matching_cells <- intersect(seurat_cells, infercnv_cells)

# 输出匹配细胞数量
length(matching_cells)
# 检查未匹配的细胞名称
unmatched_cells <- setdiff(seurat_cells, infercnv_cells)
head(unmatched_cells)
# 修正子集化
seurat_object <- subset(seurat_object, cells = matching_cells)


# 加载必要的包
library(dplyr)

# 统计 barcode 列的频率
barcode_counts <- seurat_object@meta.data %>%
  count(barcode) %>%   # 按 barcode 分组并计数
  arrange(desc(n))     # 按频率降序排列

# 查看出现最多的 barcode
head(barcode_counts)

# 设置要彩色标注的 barcodes
highlight_barcodes <- c("ATGGTGTGTT", "TTTTTGTATC", "c2_0001", "c2_0002")

# 添加新的列以区分标注
seurat_object@meta.data$highlight <- ifelse(
  seurat_object@meta.data$barcode %in% highlight_barcodes,
  seurat_object@meta.data$barcode, # 使用原始 barcode 值
  "Other" # 其他的设为 "Other"
)

# 提取 UMAP 坐标数据
umap_data <- as.data.frame(Embeddings(seurat_object, reduction = "umap"))
colnames(umap_data) <- c("UMAP_1", "UMAP_2") # 确保列名是 UMAP_1 和 UMAP_2
umap_data$highlight <- seurat_object$highlight

# 绘制 UMAP 图
library(ggplot2)

ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2)) +
    # 绘制灰色点（其他点）
    geom_point(data = subset(umap_data, highlight == "Other"),
               aes(color = highlight),
               size = 1, alpha = 0.6) + # 控制透明度和大小
    # 绘制彩色点（加边框）
    geom_point(data = subset(umap_data, highlight != "Other"),
               aes(fill = highlight),
               shape = 21, # 设置形状为带边框的点
               size = 2, # 调整点的大小
               color = "black") + # 边框颜色为黑色
    scale_fill_manual(
        values = c(
            "ATGGTGTGTT" = "#ffbbbb",
            "TTTTTGTATC" = "#9f9ffd",
            "c2_0001" = "#b3ffb3",
            "c2_0002" = "#e8c4ff"
        )
    ) +
    scale_color_manual(
        values = c("Other" = "#9d9d9d") # 灰色点颜色
    ) +
    guides(fill = guide_legend(override.aes = list(size = 5)), # 调整图例中点的大小
           color = guide_legend(override.aes = list(size = 5))) +
    theme_minimal() +
    theme(
        legend.position = "right", # 图例放在右侧
        plot.title = element_text(hjust = 0.5), # 标题居中
        panel.grid.major = element_blank(), # 移除主网格线
        panel.grid.minor = element_blank(), # 移除次网格线
        panel.background = element_blank(), # 设置背景为纯白
        axis.line = element_line(color = "black", size = 0.5), # 添加坐标轴线
        axis.title = element_text(size = 12), # 调整坐标轴标题大小
        axis.text = element_text(size = 10) # 调整坐标轴数字大小
    ) +
    labs(
        title = "UMAP with Highlighted Barcodes",
        fill = "Highlighted Barcodes",
        color = "Other"
    )


# 确保 barcode 是字符型数据
barcode_counts$barcode <- as.character(barcode_counts$barcode)

# 去掉 "NA" 和 "UNK" 的条形码
barcode_counts_filtered <- barcode_counts[!(barcode_counts$barcode %in% c("NA", "UNK")), ]

# 再次检查过滤后的数据
barcode_counts_filtered <- barcode_counts_filtered[!is.na(barcode_counts_filtered$barcode), ]

# 提取前 100 个频率最高的条形码
barcode_counts_top100 <- head(barcode_counts_filtered, 100)
# 绘制柱状图
library(ggplot2)
ggplot(barcode_counts_top100, aes(x = reorder(barcode, -n), y = n)) +
    geom_bar(stat = "identity", fill = "#ffe291") +
    labs(
        x = "Barcode",
        y = "Count",
        title = "Top 100 Frequent Barcodes"
    ) +
    theme_minimal() +
    theme(
        panel.grid.major = element_blank(), # 移除背景网格线
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8), # 横轴文本旋转为垂直，字体更小
        axis.title = element_text(size = 12),
        axis.line = element_line(color = "black"), # 添加轴线
        plot.title = element_text(hjust = 0.5, size = 14) # 标题居中并调整大小
    ) +
    scale_y_continuous(limits = c(0, 10)) # 设置纵轴范围为 0 到 10
