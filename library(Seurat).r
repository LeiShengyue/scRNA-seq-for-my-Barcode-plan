library(Seurat)
library(infercnv)
library(Matrix)

# 读取数据文件
data1 <- readRDS("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/CTRL.rds")
data2 <- readRDS("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/MPS1ip1.rds")
data3 <- readRDS("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/MPS1ip6.rds")
data4 <- readRDS("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/MPS1ip12.rds")
data5 <- readRDS("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/WGDp1.rds")
data6 <- readRDS("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/WGDp7.rds")

# 提取基因表达计数数据
rna_CTRL <- data1[["RNA"]]$`counts.Gene Expression.CTRL`
rna_MPS1ip1 <- data2[["RNA"]]$`counts.Gene Expression.MPS1ip1`
rna_MPS1ip6 <- data3[["RNA"]]$`counts.Gene Expression.MPS1ip6`
rna_MPS1ip12 <- data4[["RNA"]]$`counts.Gene Expression.MPS1ip12`
rna_WGDp1 <- data5[["RNA"]]$`counts.Gene Expression.WGDp1`
rna_WGDp7 <- data6[["RNA"]]$`counts.Gene Expression.WGDp7`

# 确认数据类型
print(class(rna_CTRL))
print(class(rna_MPS1ip1))
print(class(rna_MPS1ip6))
print(class(rna_MPS1ip12))
print(class(rna_WGDp1))
print(class(rna_WGDp7))

# 合并数据矩阵
combined_matrix <- cbind(rna_CTRL, rna_MPS1ip1, rna_MPS1ip6, rna_MPS1ip12, rna_WGDp1, rna_WGDp7)
# 删掉前缀后再看重复的再删
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


# 创建Seurat对象
seurat_object <- CreateSeuratObject(counts = combined_matrix, project = "10X_Project")

# 查看Seurat对象结构
str(seurat_object)

# 获取Seurat表达矩阵数据
seurat_expression_matrix <- GetAssayData(seurat_object, assay = "RNA", layer = "counts")

# 查看数据
head(Matrix::t(seurat_expression_matrix))

infercnv_expression_matrix <- GetAssayData(seurat_object, assay = "RNA", layer = "counts")

# 提取前缀
prefixes <- sub("_.+$", "", colnames(infercnv_expression_matrix))

# 统计前缀数量
prefix_counts <- table(prefixes)

# 打印前缀数量
print(prefix_counts)

# 删除前缀 "CTRL_" 和 "MPS1i_"
colnames(infercnv_expression_matrix) <- sub("^(CTRL_|MPS1ip1_|MPS1ip6_|MPS1ip12_|WGDp1_|WGDp7_)", "", colnames(infercnv_expression_matrix))

# 查看列名
head(colnames(infercnv_expression_matrix))

# 加载注释文件并指定列名
annotations <- read.table("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/anno/mixanno.txt", header = FALSE, sep = "\t", col.names = c("Cell_Name", "Annotation"))

# 找出注释文件中存在但在 infercnv_expression_matrix 中不存在的细胞名称
extra_cells_in_annotations <- setdiff(annotations$Cell_Name, colnames(infercnv_expression_matrix))

# 打印多余的细胞名称
print(extra_cells_in_annotations)

# 移除多余的细胞名称
annotations_filtered <- annotations[!annotations$Cell_Name %in% extra_cells_in_annotations, ]

# 检查过滤后的注释文件
print(paste("Number of cell names in filtered annotations:", nrow(annotations_filtered)))

# 确保所有过滤后的注释文件中的细胞名称都存在于 infercnv_expression_matrix 的列名中
if (!all(annotations_filtered$Cell_Name %in% colnames(infercnv_expression_matrix))) {
    stop("Some cell names in the filtered annotation file are not found in the data matrix.")
}

# 确保没有多余的空格
annotations_filtered$Cell_Name <- trimws(annotations_filtered$Cell_Name)
annotations_filtered$Annotation <- trimws(annotations_filtered$Annotation)

# 将过滤后的注释文件保存为临时文件
write.table(annotations_filtered, file = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/anno/filteredmixanno.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 重複行名の特定
anno_data <- read.table("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/anno/filteredmix.txt", header = FALSE, sep = "\t")
duplicated_rows <- anno_data[duplicated(anno_data$V1), ]
print(duplicated_rows)

# 重複行の削除
unique_anno_data <- anno_data[!duplicated(anno_data$V1), ]
write.table(unique_anno_data, "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/anno/filteredmix_unique.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 创建infercnv对象
infercnv_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = infercnv_expression_matrix,
    annotations_file = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/anno/filteredmixanno_unique.txt",
    delim = "\t",
    gene_order_file = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/hg38_gencode.txt",
    ref_group_names = "CTRL"
)


# 运行infercnv
infercnv_obj = infercnv::run(infercnv_obj, 
                             output_format="pdf",
                             num_threads = 8,
                             cutoff= 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/mix-output",  # dir is auto-created for storing outputs
                             window_length = 201, # 51 is way to small
                             cluster_by_groups=TRUE, # cluster
                             HMM=TRUE, # turn on to auto-run the HMM prediction of CNV levels
                             HMM_transition_prob=1e-6,
                             HMM_report_by = c("subcluster"), #,"cell"),
                             analysis_mode = c('subclusters'), #, 'cells'),
                             denoise=TRUE,
                             sd_amplifier=1.35)  # sets midpoint for logistic

# 应用中值过滤
infercnv_obj_medianfiltered = infercnv::apply_median_filtering(infercnv_obj)

# 绘制CNV图
infercnv::plot_cnv(infercnv_obj_medianfiltered,
                   out_dir="/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/mix-output",
                   output_format="pdf",
                   output_filename='infercnv.median_filtered',
                   x.range="auto",
                   x.center=1,
                   title = "infercnv",
                   color_safe_pal = FALSE)