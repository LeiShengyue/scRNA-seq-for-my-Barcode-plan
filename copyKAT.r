#非常非常非常慢，跑出来非常非常非常大，电脑要炸了
  library(Seurat)
  setwd(file.path("/Users/experiment/scRNC-seq/summary/QCsample"))

#  raw <- Read10X(data.dir = "MPS1ip52-data")
#   raw <- CreateSeuratObject(counts = raw, project = "copycat.test", min.cells = 0, min.features = 0)
#   exp.rawdata <- as.matrix(raw@assays$RNA@counts)

#读取数据文件
data1 <- readRDS("CTRL.rds")
data2 <- readRDS("MPS1ip1.rds")
data3 <- readRDS("MPS1ip6.rds")
data4 <- readRDS("MPS1ip12.rds")
data5 <- readRDS("MPS1ip29.rds")
data6 <- readRDS("MPS1ip52.rds")

# 提取基因表达计数数据
rna_CTRL <- data1[["RNA"]]$`counts.Gene Expression.CTRL`
rna_MPS1ip1 <- data2[["RNA"]]$`counts.Gene Expression.MPS1ip1`
rna_MPS1ip6 <- data3[["RNA"]]$`counts.Gene Expression.MPS1ip6`
rna_MPS1ip12 <- data4[["RNA"]]$`counts.Gene Expression.MPS1ip12`
rna_MPS1ip29 <- data5[["RNA"]]$`counts.Gene Expression.MPS1ip29`
rna_MPS1ip52 <- data6[["RNA"]]$`counts.Gene Expression.MPS1ip52`

# 合并数据矩阵
combined_matrix <- cbind(rna_CTRL, rna_MPS1ip1, rna_MPS1ip6, rna_MPS1ip12, rna_MPS1ip29, rna_MPS1ip52)


# 提取前缀
prefixes <- sub("_.+$", "", colnames(combined_matrix))

# 统计前缀数量
prefix_counts <- table(prefixes)

# 打印前缀数量
print(prefix_counts)

# 删除前缀 "CTRL_" 和 "MPS1i_"
colnames(combined_matrix) <- sub("^(CTRL_|MPS1ip1_|MPS1ip6_|MPS1ip12_|MPS1ip29_|MPS1ip52_)", "", colnames(combined_matrix))

# 查看列名
column_names <- colnames(combined_matrix)
head(column_names)

# 查找重复的列名
duplicated_columns <- duplicated(column_names)

# 打印重复的列名
duplicated_column_names <- column_names[duplicated_columns]
print(duplicated_column_names)

# 删除重复的列
unique_combined_matrix <- combined_matrix[, !duplicated_columns]

# 创建Seurat对象
seurat_object <- CreateSeuratObject(counts = unique_combined_matrix, project = "10X_Project")

# 查看Seurat对象结构
str(seurat_object)

# 获取Seurat表达矩阵数据
seurat_expression_matrix <- GetAssayData(seurat_object, assay = "RNA", layer = "counts")

# 查看数据
head(Matrix::t(seurat_expression_matrix))

library(copykat)
# 使用 copykat 进行 CNV 分析
copykat_results <- copykat(
  rawmat = seurat_expression_matrix,          # 原始基因表达矩阵
  id.type = "S",                 # 基因ID类型："S"表示基因符号（symbol）；可以设置为"Ensembl"等其他ID类型
  ngene.chr = 5,                 # 每条染色体至少包含的基因数量；减少噪音和提高精度
  win.size = 25,                 # 滑动窗口大小；数值越大越平滑，适合长染色体段的拷贝数分析
  KS.cut = 0.1,                  # KS检验显著性水平的截断值；用于检测肿瘤细胞 vs 正常细胞
  sam.name = "test",             # 分析的样本名称，输出文件和图表以此命名
  distance = "euclidean",        # 聚类中使用的距离度量方式；可以为“euclidean”或“correlation”
  norm.cell.names = "",          # 正常细胞的名称；默认为空表示所有细胞会参与分类
  output.seg = FALSE,            # 是否生成片段化CNV的输出文件；FALSE表示不生成
  plot.genes = TRUE,             # 是否在CNV图表中标出基因位置；TRUE表示绘制
  genome = "hg20",               # 所使用的基因组版本；这里为“hg20”，即人类基因组GRCh38
  n.cores = 4                    # 计算时使用的核心数量；可以根据计算机性能调整
)
# 从 copykat 测试结果中提取预测数据
pred.test <- data.frame(copykat_results$prediction)

# 检查数据的分类情况
table(pred.test$copykat.pred)
# 输出结果示例：
# aneuploid     diploid not.defined 
#       348         450         202 

# 提取 aneuploid（非整倍体）和 diploid（二倍体）细胞
pred.test <- pred.test[which(pred.test$copykat.pred %in% c("aneuploid", "diploid")), ]

# 提取 CNV 矩阵数据
CNA.test <- data.frame(copykat_results$CNAmat)

my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

# 染色体编号转化为奇偶性(1和2),用于区分不同的染色体。
chr <- as.numeric(CNA.test$chrom) %% 2 + 1
# 创建一个黑色到灰色的颜色渐变,用于染色体颜色标注。
rbPal1 <- colorRampPalette(c('black', 'grey'))
# 为染色体创建一个颜色向量,黑色和灰色表示不同的染色体。
CHR <- rbPal1(2)[as.numeric(chr)]
# 扩展CHR，生成一个二列矩阵,作为热图的ColSideColors参数,用于在热图的列侧边添加染色体颜色标注。
chr1 <- cbind(CHR, CHR)

# 使用Dark2调色板创建一个颜色渐变，用于预测类别的颜色标注
rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
# 从pred.test中获取拷贝数预测结果
com.preN <- pred.test$copykat.pred
# 将预测结果转化为颜色(对应于不同的预测类别)。
pred <- rbPal5(2)[as.numeric(factor(com.preN))]

cells <- rbind(pred, pred)
col_breaks <- c(seq(-1, -0.4, length = 50), 
                seq(-0.4, -0.2, length = 150),
                seq(-0.2, 0.2, length = 600), 
                seq(0.2, 0.4, length = 150), 
                seq(0.4, 1, length = 50))

heatmap.3(
  t(CNA.test[, 4:ncol(CNA.test)]),
  dendrogram = "r",
  distfun = function(x) parallelDist::parDist(x, threads = 4, method = "euclidean"),
  hclustfun = function(x) hclust(x, method = "ward.D2"),
  ColSideColors = chr1,
  RowSideColors = cells,
  Colv = NA,
  Rowv = TRUE,
  notecol = "black",
  col = my_palette,
  breaks = col_breaks,
  key = TRUE,
  keysize = 1,
  density.info = "none",
  trace = "none",
  cexRow = 0.1,
  cexCol = 0.1,
  cex.main = 1,
  cex.lab = 0.1,
  symm = FALSE,
  symkey = FALSE,
  symbreaks = TRUE,
  cex = 1,
  cex.main = 4,
  margins = c(10, 10)
)

legend(
  "topright",
  paste("pred.", names(table(com.preN)), sep = ""),
  pch = 15,
  col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1],
  cex = 0.6,
  bty = "n"
)

# 将 CopyKAT 的预测结果添加到 Seurat 对象的 meta.data 中
seurat_object$CopyKAT <- copykat_results$prediction$copykat.pred

# 使用 UMAP 降维图展示 CopyKAT 结果分组
DimPlot(seurat_object, group.by = "CopyKAT", reduction = 'umap') + 
  scale_color_manual(values = c("#F8766D", '#02BFC4', "gray"))
