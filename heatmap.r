#画一些gene的heatmap
library(Seurat)
library(stringr)   
library(dplyr)
setwd("/Users/leishengyue/Desktop/data/2024/scRNC-seq/summary/QCsample")


library(Seurat)
data <- readRDS("CTRL.rds")
data1 <- readRDS("MPS1ip6.rds")
data2 <- readRDS("MPS1ip12.rds")
data3 <- readRDS("MPS1ip29.rds")
data4 <- readRDS("MPS1ip52.rds")

data5 <- readRDS("WGDp6.rds")
data6 <- readRDS("WGDp12.rds")
data7 <- readRDS("WGDp29.rds")
data8 <- readRDS("WGDp52.rds")

rna_CTRL <- data[["RNA"]]$`counts.Gene Expression.CTRL`
rna_MPS1ip6 <- data1[["RNA"]]$`counts.Gene Expression.MPS1ip6`
rna_MPS1ip12 <- data2[["RNA"]]$`counts.Gene Expression.MPS1ip12`
rna_MPS1ip29 <- data3[["RNA"]]$`counts.Gene Expression.MPS1ip29`
rna_MPS1ip52 <- data4[["RNA"]]$`counts.Gene Expression.MPS1ip52`

rna_WGDp6 <- data5[["RNA"]]$`counts.Gene Expression.WGDp1`
rna_WGDp12 <- data6[["RNA"]]$`counts.Gene Expression.WGDp7`
rna_WGDp29 <- data7[["RNA"]]$`counts.Gene Expression.WGDp29`
rna_WGDp52 <- data8[["RNA"]]$`counts.Gene Expression.WGDp52`


seurat_CTRL <- CreateSeuratObject(counts = rna_CTRL)
seurat_CTRL$group <- "CTRL"

seurat_MPS1ip6 <- CreateSeuratObject(counts = rna_MPS1ip6)
seurat_MPS1ip6$group <- "MPS1ip6"

seurat_MPS1ip12 <- CreateSeuratObject(counts = rna_MPS1ip12)
seurat_MPS1ip12$group <- "MPS1ip12"

seurat_MPS1ip29 <- CreateSeuratObject(counts = rna_MPS1ip29)
seurat_MPS1ip29$group <- "MPS1ip29"

seurat_MPS1ip52 <- CreateSeuratObject(counts = rna_MPS1ip52)
seurat_MPS1ip52$group <- "MPS1ip52"

seurat_WGDp6 <- CreateSeuratObject(counts = rna_WGDp6)
seurat_WGDp6$group <- "WGDp6"

seurat_WGDp12 <- CreateSeuratObject(counts = rna_WGDp12)
seurat_WGDp12$group <- "WGDp12"

seurat_WGDp29 <- CreateSeuratObject(counts = rna_WGDp29)
seurat_WGDp29$group <- "WGDp29"

seurat_WGDp52 <- CreateSeuratObject(counts = rna_WGDp52)
seurat_WGDp52$group <- "WGDp52"

seurat_combined <- merge(seurat_CTRL, y = list(seurat_MPS1ip6, seurat_MPS1ip12, seurat_MPS1ip29, seurat_MPS1ip52))

seurat_combined <- merge(seurat_CTRL, y = list(seurat_WGDp6, seurat_WGDp12, seurat_WGDp29, seurat_WGDp52))

Idents(seurat_combined) <- seurat_combined$group

## 1. 三组基因
ribosome_related_genes <- c(
  "EEF1D","RPL8","MRPL36","RPL37","MRPL20","MRPS15","RPL22","RPL5",
  "EIF3I","RPF1","EIF4G3","RPL22L1","RPL24","RPL35A","RPL39L","MRPL3",
  "MRPL47","EIF2A","EIF4A2","EIF4G1","MRPL37","RPS8","RPL11","BRIX1",
  "NSA2","TAF9","RPS23","MRPS36","MRPL13","EIF3E","EIF3H","MRTO4","MZT1",
  "POMP","RPL21")

oncogene_genes <- c(
  "PTK2","PVT1","RICTOR","PIK3CB","MECOM","JUN","NOTCH2","RHOC","YBX1",
  "CUL4A","ZEB1","ITGB1","VIM")

CIN70_signature_genes <- c(
  "STAG1","RUVBL1","SMC4","CCT5","CDC20","STMN1","DHCR24","HMGN2",
  "UBE4B","UQCRH","ECT2","CKAP2","HMGB1","RNASEH2B","NRP1","RAD21",
  "PLK2","CCNB1")


## 2. 合并成一个 gene 向量，并确保都在数据里
gene <- unique(c(ribosome_related_genes,
                 oncogene_genes,
                 CIN70_signature_genes))
gene <- intersect(gene, rownames(seurat_combined))
## 3. 计算各 group 的平均表达（用 counts）
gene_cell_exp <- AggregateExpression(
  seurat_combined,
  features = gene,
  group.by = "group",
  slot = "counts"
)$RNA   # matrix

gene_cell_exp <- as.data.frame(gene_cell_exp)

# 列注释信息
df_col <- data.frame(class = colnames(gene_cell_exp))
rownames(df_col) <- colnames(gene_cell_exp)

group_colors <- c(
  'CTRL'    = '#bebebe', 
  'MPS1ip6' = "#ffc5d7",
  'MPS1ip12'= "#ffed86",
  'MPS1ip29'= "#bde6c1",
  'MPS1ip52'= "#94a9ff"
)

library(ComplexHeatmap)
library(circlize)

top_anno <- HeatmapAnnotation(
  Group = df_col$class,
  col   = list(Group = group_colors)
)

# 为每个基因标记所属 group
gene_group <- ifelse(
  rownames(gene_cell_exp) %in% ribosome_related_genes, "Ribosome-related",
  ifelse(
    rownames(gene_cell_exp) %in% oncogene_genes, "Oncogene",
    ifelse(
      rownames(gene_cell_exp) %in% CIN70_signature_genes, "CIN70-signature",
      "Other"
    )
  )
)

# 用于 row_split 的向量（factor 可以控制显示顺序）
gene_branch <- factor(
  gene_group,
  levels = c("Ribosome-related", "Oncogene", "CIN70-signature", "Other")
)

# 2. 按时间顺序排列
gene_cell_exp <- gene_cell_exp[, c("CTRL","MPS1ip6","MPS1ip12","MPS1ip29","MPS1ip52")]

# 3. log1p —— 必须做，否则 trend 不稳定
mat_log <- log1p(as.matrix(gene_cell_exp))

# 4. 行标准化（观看趋势最重要的一步）
marker_exp <- t(scale(t(mat_log), center = TRUE, scale = TRUE))

# 检查 NA / NaN / Inf
sum(is.na(marker_exp))
sum(is.nan(marker_exp))
sum(is.infinite(marker_exp))

# 把问题值换成 0
marker_exp[is.na(marker_exp)]      <- 0
marker_exp[is.nan(marker_exp)]     <- 0
marker_exp[is.infinite(marker_exp)]<- 0

Heatmap(
  marker_exp, 
  cluster_columns = FALSE, 
  show_column_names = TRUE,
  row_split = gene_branch,
  show_row_names = TRUE,
  column_title_gp = gpar(fontsize = 5),
  heatmap_legend_param = list(title = "Expression"),
  col = colorRamp2(c(-1, 0, 1), c("#7c9dff", "white", "#ff9292")),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  top_annotation = top_anno
)


#去掉oncogene和CIN70然后把ribosome按chr arm分类分组
library(Seurat)
library(stringr)   
library(dplyr)
setwd("/Users/leishengyue/Desktop/data/2024/scRNC-seq/summary/QCsample")


library(Seurat)
data <- readRDS("CTRL.rds")
data1 <- readRDS("MPS1ip6.rds")
data2 <- readRDS("MPS1ip12.rds")
data3 <- readRDS("MPS1ip29.rds")
data4 <- readRDS("MPS1ip52.rds")

data5 <- readRDS("WGDp6.rds")
data6 <- readRDS("WGDp12.rds")
data7 <- readRDS("WGDp29.rds")
data8 <- readRDS("WGDp52.rds")

rna_CTRL <- data[["RNA"]]$`counts.Gene Expression.CTRL`
rna_MPS1ip6 <- data1[["RNA"]]$`counts.Gene Expression.MPS1ip6`
rna_MPS1ip12 <- data2[["RNA"]]$`counts.Gene Expression.MPS1ip12`
rna_MPS1ip29 <- data3[["RNA"]]$`counts.Gene Expression.MPS1ip29`
rna_MPS1ip52 <- data4[["RNA"]]$`counts.Gene Expression.MPS1ip52`

rna_WGDp6 <- data5[["RNA"]]$`counts.Gene Expression.WGDp1`
rna_WGDp12 <- data6[["RNA"]]$`counts.Gene Expression.WGDp7`
rna_WGDp29 <- data7[["RNA"]]$`counts.Gene Expression.WGDp29`
rna_WGDp52 <- data8[["RNA"]]$`counts.Gene Expression.WGDp52`


seurat_CTRL <- CreateSeuratObject(counts = rna_CTRL)
seurat_CTRL$group <- "CTRL"

seurat_MPS1ip6 <- CreateSeuratObject(counts = rna_MPS1ip6)
seurat_MPS1ip6$group <- "MPS1ip6"

seurat_MPS1ip12 <- CreateSeuratObject(counts = rna_MPS1ip12)
seurat_MPS1ip12$group <- "MPS1ip12"

seurat_MPS1ip29 <- CreateSeuratObject(counts = rna_MPS1ip29)
seurat_MPS1ip29$group <- "MPS1ip29"

seurat_MPS1ip52 <- CreateSeuratObject(counts = rna_MPS1ip52)
seurat_MPS1ip52$group <- "MPS1ip52"

seurat_WGDp6 <- CreateSeuratObject(counts = rna_WGDp6)
seurat_WGDp6$group <- "WGDp6"

seurat_WGDp12 <- CreateSeuratObject(counts = rna_WGDp12)
seurat_WGDp12$group <- "WGDp12"

seurat_WGDp29 <- CreateSeuratObject(counts = rna_WGDp29)
seurat_WGDp29$group <- "WGDp29"

seurat_WGDp52 <- CreateSeuratObject(counts = rna_WGDp52)
seurat_WGDp52$group <- "WGDp52"

seurat_combined <- merge(seurat_CTRL, y = list(seurat_MPS1ip6, seurat_MPS1ip12, seurat_MPS1ip29, seurat_MPS1ip52))

seurat_combined <- merge(seurat_CTRL, y = list(seurat_WGDp6, seurat_WGDp12, seurat_WGDp29, seurat_WGDp52))

Idents(seurat_combined) <- seurat_combined$group

chr5p <- c("RPL37", "MRPL36", "BRIX1")
chr3q <- c("RPL35A", "RPL24", "RPL22L1", "RPL39L", "EIF4A2", "MRPL3", "EIF4G1", "MRPL47", "EIF2A")
chr8q <- c("RPL8", "RPL30", "EIF3E", "EIF3H", "MRPL13", "EEF1D")
chr1p <- c("RPL11", "RPL5", "EIF4G3", "RPF1", "MRPL37", "MRTO4", "RPS8", "RPL22", "MRPL20", "MRPS15", "EIF3I")
chr13q <- c("RPL21", "POMP", "MZT1")
chr5q <- c("RPS23", "NSA2", "TAF9", "MRPS36")

## 2. 合并成一个 gene 向量，并确保都在数据里
gene <- unique(c(chr5p, chr3q, chr8q,
                 chr1p, chr13q, chr5q))
gene <- intersect(gene, rownames(seurat_combined))
## 3. 计算各 group 的平均表达（用 counts）
gene_cell_exp <- AggregateExpression(
  seurat_combined,
  features = gene,
  group.by = "group",
  slot = "counts"
)$RNA   # matrix

gene_cell_exp <- as.data.frame(gene_cell_exp)

# 列注释信息
df_col <- data.frame(class = colnames(gene_cell_exp))
rownames(df_col) <- colnames(gene_cell_exp)

group_colors <- c(
  'CTRL'    = '#bebebe', 
  'MPS1ip6' = "#ffc5d7",
  'MPS1ip12'= "#ffed86",
  'MPS1ip29'= "#bde6c1",
  'MPS1ip52'= "#94a9ff"
)

library(ComplexHeatmap)
library(circlize)

top_anno <- HeatmapAnnotation(
  Group = df_col$class,
  col   = list(Group = group_colors)
)
# 为每个基因标记所属 chr arm
gene_group <- ifelse(
  rownames(gene_cell_exp) %in% chr1p,  "chr1p",
  ifelse(
    rownames(gene_cell_exp) %in% chr3q,  "chr3q",
    ifelse(
      rownames(gene_cell_exp) %in% chr5p,  "chr5p",
      ifelse(
        rownames(gene_cell_exp) %in% chr5q,  "chr5q",
        ifelse(
          rownames(gene_cell_exp) %in% chr8q,  "chr8q",
          ifelse(
            rownames(gene_cell_exp) %in% chr13q, "chr13q",
            "Other"
          )
        )
      )
    )
  )
)

# 用于 row_split 的向量（factor 可以控制显示顺序）
gene_branch <- factor(
  gene_group,
  levels = c("chr5p", "chr5q", "chr3q", "chr8q", "chr13q", "chr1p", "Other")
)

# 2. 按时间顺序排列
gene_cell_exp <- gene_cell_exp[, c("CTRL","MPS1ip6","MPS1ip12","MPS1ip29","MPS1ip52")]
gene_cell_exp <- gene_cell_exp[, c("CTRL","WGDp6","WGDp12","WGDp29","WGDp52")]

# 3. log1p —— 必须做，否则 trend 不稳定
mat_log <- log1p(as.matrix(gene_cell_exp))

# 4. 行标准化（观看趋势最重要的一步）
marker_exp <- t(scale(t(mat_log), center = TRUE, scale = TRUE))

# 检查 NA / NaN / Inf
sum(is.na(marker_exp))
sum(is.nan(marker_exp))
sum(is.infinite(marker_exp))

# 把问题值换成 0
marker_exp[is.na(marker_exp)]      <- 0
marker_exp[is.nan(marker_exp)]     <- 0
marker_exp[is.infinite(marker_exp)]<- 0

#横轴聚类，无法按照factor排的顺序来
Heatmap(
  marker_exp, 
  cluster_columns = FALSE,          # 列不聚类，保持 CTRL→MPS1ip6→... 的顺序
  show_column_names = TRUE,
  row_split = gene_branch,
  show_row_names = TRUE,
  column_title = "Ribosome-related genes",   # ← 整张图的标题
  column_title_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "Expression"),
  col = colorRamp2(c(-1, 0, 1), c("#7c9dff", "white", "#ff9292")),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  top_annotation = top_anno
)

#每个spilt内聚类，这个目前最好
Heatmap(
  marker_exp, 
  cluster_columns = FALSE,
  row_split = gene_branch,
  cluster_rows = TRUE,              # 每个 split 里面做聚类
  cluster_row_slices = FALSE,       # 但不要在 block 之间再聚类重排顺序
  show_column_names = TRUE,
  show_row_names = TRUE,
  column_title = "Ribosome-related genes",
  column_title_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "Expression"),
  col = colorRamp2(c(-2, -1, 0, 1, 2), c("#003f66", "#0073a8", "white", "#ff6666", "#cc3333")),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  top_annotation = top_anno
)

#完全不聚类，按照factor顺序来但是不方便看
Heatmap(
  marker_exp, 
  cluster_columns = FALSE,
  show_column_names = TRUE,
  cluster_rows = FALSE
  row_split = gene_branch,         
  show_row_names = TRUE,
  column_title = "Ribosome-related genes",   # ← 整张图的标题
  column_title_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "Expression"),
  col = colorRamp2(c(-1, 0, 1), c("#7c9dff", "white", "#ff9292")),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  top_annotation = top_anno
