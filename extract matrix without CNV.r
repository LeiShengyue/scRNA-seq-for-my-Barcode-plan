#抽选出没有CNV的matrix和CNV+做对比来滤掉drug影响
library(Seurat)
#先把有CNV的数据读进去，提取一下细胞名
data1 <- readRDS("/Users/leishengyue/Desktop/data/2024/scRNC-seq/summary/QCsample/CNV+_MPS1ip6.rds")
data2 <- readRDS("/Users/leishengyue/Desktop/data/2024/scRNC-seq/summary/QCsample/CNV+_MPS1ip12.rds")
data3 <- readRDS("/Users/leishengyue/Desktop/data/2024/scRNC-seq/summary/QCsample/CNV+_WGDp6.rds")
data4 <- readRDS("/Users/leishengyue/Desktop/data/2024/scRNC-seq/summary/QCsample/CNV+_WGDp12.rds")

M6CNV_cells <- colnames(data1)
M12CNV_cells <- colnames(data2)
W6CNV_cells <- colnames(data3)
W12CNV_cells <- colnames(data4)

#小看一眼
head(M6CNV_cells)
head(M12CNV_cells)
head(W6CNV_cells)
head(W12CNV_cells)

#把完整的data读入
data_MPS1ip6 <- readRDS("/Users/leishengyue/Desktop/data/2024/scRNC-seq/summary/QCsample/MPS1ip6.rds")
data_MPS1ip12 <- readRDS("/Users/leishengyue/Desktop/data/2024/scRNC-seq/summary/QCsample/MPS1ip12.rds")
data_WGDp6 <- readRDS("/Users/leishengyue/Desktop/data/2024/scRNC-seq/summary/QCsample/WGDp6.rds")
data_WGDp12 <- readRDS("/Users/leishengyue/Desktop/data/2024/scRNC-seq/summary/QCsample/WGDp12.rds")

#把count的matrix提出
rna_MPS1ip6 <- data_MPS1ip6[["RNA"]]$`counts.Gene Expression.MPS1ip6`
rna_MPS1ip12 <- data_MPS1ip12[["RNA"]]$`counts.Gene Expression.MPS1ip12`
rna_WGDp6 <- data_WGDp6[["RNA"]]$`counts.Gene Expression.WGDp1`
rna_WGDp12 <- data_WGDp12[["RNA"]]$`counts.Gene Expression.WGDp7`

#把colnames里的前缀删一下咯
prefixesM6 <- sub("_.+$", "", colnames(rna_MPS1ip6))
prefixesM12 <- sub("_.+$", "", colnames(rna_MPS1ip12))
prefixesW6 <- sub("_.+$", "", colnames(rna_WGDp6))
prefixesW12 <- sub("_.+$", "", colnames(rna_WGDp12))

pre_counts_M6 <- table(prefixesM6)
pre_counts_M12 <- table(prefixesM12)
pre_counts_W6 <- table(prefixesW6)
pre_counts_W12 <- table(prefixesW12)

print(pre_counts_M6)
print(pre_counts_M12)
print(pre_counts_W6)
print(pre_counts_W12)

# 删除前缀
colnames(rna_MPS1ip6) <- sub("^(MPS1ip6_)", "", colnames(rna_MPS1ip6))
colnames(rna_MPS1ip12) <- sub("^(MPS1ip12_)", "", colnames(rna_MPS1ip12))
colnames(rna_WGDp6) <- sub("^(WGDp1_)", "", colnames(rna_WGDp6))
colnames(rna_WGDp12) <- sub("^(WGDp7_)", "", colnames(rna_WGDp12))

#把没有CNV的部分的矩阵提出来
noCNV_MPS1ip6 <- rna_MPS1ip6[,!(colnames(rna_MPS1ip6) %in% M6CNV_cells)]
noCNV_MPS1ip12 <- rna_MPS1ip12[,!(colnames(rna_MPS1ip12) %in% M12CNV_cells)]
noCNV_WGDp6 <- rna_WGDp6[,!(colnames(rna_WGDp6) %in% W6CNV_cells)]
noCNV_WGDp12 <- rna_WGDp12[,!(colnames(rna_WGDp12) %in% W12CNV_cells)]

#本人决定做一下数量对比看一眼有多少的量哈
dim(noCNV_MPS1ip6)
dim(data1)
dim(rna_MPS1ip6)

dim(noCNV_MPS1ip12)
dim(data2)
dim(rna_MPS1ip12)

dim(noCNV_WGDp6)
dim(data3)
dim(rna_WGDp6)

dim(noCNV_WGDp12)
dim(data4)
dim(rna_WGDp12)

#将noNV的matrix保存
saveRDS(noCNV_MPS1ip6, file = "/Users/leishengyue/Desktop/data/2024/scRNC-seq/summary/QCsample/noCNV_MPS1ip6.rds")
saveRDS(noCNV_MPS1ip12, file = "/Users/leishengyue/Desktop/data/2024/scRNC-seq/summary/QCsample/noCNV_MPS1ip12.rds")
saveRDS(noCNV_WGDp6, file = "/Users/leishengyue/Desktop/data/2024/scRNC-seq/summary/QCsample/noCNV_WGDp6.rds")
saveRDS(noCNV_WGDp12, file = "/Users/leishengyue/Desktop/data/2024/scRNC-seq/summary/QCsample/noCNV_WGDp12.rds")
