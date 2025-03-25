#把barcode的link信息提出来
#读一下SummarizedExperiment的data
CTRLp52_BCD <- readRDS("/Users/Desktop/がん研/experiment/scRNC-seq/BCD/se_CTRLp52.rds")
MPS1ip52_BCD <- readRDS("/Users/Desktop/がん研/experiment/scRNC-seq/BCD/se_MPS1ip52.rds")
WGDp52_BCD <- readRDS("/Users/Desktop/がん研/experiment/scRNC-seq/BCD/se_WGDp52.rds")

#从里面提出barcode的那一列
CTRLp52_BCD_info <- CTRLp52_BCD@colData
MPS1ip52_BCD_info <- MPS1ip52_BCD@colData
WGDp52_BCD_info <- WGDp52_BCD@colData

#转换一下格式
CTRLp52_BCD_info <- as.data.frame(CTRLp52_BCD_info)
MPS1ip52_BCD_info <- as.data.frame(MPS1ip52_BCD_info)
WGDp52_BCD_info <- as.data.frame(WGDp52_BCD_info)

#看一眼
head(CTRLp52_BCD_info)
head(MPS1ip52_BCD_info)
head(WGDp52_BCD_info)

#保存一下
write.table(CTRLp52_BCD_info ,file = "/Users/Desktop/がん研/experiment/scRNC-seq/BCD/CTRLp52_BCD.txt", row.names = T, sep = "\t")
write.table(MPS1ip52_BCD_info ,file = "/Users/Desktop/がん研/experiment/scRNC-seq/BCD/MPS1ip52_BCD.txt", row.names = T, sep = "\t")
write.table(WGDp52_BCD_info ,file = "/Users/Desktop/がん研/experiment/scRNC-seq/BCD/WGDp52_BCD.txt", row.names = T, sep = "\t")

#提取BCD的矩阵
#读取已处理的对象
setwd("/Users/Desktop/がん研/experiment/scRNC-seq/summary/QCsample")

data1 <- readRDS("MPS1ip52.rds")
data2 <- readRDS("WGDp52.rds")

counts_MPS1ip52 <- data1[["RNA"]]$`counts.Gene Expression.MPS1ip52`
counts_WGDp52 <- data2[["RNA"]]$`counts.Gene Expression.WGDp52`

# 提取前缀
prefixes1 <- sub("_.+$", "", colnames(counts_MPS1ip52))
prefixes2 <- sub("_.+$", "", colnames(counts_WGDp52))

# 统计前缀数量
prefix1_counts <- table(prefixes1)
prefix2_counts <- table(prefixes2)

# 打印前缀数量
print(c(prefix1_counts, prefix2_counts))

colnames(counts_MPS1ip52) <- sub("^(MPS1ip52_)", "", colnames(counts_MPS1ip52))
colnames(counts_WGDp52) <- sub("^(WGDp52_)", "", colnames(counts_WGDp52))

#开始读BCD的clone的细胞
data_MPS1ip52 <- read.table("/Users/Desktop/がん研/experiment/scRNC-seq/BCD/p52_BCDlink/MPS1ip52_topBCD_cell.txt", header = F)
data_WGDp52 <- read.table("/Users/Desktop/がん研/experiment/scRNC-seq/BCD/p52_BCDlink/WGDp52_topBCD_cell.txt", header = F)

MPS1ip52_cell <- data_MPS1ip52[, 1]
WGDp52_cell <- data_WGDp52[, 1]

head(MPS1ip52_cell)
head(WGDp52_cell)

# 打印查看原始matrix列名
head(colnames(counts_MPS1ip52))
head(colnames(counts_WGDp52))

topBCD_MPS1ip52 <- counts_MPS1ip52[, colnames(counts_MPS1ip52) %in% MPS1ip52_cell]
topBCD_WGDp52 <- counts_WGDp52[, colnames(counts_WGDp52) %in% WGDp52_cell]

#看一下格式
dim(topBCD_MPS1ip52)
dim(topBCD_WGDp52)

#保存一下
saveRDS(topBCD_MPS1ip52, file = "/Users/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/topBCD_MPS1ip52.rds")
#保存一下
saveRDS(topBCD_WGDp52, file = "/Users/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/topBCD_WGDp52.rds")
