#抽选出CNV+的矩阵
#读取已处理的对象
data1 <- readRDS("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/MPS1ip6.rds")
data2 <- readRDS("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/MPS1ip1.rds")

# 提取特定层的数据
counts_MPS1ip6 <- data1[["RNA"]]$`counts.Gene Expression.MPS1ip6`
counts_MPS1ip1 <- data2[["RNA"]]$`counts.Gene Expression.MPS1ip1`

# 提取前缀
prefixes1 <- sub("_.+$", "", colnames(counts_MPS1ip6))
prefixes2 <- sub("_.+$", "", colnames(counts_MPS1ip1))

# 统计前缀数量
prefix1_counts <- table(prefixes1)
prefix2_counts <- table(prefixes2)

# 打印前缀数量
print(c(prefix1_counts, prefix2_counts))

# 删除前缀 "MPS1ip6_" 和 "MPS1i_"
colnames(counts_MPS1ip6) <- sub("^(MPS1ip6_)", "", colnames(counts_MPS1ip6))
# 删除前缀 "MPS1ip6_" 和 "MPS1i_"
colnames(counts_MPS1ip1) <- sub("^(MPS1ip1_)", "", colnames(counts_MPS1ip1))

str(counts_MPS1ip6)
head(counts_MPS1ip6)

str(counts_MPS1ip1)
head(counts_MPS1ip1)


# 读取txt文件并提取第二列
#c_data1 <- read.table("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ subcluster-MPS1i/MPS1ip6-s23.txt", header = TRUE)
#second_column <- txt_data[, 2]
#library(readxl)

c_data1 <- read.table("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ subcluster-MPS1i/MPS1ip6-s23.txt", header = TRUE)
c_data2 <- read.table("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ subcluster-MPS1i/MPS1ip6-s5.txt", header = TRUE)
c_data3 <- read.table("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ subcluster-MPS1i/MPS1ip1-s10.txt", header = TRUE)
c_data4 <- read.table("/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ subcluster-MPS1i/MPS1ip1-s20.txt", header = TRUE)

second_column1 <- c_data1[, 2]
second_column2 <- c_data2[, 2]
second_column3 <- c_data3[, 2]
second_column4 <- c_data4[, 2]

# 打印查看提取的第二列文字
head(second_column1)
head(second_column2)
head(second_column3)
head(second_column4)

# 打印查看原始matrix列名
head(colnames(counts_MPS1ip6))
head(colnames(counts_MPS1ip1))

# 去除列名和second_column中的空格和其他特殊字符
clean_colnames6 <- gsub("[[:space:]]", "", colnames(counts_MPS1ip6))
clean_colnames1 <- gsub("[[:space:]]", "", colnames(counts_MPS1ip1))
clean_second_column1 <- gsub("[[:space:]]", "", second_column1)
clean_second_column2 <- gsub("[[:space:]]", "", second_column2)
clean_second_column3 <- gsub("[[:space:]]", "", second_column3)
clean_second_column4 <- gsub("[[:space:]]", "", second_column4)

# 保留列名中含有第二列文字的列
s23_MPS1ip6 <- counts_MPS1ip6[, clean_colnames6 %in% clean_second_column1]
s5_MPS1ip6 <- counts_MPS1ip6[, clean_colnames6 %in% clean_second_column2]
s10_MPS1ip1 <- counts_MPS1ip1[, clean_colnames1 %in% clean_second_column3]
s20_MPS1ip1 <- counts_MPS1ip1[, clean_colnames1 %in% clean_second_column4]

# 打印查看过滤后的matrix
head(s23_MPS1ip6)
head(s5_MPS1ip6)
head(s10_MPS1ip1)
head(s20_MPS1ip1)

# 检查过滤后的matrix结构
str(s23_MPS1ip6)
str(s5_MPS1ip6)
str(s10_MPS1ip1)
str(s20_MPS1ip1)

#保存含有CNV的rds矩阵文件
saveRDS(s23_MPS1ip6, file = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ matrix data/CNV+ subcluster matrix/s23-MPS1ip6.rds")
saveRDS(s5_MPS1ip6, file = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ matrix data/CNV+ subcluster matrix/s5-MPS1ip6.rds")
saveRDS(s10_MPS1ip1, file = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ matrix data/CNV+ subcluster matrix/s10-MPS1ip1.rds")
saveRDS(s20_MPS1ip1, file = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/CNV+ matrix data/CNV+ subcluster matrix/s20-MPS1ip1.rds")