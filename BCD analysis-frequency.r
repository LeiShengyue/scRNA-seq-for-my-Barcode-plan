#BCD analysis
# 读取数据
setwd("/Users/scRNC-seq/BCD/Passage6+12 紐付け/")

MPS1ip6 <- read.table("MPS1ip6-BCD.txt", header = FALSE, stringsAsFactors = FALSE)
MPS1ip12 <- read.table("MPS1ip12-BCD.txt", header = FALSE, stringsAsFactors = FALSE)
WGDp6 <- read.table("WGDp6-BCD.txt", header = FALSE, stringsAsFactors = FALSE)
WGDp12 <- read.table("WGDp12-BCD.txt", header = FALSE, stringsAsFactors = FALSE)

# 设定列名称
colnames(MPS1ip6) <- c("Cell", "BCD")
colnames(MPS1ip12) <- c("Cell", "BCD")
colnames(WGDp6) <- c("Cell", "BCD")
colnames(WGDp12) <- c("Cell", "BCD")

# 计算每个类别的频数
frequen_MPS1ip6 <- table(MPS1ip6$BCD)
frequen_MPS1ip12 <- table(MPS1ip12$BCD)
frequen_WGDp6 <- table(WGDp6$BCD)
frequen_WGDp12 <- table(WGDp12$BCD)

# 将频数转换为频率
frequency_MPS1ip6 <- as.data.frame(frequen_MPS1ip6)
frequency_MPS1ip6$Frequency <- frequency_MPS1ip6$Freq / sum(frequency_MPS1ip6$Freq)

# 将频数转换为频率
frequency_MPS1ip12 <- as.data.frame(frequen_MPS1ip12)
frequency_MPS1ip12$Frequency <- frequency_MPS1ip12$Freq / sum(frequency_MPS1ip12$Freq)

# 将频数转换为频率
frequency_WGDp6 <- as.data.frame(frequen_WGDp6)
frequency_WGDp6$Frequency <- frequency_WGDp6$Freq / sum(frequency_WGDp6$Freq)

# 将频数转换为频率
frequency_WGDp12 <- as.data.frame(frequen_WGDp12)
frequency_WGDp12$Frequency <- frequency_WGDp12$Freq / sum(frequency_WGDp12$Freq)


library(ggplot2)
library(dplyr)

# 读取 CSV 文件（如果你有文件）
# df <- read.csv("你的文件路径.csv")
# df <- frequency_MPS1ip6

# 统计每个 Freq 值有多少个 barcode
freq_distribution <- df %>%
  group_by(Freq) %>%
  summarise(Barcode_Count = n())

# 画柱状图
ggplot(freq_distribution, aes(x = Freq, y = Barcode_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
        geom_text(aes(label=Barcode_Count), vjust=0)+
  # scale_y_log10() +  # 如果数据跨度大，建议用 log 轴
  theme_minimal() +
  labs(x = "Barcode frequency", 
       y = "Barcode numbers",
       title = "Barcode heterogeneity")

library(ggplot2)
library(dplyr)

# 计算各 Freq 占比
freq_distribution <- df %>%
  group_by(Freq) %>%
  summarise(Barcode_Count = n()) %>%
  mutate(Percentage = Barcode_Count / sum(Barcode_Count) * 100)  # 计算百分比

# 画饼图
ggplot(freq_distribution, aes(x = "", y = Percentage, fill = as.factor(Freq))) +
  geom_bar(stat = "identity", width = 1) +geom_text(aes(label=Percentage), vjust=0)+
  coord_polar(theta = "y") +  # 变成饼图
  theme_void() +  # 移除背景和坐标轴
  labs(fill = "Barcode Frequency", title = "Barcode Frequency Distribution")

