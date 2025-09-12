#画dotplot

library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(gghighlight)

gainorloss <- read.table("/Users/leishengyue/Desktop/data/2025/genes on CNV/CNV gainorloss.txt", sep = "\t",
  header = TRUE,
  row.names = 1)

CNVpro <- read.table("/Users/leishengyue/Desktop/data/2025/genes on CNV/topCNV pro.txt", sep = "\t",
  header = TRUE,
  row.names = 1)
library(dplyr)
library(tidyr)
library(ggplot2)

# 假设 gainorloss 和 CNVpro 行列名完全一致
# 行 = 基因 / 通路, 列 = 样本 / 分组

# 1) 转成长表格
df_gain <- gainorloss %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "variation")

df_pro <- CNVpro %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "proportion")

# 2) 合并两个表
df <- left_join(df_gain, df_pro, by = c("gene", "sample"))

# 3) 画 dotplot
ggplot(df, aes(x = sample, y = gene)) +
  geom_point(aes(size = proportion,
                 color = factor(variation))) +
  scale_color_manual(values = c(
    "-1" = "#6BAED6",
     "0" = "#ffffff",
     "1" = "#F4A582")) +
  scale_size_continuous(range = c(1, 6)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color="black"),
        axis.text.y = element_text(face = "bold", color="black"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "bottom",
        legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=10, face="bold"))
