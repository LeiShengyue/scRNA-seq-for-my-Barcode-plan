#画dotplot

library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(gghighlight)
library(RColorBrewer)


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
  rownames_to_column("Chr") %>%
  pivot_longer(-Chr, names_to = "Clone", values_to = "variation")

df_pro <- CNVpro %>%
  rownames_to_column("Chr") %>%
  pivot_longer(-Chr, names_to = "Clone", values_to = "proportion")

# 2) 合并两个表
df <- left_join(df_gain, df_pro, by = c("Chr", "Clone"))

# # 3) 画 dotplot
# ggplot(df, aes(x = Clone, y = Chr)) +
#   geom_point(aes(size = proportion,
#                  color = factor(variation))) +
#   scale_color_manual(values = c(
#     "-1" = "#6BAED6",
#      "0" = "#dbdbdb",
#      "1" = "#F4A582")) +
#   scale_size_continuous(range = c(1, 6)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color="black"),
#         axis.text.y = element_text(face = "bold", color="black"),
#         axis.title = element_text(size=12, face="bold"),
#         legend.position = "bottom",
#         legend.title = element_text(size=10, face="bold"),
#         legend.text = element_text(size=10, face="bold"))


df <- df %>%
  mutate(proportion_scaled = (proportion - min(proportion, na.rm=TRUE)) / 
                               (max(proportion, na.rm=TRUE) - min(proportion, na.rm=TRUE)))

# 假设 gainorloss 和 CNVpro 行列名完全一致
# 行 = 基因 / 通路, 列 = 样本 / 分组

# 1) 转成长表格
df_gain <- gainorloss %>%
  rownames_to_column("Chr") %>%
  pivot_longer(-Chr, names_to = "Clone", values_to = "variation")

df_pro <- CNVpro %>%
  rownames_to_column("Chr") %>%
  pivot_longer(-Chr, names_to = "Clone", values_to = "proportion")

# 2) 合并两个表
df <- left_join(df_gain, df_pro, by = c("Chr", "Clone"))

# # 3) 画 dotplot
# ggplot(df, aes(x = Clone, y = Chr)) +
#   geom_point(aes(size = proportion,
#                  color = factor(variation))) +
#   scale_color_manual(values = c(
#     "-1" = "#6BAED6",
#      "0" = "#dbdbdb",
#      "1" = "#F4A582")) +
#   scale_size_continuous(range = c(1, 6)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color="black"),
#         axis.text.y = element_text(face = "bold", color="black"),
#         axis.title = element_text(size=12, face="bold"),
#         legend.position = "bottom",
#         legend.title = element_text(size=10, face="bold"),
#         legend.text = element_text(size=10, face="bold"))


df <- df %>%
  mutate(proportion_scaled = (proportion - min(proportion, na.rm=TRUE)) / 
                               (max(proportion, na.rm=TRUE) - min(proportion, na.rm=TRUE)))

# library(paletteer)

# ggplot(df, aes(x = Clone, y = Chr)) +
#   geom_point(aes(size = proportion_scaled, color = variation)) +
#   scale_color_gradientn(
#     colours = rev(paletteer::paletteer_d("MetBrewer::Hiroshige")),  # 反转
#     limits = c(1.5, 4.5),  
#     oob = scales::squish
#   ) +
#   scale_size_continuous(
#     range = c(0.5, 8),
#     breaks = c(0, 0.25, 0.50, 0.75, 1.00),
#     labels = scales::label_number(accuracy = 0.01)
#   ) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color="black"),
#     axis.text.y = element_text(face = "bold", color="black"),
#     axis.title  = element_text(size=12, face="bold"),
#     legend.position = "bottom",
#     legend.title    = element_text(size=10, face="bold"),
#     legend.text     = element_text(size=10, face="bold")
#   )

#中间改成白色
ggplot(df, aes(x = Clone, y = Chr)) +
    geom_point(aes(size = proportion_scaled, color = variation)) +
    scale_color_gradientn(
        colours = rev(c("#A50026FF", "#D73027FF", "#F46D43FF", "#FDAE61FF", "#FEE090FF", "white", "#E0F3F8FF", "#ABD9E9FF", "#74ADD1FF", "#4575B4FF", "#313695FF")),  # 反转
        limits = c(1.5, 4.5),  
        oob = scales::squish
    ) +
    scale_size_continuous(
        range = c(0.5, 8),
        breaks = c(0, 0.25, 0.50, 0.75, 1.00),
        labels = scales::label_number(accuracy = 0.01)
    ) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color="black"),
        axis.text.y = element_text(face = "bold", color="black"),
        axis.title  = element_text(size=12, face="bold"),
        legend.position = "bottom",
        legend.title    = element_text(size=10, face="bold"),
        legend.text     = element_text(size=10, face="bold")
    )
