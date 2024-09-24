library(readxl)
library(ggplot2)
library(RColorBrewer)

# Set working directory
setwd("/experiment/scRNC-seq/summary/CNV+ DEG")

# Read data
data = read.table("DEG_CNV+MPS1i.txt", header = TRUE)
data = as.data.frame(data)

# Set factor levels for Cell_type
data$Cell_type = factor(data$Cell_type, levels = c('MPS1ip1', 'MPS1ip6', 'MPS1ip12'))

# Filter the data based on gene list
data <- data[data$Gene %in% c("CCL1", "POLD1", "PCNA", "MCM6", "MCM7", "GSTO2", 
                              "ESPL1", "DBF4B", "CCNE1", "CDK1", "NCL"), ]

# Set color palette
cols <- rev(brewer.pal(11, 'PiYG'))

# Plot
P6 <- ggplot(data, aes(y = Gene, x = Cell_type)) +
  geom_point(aes(size = p_val_adj, fill = avg_log2FC), shape = 21, color = "black") +   # Point size represents p-value and fill color represents log2FC
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Color gradient for log2FC values
  scale_size_continuous(range = c(4, 8)) +  # Adjust point size range
  labs(x = "Group", y = "Gene", title = "",
       color = expression(log2FC, size = "FDR")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10, color = "black"),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        panel.grid = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

P6

# Save plot to PDF
pdf("dot.pdf", width = 8, height = 8)
P6
dev.off()
