# BCD analysis

# Load data files
setwd("/path/to/your/BCD/Passage6+12_linkage/")  # Set your working directory

MPS1ip6 <- read.table("MPS1ip6-BCD.txt", header = FALSE, stringsAsFactors = FALSE)
MPS1ip12 <- read.table("MPS1ip12-BCD.txt", header = FALSE, stringsAsFactors = FALSE)
WGDp6 <- read.table("WGDp6-BCD.txt", header = FALSE, stringsAsFactors = FALSE)
WGDp12 <- read.table("WGDp12-BCD.txt", header = FALSE, stringsAsFactors = FALSE)

# Set column names
colnames(MPS1ip6) <- c("Cell", "BCD")
colnames(MPS1ip12) <- c("Cell", "BCD")
colnames(WGDp6) <- c("Cell", "BCD")
colnames(WGDp12) <- c("Cell", "BCD")

# Calculate frequency of each category
frequen_MPS1ip6 <- table(MPS1ip6$BCD)
frequen_MPS1ip12 <- table(MPS1ip12$BCD)
frequen_WGDp6 <- table(WGDp6$BCD)
frequen_WGDp12 <- table(WGDp12$BCD)

# Convert frequency counts to frequencies
frequency_MPS1ip6 <- as.data.frame(frequen_MPS1ip6)
frequency_MPS1ip6$Frequency <- frequency_MPS1ip6$Freq / sum(frequency_MPS1ip6$Freq)

frequency_MPS1ip12 <- as.data.frame(frequen_MPS1ip12)
frequency_MPS1ip12$Frequency <- frequency_MPS1ip12$Freq / sum(frequency_MPS1ip12$Freq)

frequency_WGDp6 <- as.data.frame(frequen_WGDp6)
frequency_WGDp6$Frequency <- frequency_WGDp6$Freq / sum(frequency_WGDp6$Freq)

frequency_WGDp12 <- as.data.frame(frequen_WGDp12)
frequency_WGDp12$Frequency <- frequency_WGDp12$Freq / sum(frequency_WGDp12$Freq)

# View results
head(frequency_MPS1ip6)
head(frequency_MPS1ip12)
head(frequency_WGDp6)
head(frequency_WGDp12)

BCD_MPS1ip6 <- frequency_MPS1ip6[[1]]
BCD_MPS1ip12 <- frequency_MPS1ip12[[1]]
BCD_WGDp6 <- frequency_WGDp6[[1]]
BCD_WGDp12 <- frequency_WGDp12[[1]]

# Find common gene names in Passage 12 and Passage 6
common_MPS1i <- BCD_MPS1ip12[BCD_MPS1ip12 %in% BCD_MPS1ip6]
common_WGD <- BCD_WGDp12[BCD_WGDp12 %in% BCD_WGDp6]

# Extract rows containing common genes
filtered_MPS1ip12 <- frequency_MPS1ip12[frequency_MPS1ip12[[1]] %in% common_MPS1i, ]
filtered_WGDp12 <- frequency_WGDp12[frequency_WGDp12[[1]] %in% common_WGD, ]

head(filtered_MPS1ip12)
head(filtered_WGDp12)

# Save frequency data to files
write.table(frequency_MPS1ip6, quote=FALSE, file="frequency_MPS1ip6.txt", sep="\t")
write.table(frequency_MPS1ip12, quote=FALSE, file="frequency_MPS1ip12.txt", sep="\t")
write.table(filtered_MPS1ip12, quote=FALSE, file="frequency_MPS1ip12_filtered.txt", sep="\t")

write.table(frequency_WGDp6, quote=FALSE, file="frequency_WGDp6.txt", sep="\t")
write.table(frequency_WGDp12, quote=FALSE, file="frequency_WGDp12.txt", sep="\t")
write.table(filtered_WGDp12, quote=FALSE, file="frequency_WGDp12_filtered.txt", sep="\t")

# Load required libraries
library(MullerPlot)

# Set working directory for Muller plot generation
setwd("/path/to/your/MullerPlot-Generation/")  # Set your working directory

# Load data files for Muller plot
PopulationDataList <- read.table("MPS1i-PopulationDataList.txt", header = TRUE)
Attributes <- read.table("MPS1i-attributes.txt", header = TRUE, sep = "\t")

# Verify data loading
head(PopulationDataList)
head(Attributes)

# Generate Muller plot
layout(matrix(c(1, 2), nrow = 1), widths = c(0.8, 0.2))
par(mar = c(3, 3, 1, 0))
m.plot <- Muller.plot(attributes = Attributes, population.data = PopulationDataList,
                      data.method = "list", time.interval.method = "equal", mgp = c(2, 0.5, 0))
