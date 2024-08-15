# Extract the CNV+ matrix
# Load processed objects
data1 <- readRDS("/path/to/MPS1ip6.rds")
data2 <- readRDS("/path/to/MPS1ip1.rds")

# Extract specific layers of data
counts_MPS1ip6 <- data1[["RNA"]]$counts
counts_MPS1ip1 <- data2[["RNA"]]$counts

# Remove prefixes from column names
colnames(counts_MPS1ip6) <- sub("^(MPS1ip6_)", "", colnames(counts_MPS1ip6))
colnames(counts_MPS1ip1) <- sub("^(MPS1ip1_)", "", colnames(counts_MPS1ip1))

str(counts_MPS1ip6)
head(counts_MPS1ip6)

str(counts_MPS1ip1)
head(counts_MPS1ip1)

# Load text files and extract the second column
c_data1 <- read.table("/path/to/MPS1ip6-s23.txt", header = TRUE)
c_data2 <- read.table("/path/to/MPS1ip6-s5.txt", header = TRUE)
c_data3 <- read.table("/path/to/MPS1ip1-s10.txt", header = TRUE)
c_data4 <- read.table("/path/to/MPS1ip1-s20.txt", header = TRUE)

second_column1 <- c_data1[, 2]
second_column2 <- c_data2[, 2]
second_column3 <- c_data3[, 2]
second_column4 <- c_data4[, 2]

# Clean column names and second column values by removing spaces and special characters
clean_colnames6 <- gsub("[[:space:]]", "", colnames(counts_MPS1ip6))
clean_colnames1 <- gsub("[[:space:]]", "", colnames(counts_MPS1ip1))
clean_second_column1 <- gsub("[[:space:]]", "", second_column1)
clean_second_column2 <- gsub("[[:space:]]", "", second_column2)
clean_second_column3 <- gsub("[[:space:]]", "", second_column3)
clean_second_column4 <- gsub("[[:space:]]", "", second_column4)

# Retain columns where column names match the second column values
s23_MPS1ip6 <- counts_MPS1ip6[, clean_colnames6 %in% clean_second_column1]
s5_MPS1ip6 <- counts_MPS1ip6[, clean_colnames6 %in% clean_second_column2]
s10_MPS1ip1 <- counts_MPS1ip1[, clean_colnames1 %in% clean_second_column3]
s20_MPS1ip1 <- counts_MPS1ip1[, clean_colnames1 %in% clean_second_column4]

# Check the filtered matrices
head(s23_MPS1ip6)
head(s5_MPS1ip6)
head(s10_MPS1ip1)
head(s20_MPS1ip1)

# Check the structure of the filtered matrices
str(s23_MPS1ip6)
str(s5_MPS1ip6)
str(s10_MPS1ip1)
str(s20_MPS1ip1)

# Save the CNV+ matrices as RDS files
saveRDS(s23_MPS1ip6, file = "/path/to/s23-MPS1ip6.rds")
saveRDS(s5_MPS1ip6, file = "/path/to/s5-MPS1ip6.rds")
saveRDS(s10_MPS1ip1, file = "/path/to/s10-MPS1ip1.rds")
saveRDS(s20_MPS1ip1, file = "/path/to/s20-MPS1ip1.rds")
