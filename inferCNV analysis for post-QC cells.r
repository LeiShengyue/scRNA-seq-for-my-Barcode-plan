# Data after QC
library(Seurat)
library(infercnv)
library(Matrix)

# Load data files
data1 <- readRDS("/path/to/CTRL.rds")
data2 <- readRDS("/path/to/MPS1ip1.rds")
data3 <- readRDS("/path/to/MPS1ip6.rds")
data4 <- readRDS("/path/to/MPS1ip12.rds")
data5 <- readRDS("/path/to/WGDp1.rds")
data6 <- readRDS("/path/to/WGDp7.rds")

# Extract gene expression count data
rna_CTRL <- data1[["RNA"]]$counts
rna_MPS1ip1 <- data2[["RNA"]]$counts
rna_MPS1ip6 <- data3[["RNA"]]$counts
rna_MPS1ip12 <- data4[["RNA"]]$counts
rna_WGDp1 <- data5[["RNA"]]$counts
rna_WGDp7 <- data6[["RNA"]]$counts

# Confirm data types
print(class(rna_CTRL))
print(class(rna_MPS1ip1))
print(class(rna_MPS1ip6))
print(class(rna_MPS1ip12))
print(class(rna_WGDp1))
print(class(rna_WGDp7))

# Combine data matrices
combined_matrix <- cbind(rna_CTRL, rna_MPS1ip1, rna_MPS1ip6, rna_MPS1ip12, rna_WGDp1, rna_WGDp7)

# Check column names
column_names <- colnames(combined_matrix)
print(column_names)

# Find duplicated column names
duplicated_columns <- duplicated(column_names)

# Print duplicated column names
duplicated_column_names <- column_names[duplicated_columns]
print(duplicated_column_names)

# Remove duplicated columns
unique_combined_matrix <- combined_matrix[, !duplicated_columns]

# Create Seurat object
seurat_object <- CreateSeuratObject(counts = unique_combined_matrix, project = "10X_Project")

# View Seurat object structure
str(seurat_object)

# Get Seurat expression matrix data
seurat_expression_matrix <- GetAssayData(seurat_object, assay = "RNA", slot = "counts")

# View data
head(Matrix::t(seurat_expression_matrix))

infercnv_expression_matrix <- GetAssayData(seurat_object, assay = "RNA", slot = "counts")

# Extract prefixes
prefixes <- sub("_.+$", "", colnames(infercnv_expression_matrix))

# Count prefixes
prefix_counts <- table(prefixes)

# Print prefix counts
print(prefix_counts)

# Remove prefixes "CTRL_", "MPS1ip1_", "MPS1ip6_", "MPS1ip12_", "WGDp1_", and "WGDp7_"
colnames(infercnv_expression_matrix) <- sub("^(CTRL_|MPS1ip1_|MPS1ip6_|MPS1ip12_|WGDp1_|WGDp7_)", "", colnames(infercnv_expression_matrix))

# View column names
head(colnames(infercnv_expression_matrix))

# Load annotations file and specify column names
annotations <- read.table("/path/to/mixanno.txt", header = FALSE, sep = "\t", col.names = c("Cell_Name", "Annotation"))

# Find cells in annotations file that are not present in infercnv_expression_matrix
extra_cells_in_annotations <- setdiff(annotations$Cell_Name, colnames(infercnv_expression_matrix))

# Print extra cell names
print(extra_cells_in_annotations)

# Remove extra cell names
annotations_filtered <- annotations[!annotations$Cell_Name %in% extra_cells_in_annotations, ]

# Check filtered annotations file
print(paste("Number of cell names in filtered annotations:", nrow(annotations_filtered)))

# Ensure all cell names in filtered annotations file are present in infercnv_expression_matrix
if (!all(annotations_filtered$Cell_Name %in% colnames(infercnv_expression_matrix))) {
    stop("Some cell names in the filtered annotation file are not found in the data matrix.")
}

# Remove leading and trailing spaces
annotations_filtered$Cell_Name <- trimws(annotations_filtered$Cell_Name)
annotations_filtered$Annotation <- trimws(annotations_filtered$Annotation)

# Save filtered annotations file as a temporary file
write.table(annotations_filtered, file = "/path/to/filteredmixanno.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Identify duplicate row names
anno_data <- read.table("/path/to/filteredmix.txt", header = FALSE, sep = "\t")
duplicated_rows <- anno_data[duplicated(anno_data$V1), ]
print(duplicated_rows)

# Remove duplicated rows
unique_anno_data <- anno_data[!duplicated(anno_data$V1), ]
write.table(unique_anno_data, "/path/to/filteredmix_unique.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Create infercnv object
infercnv_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = infercnv_expression_matrix,
    annotations_file = "/path/to/filteredmixanno_unique.txt",
    delim = "\t",
    gene_order_file = "/path/to/hg38_gencode.txt",
    ref_group_names = "CTRL"
)

# Run infercnv
infercnv_obj = infercnv::run(infercnv_obj, 
                             output_format = "pdf",
                             num_threads = 8,
                             cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir = "/path/to/mix-output",  # dir is auto-created for storing outputs
                             window_length = 201, # 51 is way too small
                             cluster_by_groups = TRUE, # cluster by groups
                             HMM = TRUE, # turn on to auto-run the HMM prediction of CNV levels
                             HMM_transition_prob = 1e-6,
                             HMM_report_by = c("subcluster"), # or "cell"
                             analysis_mode = c('subclusters'), # or 'cells'
                             denoise = TRUE,
                             sd_amplifier = 1.35)  # sets midpoint for logistic

# Apply median filtering
infercnv_obj_medianfiltered = infercnv::apply_median_filtering(infercnv_obj)

# Plot CNV
infercnv::plot_cnv(infercnv_obj_medianfiltered,
                   out_dir = "/path/to/mix-output",
                   output_format = "pdf",
                   output_filename = 'infercnv.median_filtered',
                   x.range = "auto",
                   x.center = 1,
                   title = "infercnv",
                   color_safe_pal = FALSE)

