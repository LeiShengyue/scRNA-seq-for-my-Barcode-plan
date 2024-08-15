library(fgsea)
library(readxl)

# Load the gene expression data
Tumor1 <- read.delim("/path/to/MPS1ip6.txt") # 15. Tumor_Slow_Aneuploid_Euploid.xlsx

# Convert matrix to list function
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

# Check the first few rows of the data
head(Tumor1)

# Calculate rankings based on avg_log2FC and FDR
rankings <- Tumor1$avg_log2FC * Tumor1$FDR
names(rankings) <- Tumor1$gene

# Check the first few rankings
head(rankings)

# Sort the gene list in descending order
Tumor1_genelist <- sort(rankings, decreasing = TRUE)
head(Tumor1_genelist)

# Define the path to the GMT file
gmt_file <- "/path/to/HALLMARK_NOTCH_SIGNALING.v2023.2.Hs.gmt"

# Function to prepare GMT file
prepare_gmt <- function(gmt_file, genes, savefile = FALSE) {
    gmt <- gmtPathways(gmt_file)
    hidden <- unique(unlist(gmt))
    mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                  nrow = length(hidden), ncol = length(gmt))
    for (i in 1:dim(mat)[2]) {
        mat[, i] <- as.numeric(hidden %in% gmt[[i]])
    }
    hidden1 <- intersect(genes, hidden)
    valid_cols <- which(colSums(mat[hidden1, , drop = FALSE]) > 5)
    mat <- mat[hidden1, valid_cols, drop = FALSE] # filter for gene sets with more than 5 genes annotated
    final_list <- matrix_to_list(mat) # use the function we previously defined
    print('Wohoo! .gmt conversion successful! :)')
    return(final_list)
}

# Prepare the GMT file and filter the gene sets
genes <- Tumor1$gene
bg_genes <- prepare_gmt(gmt_file, genes, savefile = FALSE)

# Perform GSEA analysis
GSEA_Tumor1 <- fgsea(bg_genes, rankings)
head(GSEA_Tumor1)
