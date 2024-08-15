suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
# suppressMessages(require(DoubletFinder))
# There are tetraploid samples, so avoid removing tetraploid cells by using DoubletFinder
library(hdf5r)

# Read h5 files using the Read10X_h5 function
CTRL <- Seurat::Read10X_h5(filename = "/path/to/CTRL.h5", use.names = TRUE)
MPS1ip1 <- Seurat::Read10X_h5(filename = "/path/to/MPS1ip1.h5", use.names = TRUE)
MPS1ip6 <- Seurat::Read10X_h5(filename = "/path/to/MPS1ip6.h5", use.names = TRUE)
MPS1ip12 <- Seurat::Read10X_h5(filename = "/path/to/MPS1ip12.h5", use.names = TRUE)
WGDp1 <- Seurat::Read10X_h5(filename = "/path/to/WGDp1.h5", use.names = TRUE)
WGDp7 <- Seurat::Read10X_h5(filename = "/path/to/WGDp7.h5", use.names = TRUE)

# Create Seurat objects using the CreateSeuratObject function
sdata.CTRL <- CreateSeuratObject(CTRL, project = "CTRL")
sdata.MPS1ip1 <- CreateSeuratObject(MPS1ip1, project = "MPS1ip1")
sdata.MPS1ip6 <- CreateSeuratObject(MPS1ip6, project = "MPS1ip6")
sdata.MPS1ip12 <- CreateSeuratObject(MPS1ip12, project = "MPS1ip12")
sdata.WGDp1 <- CreateSeuratObject(WGDp1, project = "WGDp1")
sdata.WGDp7 <- CreateSeuratObject(WGDp7, project = "WGDp7")

# Add sample grouping information
sdata.CTRL$type <- "CTRL"
sdata.WGDp1$type <- "WGDp1"
sdata.MPS1ip1$type <- "MPS1ip1"
sdata.WGDp7$type <- "WGDp7"
sdata.MPS1ip6$type <- "MPS1ip6"
sdata.MPS1ip12$type <- "MPS1ip12"

# Merge all samples
alldata <- merge(x = sdata.CTRL, 
                 y = c(sdata.WGDp1, sdata.MPS1ip1, sdata.MPS1ip12, sdata.MPS1ip6, sdata.WGDp7), 
                 add.cell.ids = c("CTRL", "WGDp1", "MPS1ip1", "MPS1ip12", "MPS1ip6", "WGDp7"))

# Calculate mitochondrial gene content using PercentageFeatureSet
alldata <- PercentageFeatureSet(alldata, "^MT-", col.name = "percent_mito")

# Calculate ribosomal gene content using PercentageFeatureSet
alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")

# Calculate hemoglobin gene content using PercentageFeatureSet
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")

# Visualize QC results
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")

# Use violin plots for QC results
VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + 
  NoLegend()

# Relationship between gene counts and cell counts
FeatureScatter(alldata, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)

# Relationship between mitochondrial gene counts and cell counts
FeatureScatter(alldata, "nCount_RNA", "percent_mito", group.by = "orig.ident", pt.size = 0.5)

# Select cells that meet both mitochondrial and ribosomal gene criteria
selected_cells <- WhichCells(alldata, expression = percent_mito < 10 & percent_ribo > 0.1)

# Subset the object based on selected cells
alldata <- subset(alldata, cells = selected_cells)

# Verify the subsetted object
print(dim(alldata))
head(alldata@meta.data)

# Subset Control group
CTRL_cells <- WhichCells(alldata, idents = "CTRL")
CTRL_data <- subset(alldata, cells = CTRL_cells)
saveRDS(CTRL_data, file = "/path/to/CTRL.rds")

# Subset WGD group
WGDp1_cells <- WhichCells(alldata, idents = "WGDp1")
WGDp1_data <- subset(alldata, cells = WGDp1_cells)
saveRDS(WGDp1_data, file = "/path/to/WGDp1.rds")

# Subset MPS1i group
MPS1ip1_cells <- WhichCells(alldata, idents = "MPS1ip1")
MPS1ip1_data <- subset(alldata, cells = MPS1ip1_cells)
saveRDS(MPS1ip1_data, file = "/path/to/MPS1ip1.rds")

MPS1ip6_cells <- WhichCells(alldata, idents = "MPS1ip6")
MPS1ip6_data <- subset(alldata, cells = MPS1ip6_cells)
saveRDS(MPS1ip6_data, file = "/path/to/MPS1ip6.rds")

MPS1ip12_cells <- WhichCells(alldata, idents = "MPS1ip12")
MPS1ip12_data <- subset(alldata, cells = MPS1ip12_cells)
saveRDS(MPS1ip12_data, file = "/path/to/MPS1ip12.rds")

WGDp7_cells <- WhichCells(alldata, idents = "WGDp7")
WGDp7_data <- subset(alldata, cells = WGDp7_cells)
saveRDS(WGDp7_data, file = "/path/to/WGDp7.rds")
