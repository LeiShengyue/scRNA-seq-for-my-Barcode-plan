suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
#suppressMessages(require(DoubletFinder))
library(hdf5r)

# 使用Read10X_h5函数读取h5文件
CTRL <- Seurat::Read10X_h5(filename = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/20240415/per_sample_outs/1/count/1.h5", use.names = TRUE)
MPS1ip1 <- Seurat::Read10X_h5(filename = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/20240415/per_sample_outs/3/count/3.h5", use.names = TRUE)
MPS1ip6 <- Seurat::Read10X_h5(filename = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/20240612/per_sample_outs/3/count/3.h5", use.names = TRUE)
MPS1ip12 <- Seurat::Read10X_h5(filename = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/20240612/per_sample_outs/6/count/6.h5", use.names = TRUE)
WGDp1 <- Seurat::Read10X_h5(filename = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/20240612/per_sample_outs/2/count/2.h5", use.names = TRUE)
WGDp7 <- Seurat::Read10X_h5(filename = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/20240612/per_sample_outs/5/count/5.h5", use.names = TRUE)

# 使用CreateSeuratObject函数构建Seurat对象
sdata.CTRL <- CreateSeuratObject(CTRL, project = "CTRL")
sdata.MPS1ip1<- CreateSeuratObject(MPS1ip1, project = "MPS1ip1")
sdata.MPS1ip6 <- CreateSeuratObject(MPS1ip6, project = "MPS1ip6")
sdata.MPS1ip12 <- CreateSeuratObject(MPS1ip12, project = "MPS1ip12")
sdata.WGDp1 <- CreateSeuratObject(WGDp1, project = "WGDp1")
sdata.WGDp7 <- CreateSeuratObject(WGDp7, project = "WGDp7")


# 添加样本分组信息
sdata.CTRL$type <- "CTRL"
sdata.WGDp1$type <- "WGDp1"
sdata.MPS1ip1$type <- "MPS1ip1"
sdata.WGDp7$type <- "WGDp7"
sdata.MPS1ip6$type <- "MPS1ip6"
sdata.MPS1ip12$type <- "MPS1ip12"

# 合并所有样本
alldata <- merge(x = sdata.CTRL, 
                 y = c(sdata.WGDp1, sdata.MPS1ip1, sdata.MPS1ip12, sdata.MPS1ip6, sdata.WGDp7), 
                 add.cell.ids = c("CTRL", "WGDp1", "MPS1ip1", "MPS1ip12", "MPS1ip6", "WGDp7"))

# 线粒体相关基因含量计算，计算线粒体相关基因百分比
# 使用Seurat自带函数计算
# 使用PercentageFeatureSet函数计算线粒体基因含量
alldata <- PercentageFeatureSet(alldata, "^MT-", col.name = "percent_mito")

# 核糖体相关基因计算，也和上面线粒体计算一样，可以直接使用Seurat中的自带函数进行计算，也可以自行计算，结果是一致的
# 自带函数PercentageFeatureSet
alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")

# 血红细胞相关基因计算
# 自带函数PercentageFeatureSet
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")

# 质控结果的可视化
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")

# 一般用小提琴图对结果进行质控
VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + 
  NoLegend()

# 基因数量随细胞数量的变化关系
FeatureScatter(alldata, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)

# 线粒体基因数量随细胞数量的变化关系
FeatureScatter(alldata, "nCount_RNA", "percent_mito", group.by = "orig.ident", pt.size = 0.5)

# 同时选择符合线粒体和核糖体基因标准的细胞
selected_cells <- WhichCells(alldata, expression = percent_mito < 10 & percent_ribo > 0.1)

# 根据选择的细胞子集化对象
alldata <- subset(alldata, cells = selected_cells)

# 验证子集化后的对象
print(dim(alldata))
head(alldata@meta.data)

# 子集化 Control 组
CTRL_cells <- WhichCells(alldata, idents = "CTRL")
CTRL_data <- subset(alldata, cells = CTRL_cells)
saveRDS(CTRL_data, file = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/CTRL.rds")

# 子集化 WGD 组
WGDp1_cells <- WhichCells(alldata, idents = "WGDp1")
WGDp1_data <- subset(alldata, cells = WGDp1_cells)
saveRDS(WGDp1_data, file = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/WGDp1.rds")

# 子集化 MPS1i 组
MPS1ip1_cells <- WhichCells(alldata, idents = "MPS1ip1")
MPS1ip1_data <- subset(alldata, cells = MPS1ip1_cells)
saveRDS(MPS1ip1_data, file = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/MPS1ip1.rds")

# 子集化 MPS1i 组
MPS1ip6_cells <- WhichCells(alldata, idents = "MPS1ip6")
MPS1ip6_data <- subset(alldata, cells = MPS1ip6_cells)
saveRDS(MPS1ip6_data, file = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/MPS1ip6.rds")

# 子集化 MPS1i 组
MPS1ip12_cells <- WhichCells(alldata, idents = "MPS1ip12")
MPS1ip12_data <- subset(alldata, cells = MPS1ip12_cells)
saveRDS(MPS1ip12_data, file = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/MPS1ip12.rds")

# 子集化 MPS1i 组
WGDp7_cells <- WhichCells(alldata, idents = "WGDp7")
WGDp7_data <- subset(alldata, cells = WGDp7_cells)
saveRDS(WGDp7_data, file = "/Users/asakawayuuki/Desktop/がん研/experiment/scRNC-seq/summary/QCsample/WGDp7.rds")
