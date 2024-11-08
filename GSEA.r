library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggridges)
setwd("/Users/experiment/scRNC-seq/summary/DEG/")
data = read.table("MPS1ip52/DEG-MPS1ip52.txt",header=TRUE,sep="\t")
colnames(data)[1]="SYMBOL"
head(data)
gene = data$SYMBOL
gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
gene = dplyr::distinct(gene,SYMBOL,.keep_all=T)
data_all <- data %>% 
  inner_join(gene,by="SYMBOL")
data_all_sort <- data_all %>% 
  arrange(desc(avg_log2FC))
geneList = data_all_sort$avg_log2FC 
names(geneList) <- data_all_sort$ENTREZID 

KEGG_database="hsa"
gsea_kegg <- gseKEGG(geneList, organism = KEGG_database, pvalueCutoff = 0.05)

# library(ReactomePA)
# # 假设 geneList 是您的基因列表，并且基因名是 ENTREZ ID
# # 使用 Reactome 富集分析，适用于显著基因集（如阈值筛选后得到的基因子集），可以分析特定子集的富集情况
# gsea_reactome <- enrichPathway(gene = names(geneList), organism = "human", pvalueCutoff = 0.05)

#基因数据是连续排序的数值数据（如 log2 fold change）
data(geneList,package='DOSE')
gsea_reactome <- gsePathway(geneList, nPerm=10000,pvalueCutoff=0.05, pAdjustMethod="BH", verbose=FALSE)
gsea_reactome <- as.data.frame(gsea_reactome)


# # 使用 GO BP 数据库进行 GSEA 分析
# gsea_GO <- gseGO(geneList, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)

# # 设置结果为可读形式
# gsea_GO <- setReadable(gsea_GO, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')


barplot(KEGG,showCategory= 40, title="KEGG pathway")


dotplot(gsea_)
ridgeplot(gsea,label_format = 100)
gseaplot2(gsea,1,pvalue_table = T) 

# 将GSEA结果转换为数据框
gsea_kegg <- as.data.frame(gsea_kegg@result)

gsea_reactome <- as.data.frame(gsea_reactome@result)

# 导出为TXT文件
write.table(gsea_kegg, file = "kegg_MPS1ip52.txt", sep = "\t", row.names = FALSE, quote = FALSE)

write.table(gsea_reactome, file = "reactome_MPS1ip52.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# 绘制与“cell migration”相关的 GSEA 轨迹图
gseaplot2(gsea_GO, geneSetID = "GO:0016477", pvalue_table = TRUE, 
          title = "cell migration")

# 检查 GO ID 是否存在于 GSEA 结果中
any(gsea_GO@result$ID == "GO:0010634")

# 使用更广泛的关键词进行筛选
repair_related <- gsea_GO@result %>%
  filter(grepl("repair|DNA damage|double-strand break|response to DNA damage", Description, ignore.case = TRUE)) %>%
  arrange(p.adjust)

# 检查是否找到相关通路
print(repair_related)

