library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggridges)
setwd("/Users/scRNC-seq/GSEA/")
data = read.table("MPS1ip6_DEG.txt",header=TRUE,sep="\t")
colnames(data)[1]="SYMBOL"
head(data)
gene = data$SYMBOL
gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
gene = dplyr::distinct(gene,SYMBOL,.keep_all=T)
data_all <- data %>% 
  inner_join(gene,by="SYMBOL")
data_all_sort <- data_all %>% 
  arrange(desc(log2FC))
geneList = data_all_sort$log2FC 
names(geneList) <- data_all_sort$ENTREZID 
KEGG_database="hsa"
gsea <- gseKEGG(geneList, organism = KEGG_database, pvalueCutoff = 0.05)
gsea<- setReadable(gsea, OrgDb=org.Hs.eg.db,keyType = 'ENTREZID')
dotplot(gsea)
ridgeplot(gsea,label_format = 100)
gseaplot2(gsea,1,pvalue_table = T) 
