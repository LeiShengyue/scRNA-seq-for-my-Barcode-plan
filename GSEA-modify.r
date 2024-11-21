rm(list = ls())
library(tidyverse)
library(clusterProfiler)
library(corrplot)
library(ggforce)
library(RColorBrewer)

#读取数据
setwd("/Users/experiment/scRNC-seq/summary/DEG/")
GSEAdata <- data.table::fread('Passage52-DEG/GSEAmodify-KEGG-diffe.csv',data.table = F)


data<-as.data.frame(GSEAdata)%>%
  dplyr::select('Description','ONTOLOGY','FoldEnrichment','qvalue','ID')%>%
  mutate(len=str_length(Description))%>%
  dplyr::filter(len<120)%>%
  rowwise()%>%
  mutate(FoldEnrichment=round(eval(parse(text=FoldEnrichment)),3))%>%
  arrange(qvalue)%>%
  # group_by(ONTOLOGY)%>%
  # mutate(ID=1:n())%>%
  # top_n(5,wt=-ID)%>%

  mutate(ONTOLOGY=factor(ONTOLOGY,levels=c("KEGG-Up-MPS1i","KEGG-Up-WGD","KEGG-Down-MPS1i", "KEGG-Down-WGD")),
         Description=str_wrap(Description,width=60),
         Description=factor(Description,levels=rev(Description)))

color <- c("#EFA39F","#f8af76","#58b1ff",  "#9ee2b0")
ggplot(data,aes(-log10(qvalue),ID))+
  geom_col(aes(y=Description,fill=ONTOLOGY),alpha=0.5,show.legend=F)+
  geom_line(aes(x=FoldEnrichment,y=ID,group=1),col='black',size = 1,orientation="y",show.legend = F)+
  geom_point(aes(x=FoldEnrichment,y=ID,fill = ONTOLOGY),size=4,color='black',shape = 21,show.legend=F)+
  facet_wrap("ONTOLOGY",scales="free", nrow = 2)+
  scale_fill_manual(values=color)+
  labs(x="FoldEnrichment and -log10(FDR)",y=NULL)+
  scale_x_continuous(limits=c(0,max(data$FoldEnrichment)))+
  theme(strip.text = element_text(size = 12,face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 12,color="black"),
        axis.title.x = element_text(size = 14,color="black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color="black",
                                    size = 0.8,linetype = 1)) +
  scale_y_discrete(labels=function(y) str_wrap(y, width=30))

# 范围那里也可以改成scale_x_continuous(limits=c(0,max(-log10(data$qvalue)+1)))
                   
