library(ggplot2)
library(ggpubr)
library(ggsci)
library(dplyr)
library(maftools)
library(xlsx2dfs)
library(ComplexHeatmap)
library(circlize)

####CDK4/6
rm(list = ls())
graphics.off()

load("FUSCC_SNF_RW_Validation_20250401.RData")

GSEA_path <- Enrichment_dataframe_CDK
GSEA_path <- GSEA_path [,-1]
GSEA_path$GROUP <- as.character(GSEA_path$GROUP)
GSEA_path <- GSEA_path[order(GSEA_path$NES,decreasing = T),]
ggbarplot(GSEA_path , x = "NAME" , y = "NES",color = "#8BC37C")+
  ggtitle("Changes of pathway in PDO after drug treatment ")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(angle=90, size=8))+
  labs(x="KEGG signaling pathway",y="Normalization enrichment score")

####PARPi
rm(list = ls())
graphics.off()

load("FUSCC_SNF_RW_Validation_20250401.RData")

GSEA_path <- Enrichment_dataframe_PARP
GSEA_path <- GSEA_path [,-1]
GSEA_path$GROUP <- as.character(GSEA_path$GROUP)
GSEA_path <- GSEA_path[order(GSEA_path$NES,decreasing = T),]
ggbarplot(GSEA_path , x = "NAME" , y = "NES",color = "#4E9ACF")+
  ggtitle("Changes of pathway in PDO after drug treatment ")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(angle=90, size=8))+
  labs(x="KEGG signaling pathway",y="Normalization enrichment score")
