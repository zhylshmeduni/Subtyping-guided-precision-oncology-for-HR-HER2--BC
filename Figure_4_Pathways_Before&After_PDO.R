library(ggplot2)
library(ggpubr)
library(ggsci)
library(dplyr)
library(maftools)
library(xlsx2dfs)
library(ComplexHeatmap)
library(circlize)

rm(list = ls())
graphics.off()

load("FUSCC_SNF_RW_Validation_20250401.RData")

GSEA_path <- Enrichment_dataframe_CO

###CO BarPlot###
GSEA_path <- GSEA_path [,-1]
GSEA_path$GROUP <- as.character(GSEA_path$GROUP)
GSEA_path <- GSEA_path[order(GSEA_path$NES,decreasing = T),]
ggbarplot(GSEA_path , x = "NAME" , y = "NES",color = "GROUP")+
  ggtitle("Changes of pathway in PDO after drug treatment ")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(angle=90, size=8),
        #axis.title.x=element_text(angle=10, color='red'),
        # axis.title.y=element_text(angle=360, color='blue', face='bold', size=14,vjust=0.5)
  )+
  scale_color_manual("#E15759", "#F1CE63", "#4E79A7", "#59A14F")+
  labs(x="KEGG signaling pathway",y="Normalization enrichment score")