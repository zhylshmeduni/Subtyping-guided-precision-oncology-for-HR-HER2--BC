library(GSVA)
library(GSEABase)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pheatmap)

####Endo_Sen&Resist_Heatmap
rm(list = ls())
graphics.off()
load("CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20221010.Rdata")

pathway = getGmt("./genesets_endo_sen&resist.gmt")

gsva_matrix1 <- gsva(as.matrix(exp.fpkm.TT), pathway, method='ssgsea')  
rownames(gsva_matrix1) = c("Aromatase_inhibitor_responsive_module",
                           "Tamoxifen_or_SERD_responsive_module",
                           "RTK_pathway","EMT_pathway","WNT_pathway","CellCycle_pathway","PI3K_pathway")

pathway_score = as.data.frame(t(gsva_matrix1))
pathway_score$SNF = paste0("SNF",SNF_Cluster[rownames(pathway_score)])
pathway_score$SNF = factor(pathway_score$SNF,levels = c("SNF1","SNF2","SNF3","SNF4"))

summary(aov(Aromatase_inhibitor_responsive_module~SNF,data = pathway_score))
summary(aov(Tamoxifen_or_SERD_responsive_module~SNF,data = pathway_score))
summary(aov(RTK_pathway~SNF,data = pathway_score))
summary(aov(EMT_pathway~SNF,data = pathway_score))
summary(aov(WNT_pathway~SNF,data = pathway_score))
summary(aov(CellCycle_pathway~SNF,data = pathway_score))
summary(aov(PI3K_pathway~SNF,data = pathway_score))


col <- c("#2378b3", "#1cb038","#f8a900","#d5271a")
path <- c("Aromatase_inhibitor_responsive_module",
          "Tamoxifen_or_SERD_responsive_module",
          "RTK_pathway","EMT_pathway","WNT_pathway","CellCycle_pathway","PI3K_pathway")

i <- path[7]
pdf(paste0("./",i,"_.pdf"),height = 8, width = 8)
ggplot(pathway_score,aes(x=SNF,y=pathway_score[,i],color=SNF))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  #stat_compare_means()+
  xlab("")+ylab(i)+
  geom_boxplot(outlier.colour="gray", outlier.shape=16,outlier.size=1.5)+
  #stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)+
  geom_jitter(shape=16, position = position_jitter(0.2))+
  scale_color_manual(values = col)+
  stat_compare_means(method = "kruskal.test", label = "p.format", label.x = 1.5, label.y = max(pathway_score[,i]) * 1.1)
dev.off()
