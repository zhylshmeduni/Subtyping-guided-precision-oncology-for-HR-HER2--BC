library(xlsx2dfs)
library(dplyr)
library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)

rm(list = ls())
graphics.off()
load("CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20221010.Rdata")
rm(list = ls()[which(ls() != "PD_L1_dataframe")])

color_SNF <- c("SNF1" = "#3D76AE", "SNF2" = "#53AD4A",
               "SNF3" = "#EDAB3C", "SNF4" = "#C3392A")

colnames(PD_L1_dataframe)[2] <- "SNF_Subtype"
PD_L1_dataframe$SNF_Subtype <- factor(PD_L1_dataframe$SNF_Subtype, 
                                       levels = c("SNF1", "SNF2", "SNF3", "SNF4"))

stat.test <- PD_L1_dataframe %>%
  t_test(CPS ~ SNF_Subtype)
stat.test <- filter(stat.test, p.adj.signif != "ns")
stat.test <- arrange(stat.test, group2)

max_point <- max(PD_L1_dataframe$CPS)

p <- ggboxplot(PD_L1_dataframe, x = "SNF_Subtype", y = "CPS",
               color = "SNF_Subtype", add = "jitter")+
  theme_bw()+
  theme(panel.grid.minor=element_blank())+
  scale_y_continuous(breaks = c(0, 10, 25, 50, 75, 100))+
  geom_hline(yintercept = 10, linetype = "dashed", color = "black")+
  stat_pvalue_manual(stat.test, 
                     y.position = seq(from = 1.05*max_point, to = (nrow(stat.test)*0.05+1)*max_point, by = 0.05*max_point),
                     label = "{p.adj.signif}", tip.length = 0)+
  scale_color_manual(values = color_SNF)