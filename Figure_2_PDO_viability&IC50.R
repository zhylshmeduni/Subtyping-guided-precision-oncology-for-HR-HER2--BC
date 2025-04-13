library(readxl)
library(ggpubr)
library(pheatmap)
library(reshape2)

####Endocrine_Viability
###Heatmap
rm(list = ls())
graphics.off()
load("FUSCC_SNF_RW_Validation_20250401.RData")

snf1 <- subtset(Endo_1, Endo_1$SNF_byDP_Subtype=="SNF1")[,2:4]
snf2 <- subtset(Endo_1, Endo_1$SNF_byDP_Subtype=="SNF2")[,2:4]
snf3 <- subtset(Endo_1, Endo_1$SNF_byDP_Subtype=="SNF3")[,2:4]
snf4 <- subtset(Endo_1, Endo_1$SNF_byDP_Subtype=="SNF4")[,2:4]

SNF1_mean <- colMeans(snf1,na.rm = T)
SNF2_mean <- colMeans(snf2,na.rm = T)
SNF3_mean <- colMeans(snf3,na.rm = T)
SNF4_mean <- colMeans(snf4,na.rm = T)

Data_Viability_mean <- data.frame(
  row.names = colnames(Endo_1)[2:4],
  SNF1=SNF1_mean,
  SNF2=SNF2_mean,
  SNF3=SNF3_mean,
  SNF4=SNF4_mean
)

pheatmap(Data_Viability_mean,scale = "row",cluster_rows = T,cluster_cols = F)
write.csv(Data_Viability_mean,file = "Endo_Viability_mean.csv")

###Stat
tmp_data_1 <- melt(Endo_1, "SNF_byDP_Subtype",variable_name = "PDO_Viability")
tmp_data_1 <- tmp_data_1[!is.na(tmp_data_1$value),]
tmp_WT <- compare_means(value~SNF_byDP_Subtype,tmp_data_1,group.by = "variable",method = "wilcox.test")
tmp_KW <- compare_means(value~SNF_byDP_Subtype,tmp_data_1,group.by = "variable",method = "kruskal.test")
write.csv(tmp_WT,file = "Endo_Viability_Wilcox.csv")
write.csv(tmp_KW,file = "Endo_Viability_kruskal.csv")

####TKI_IC50
###Heatmap
rm(list = ls())
graphics.off()
load("FUSCC_SNF_RW_Validation_20250401.RData")

snf1 <- subset(RTK_2,RTK_2$SNF_byDP_Subtype=="SNF1")[,2:12]
snf2 <- subset(RTK_2,RTK_2$SNF_byDP_Subtype=="SNF2")[,2:12]
snf3 <- subset(RTK_2,RTK_2$SNF_byDP_Subtype=="SNF3")[,2:12]
snf4 <- subset(RTK_2,RTK_2$SNF_byDP_Subtype=="SNF4")[,2:12]

SNF1_mean <- colMeans(snf1,na.rm = T)
SNF2_mean <- colMeans(snf2,na.rm = T)
SNF3_mean <- colMeans(snf3,na.rm = T)
SNF4_mean <- colMeans(snf4,na.rm = T)

Data_IC50_mean <- data.frame(
  row.names = colnames(RTK_2)[2:12],
  SNF1=SNF1_mean,
  SNF2=SNF2_mean,
  SNF3=SNF3_mean,
  SNF4=SNF4_mean
)

pheatmap(Data_IC50_mean,scale = "row",cluster_rows = T,cluster_cols = F)
write.csv(Data_IC50_mean,file = "TKI_IC50_mean.csv") # Supplementary Table 6

###Stat
tmp_data_2 <- melt(RTK_2, "SNF_byDP_Subtype",variable_name = "IC50")
tmp_data_2 <- tmp_data_2[!is.na(tmp_data_2$value),]
tmp_WT <- compare_means(value~SNF_byDP_Subtype,tmp_data_2,group.by = "variable",method = "wilcox.test")
tmp_KW <- compare_means(value~SNF_byDP_Subtype,tmp_data_2,group.by = "variable",method = "kruskal.test")
write.csv(tmp_WT,file = "RTK_IC50_Wilcox.csv")
write.csv(tmp_KW,file = "RTK_IC50_kruskal.csv")