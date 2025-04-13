library(openxlsx)
library(survminer) 
library(survival) 
library(dplyr)

rm(list = ls())
graphics.off()

load("FUSCC_SNF_RW_Validation_20250401.RData")

rm(list = ls()[which(ls() != "RW_whole_info_advanced")])

color <- c("#3D76AE", "#53AD4A", "#EDAB3C", "#C3392A")

clinical_dataframe <- RW_whole_info_advanced
clinical_dataframe <- filter(clinical_dataframe, WhetherStageIV == "R")
clinical_dataframe <- select(clinical_dataframe, ID_short, SNF_Subtype,
                             PD_2, DFS_2)
clinical_dataframe <- filter(clinical_dataframe, DFS_2 != "/")
clinical_dataframe$DFS_2 <- as.numeric(clinical_dataframe$DFS_2)
clinical_dataframe$PD_2 <- as.numeric(clinical_dataframe$PD_2)

surv_object <- Surv(time = clinical_dataframe$DFS_2, event = clinical_dataframe$PD_2)
fit1 <- survfit(surv_object ~ SNF_Subtype, data = clinical_dataframe)

p <- ggsurvplot(fit1, data = clinical_dataframe, pval = TRUE, palette = color,
                legend = "right",
                legend.title = "SNF Subtype",
                legend.labs = c("SNF1","SNF2","SNF3","SNF4"),
                risk.table = TRUE)
