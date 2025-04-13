library(openxlsx)
library(survminer) 
library(survival)
library(dplyr)
library(ggplot2)
library(tidyverse)

rm(list = ls())
graphics.off()

load("FUSCC_SNF_RW_Validation_20250401.RData")

rm(list = ls()[which(ls() != "RW_whole_info_advanced")])

clinical_dataframe <- RW_whole_info_advanced

dataframe_forplot <- as.data.frame(table(clinical_dataframe$SNF_Subtype))
colnames(dataframe_forplot) <- c("SNF_Subtype", "Count")

dataframe_forplot <- rbind(primary_percentage,
                           dataframe_forplot,
                           external_percentage)

dataframe_forplot$Cohort <- c(rep("Multi-omics Cohort", 4),
                              rep("Real-world Cohort", 4),
                              rep("Multi-center Cohort", 4))

dataframe_forplot <- plyr::ddply(dataframe_forplot, "Cohort", transform,
                                 Percentage = Count / sum(Count))
dataframe_forplot$SNF_Subtype <- factor(dataframe_forplot$SNF_Subtype,
                                        levels = c("SNF1", "SNF2", "SNF3", "SNF4"))
dataframe_forplot$Cohort <- factor(dataframe_forplot$Cohort,
                                   levels = c("Multi-omics Cohort", "Real-world Cohort", "Multi-center Cohort"))


p <- ggplot(data=dataframe_forplot, aes(x=Cohort,y=Percentage,fill=SNF_Subtype))+
  geom_bar(stat = "identity", position = "fill", width = 0.7)+
  scale_fill_manual(values = c("#3D76AE", "#53AD4A", "#EDAB3C", "#C3392A"))+
  scale_y_continuous(expand = expansion(mult = c(0.01,0.1)),
                     labels = scales::percent_format())+
  labs(x=NULL,y="Percentages of SNF Subtypes")+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5))