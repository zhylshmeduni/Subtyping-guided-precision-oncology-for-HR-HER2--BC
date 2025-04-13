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

####Percentage of SNF Subtypes
clinical_dataframe <- RW_whole_info_early_stage

dataframe_forplot <- as.data.frame(table(clinical_dataframe$SNF_Subtype))
colnames(dataframe_forplot) <- c("SNF_Subtype", "Count")
dataframe_forplot <- rbind(dataframe_forplot, primary_percentage)
dataframe_forplot$Cohort <- c(rep("Real-world", 4),
                              rep("Multi-omics", 4))

dataframe_forplot <- plyr::ddply(dataframe_forplot, "Cohort", transform,
                                 Percentage = Count / sum(Count))
dataframe_forplot$SNF_Subtype <- factor(dataframe_forplot$SNF_Subtype,
                                        levels = c("SNF1", "SNF2", "SNF3", "SNF4"))
dataframe_forplot$Cohort <- factor(dataframe_forplot$Cohort,
                                   levels = c("Multi-omics", "Real-world"))

dataframe_fortest <- data.frame(
  Multi = dataframe_forplot$Count[1:4],
  RW = dataframe_forplot$Count[5:8]
)
rownames(dataframe_fortest) <- dataframe_forplot$SNF_Subtype[1:4]
set.seed(2022)
test <- fisher.test(dataframe_fortest)
p_value <- test$p.value

p <- ggplot(data=dataframe_forplot, aes(x=Cohort,y=Percentage,fill=SNF_Subtype))+
  geom_bar(stat = "identity", position = "fill", width = 0.7)+
  scale_fill_manual(values = c("#3D76AE", "#53AD4A", "#EDAB3C", "#C3392A"))+
  scale_y_continuous(expand = expansion(mult = c(0.01,0.1)),
                     labels = scales::percent_format())+
  labs(x=NULL,y="Percentages of SNF Subtypes")+
  ggtitle(paste0("p = ",as.character(round(p_value,2))))+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5))
print(p) #Supplementary Figure 8D

####Proportions of ER-Positive Cells
clinical_dataframe <- RW_whole_info_early_stage

clinical_dataframe <- clinical_dataframe[,c("SNF_Subtype", "ER_Percentage")]
clinical_dataframe <- clinical_dataframe %>%
  filter(ER_Percentage != "NA") %>%
  mutate(ER_Level = ifelse(ER_Percentage <= 80, "≤80%", ">80%"))

dataframe_forplot <- as.data.frame(table(clinical_dataframe$SNF_Subtype, clinical_dataframe$ER_Level))
colnames(dataframe_forplot) <- c("SNF_Subtype", "ER_Percentage", "Count")
dataframe_forplot <- plyr::ddply(dataframe_forplot, "SNF_Subtype", transform, Percentage = Count / sum(Count))

dataframe_forplot$SNF_Subtype <- factor(dataframe_forplot$SNF_Subtype, levels = c("SNF1", "SNF2", "SNF3", "SNF4"))
dataframe_forplot$ER_Percentage <- factor(dataframe_forplot$ER_Percentage, levels = c(">80%", "≤80%"))

dataframe_fortest <- data.frame(t(matrix(dataframe_forplot$Count, nrow = 2)))
test <- fisher.test(dataframe_fortest)
p_value <- test$p.value

p <- ggplot(data=dataframe_forplot, aes(x=SNF_Subtype, y=Percentage, fill=ER_Percentage))+
  geom_bar(stat = "identity", position = "fill", width = 0.85)+
  scale_fill_manual(values = c("≤80%" = "#D7EFF4", ">80%" = "#4B90BB"))+
  scale_y_continuous(expand = expansion(mult = c(0.01,0.1)),
                     labels = scales::percent_format())+
  labs(x=NULL,y="Percentage")+
  ggtitle(paste0("p = ", as.character(round(p_value, 4))))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill=NA, colour = "black", size=2),
        legend.position = "bottom",
        axis.text.x=element_text(size = 16),
        axis.text.y=element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) #Supplementary Figure 8B

####Proportions of Ki67-Positive Cells
clinical_dataframe <- RW_whole_info_early_stage

clinical_dataframe <- clinical_dataframe[,c("SNF_Subtype", "Ki67_Status")]
clinical_dataframe <- clinical_dataframe %>%
  filter(Ki67_Status != "NA") %>%
  mutate(Ki67_Level = ifelse(Ki67_Status >= 30, "≥30%", "<30%")) %>%
  mutate(Ki67_Level = ifelse(Ki67_Status < 30 & Ki67_Status >= 15, "≥15%, <30%", Ki67_Level)) %>%
  mutate(Ki67_Level = ifelse(Ki67_Status < 15, "<15%", Ki67_Level))

dataframe_forplot <- as.data.frame(table(clinical_dataframe$SNF_Subtype, clinical_dataframe$Ki67_Level))
colnames(dataframe_forplot) <- c("SNF_Subtype", "Ki67_Status", "Count")
dataframe_forplot <- plyr::ddply(dataframe_forplot, "SNF_Subtype", transform, Percentage = Count / sum(Count))

dataframe_forplot$SNF_Subtype <- factor(dataframe_forplot$SNF_Subtype, levels = c("SNF1", "SNF2", "SNF3", "SNF4"))
dataframe_forplot$Ki67_Status <- factor(dataframe_forplot$Ki67_Status, levels = c("<15%", "≥15%, <30%", "≥30%"))

dataframe_fortest <- data.frame(t(matrix(dataframe_forplot$Count, nrow = 3)))
test <- chisq.test(dataframe_fortest)
p_value <- test$p.value

p <- ggplot(data=dataframe_forplot, aes(x=SNF_Subtype, y=Percentage, fill=Ki67_Status))+
  geom_bar(stat = "identity", position = "fill", width = 0.85)+
  scale_fill_manual(values = c("#d9f0a3", "#66bd63", "#1a9850"))+
  scale_y_continuous(expand = expansion(mult = c(0.01,0.1)),
                     labels = scales::percent_format())+
  labs(x=NULL,y="Percentage")+
  ggtitle(paste0("p = ", as.character(round(p_value, 4))))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill=NA, colour = "black", size=2),
        legend.position = "bottom",
        axis.text.x=element_text(size = 16),
        axis.text.y=element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) #Supplementary Figure 8C

####T Stage
clinical_dataframe <- RW_whole_info_early_stage

clinical_dataframe <- select(clinical_dataframe, SNF_Subtype, pT)
clinical_dataframe <- na.omit(clinical_dataframe)

dataframe_forplot <- as.data.frame(table(clinical_dataframe$SNF_Subtype,
                                         clinical_dataframe$pT))
colnames(dataframe_forplot) <- c("SNF_Subtype", "pT", "Count")
dataframe_forplot <- plyr::ddply(dataframe_forplot, "SNF_Subtype", transform,
                                 Percentage = Count / sum(Count))

dataframe_forplot$SNF_Subtype <- factor(dataframe_forplot$SNF_Subtype,
                                        levels = c("SNF1", "SNF2", "SNF3", "SNF4"))

#dataframe_forplot$pT <- factor(dataframe_forplot$pT,
#levels = c("pT0", "pT1", "pT2", "pT3"))
dataframe_forplot$pT <- factor(dataframe_forplot$pT,
                               levels = c("pT3", "pT2", "pT1", "pT0"))

dataframe_fortest <- data.frame(t(matrix(dataframe_forplot$Count, nrow = 4)))
set.seed(2022)
test <- fisher.test(dataframe_fortest)
p_value <- test$p.value

p <- ggplot(data=dataframe_forplot, aes(x=SNF_Subtype, y=Percentage, fill=pT))+
  geom_bar(stat = "identity", position = "fill", width = 0.8)+
  scale_fill_manual(values = c("#FF0000", "#FF6600", "#3399FF", "#0066CC"))+
  #scale_fill_manual(values = c("#0066CC", "#3399FF", "#FF6600", "#FF0000"))+
  scale_y_continuous(expand = expansion(mult = c(0.01,0.1)),
                     labels = scales::percent_format())+
  labs(x=NULL,y="Percentage")+
  ggtitle(paste0("p = ", as.character(round(p_value, 4))))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        panel.border = element_rect(fill=NA, colour = "black", size=2),
        axis.ticks = element_line(size = 2),
        axis.text.x=element_text(size = 16),
        axis.text.y=element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) #Supplementary Figure 8D

####N Stage
clinical_dataframe <- RW_whole_info_early_stage

clinical_dataframe <- select(clinical_dataframe, SNF_Subtype, pN)
clinical_dataframe <- na.omit(clinical_dataframe)

dataframe_forplot <- as.data.frame(table(clinical_dataframe$SNF_Subtype,
                                         clinical_dataframe$pN))
colnames(dataframe_forplot) <- c("SNF_Subtype", "pN", "Count")
dataframe_forplot <- plyr::ddply(dataframe_forplot, "SNF_Subtype", transform,
                                 Percentage = Count / sum(Count))

dataframe_forplot$SNF_Subtype <- factor(dataframe_forplot$SNF_Subtype,
                                        levels = c("SNF1", "SNF2", "SNF3", "SNF4"))
#dataframe_forplot$pN <- factor(dataframe_forplot$pN,
#levels = c("pN0", "pN1", "pN2", "pN3"))

dataframe_forplot$pN <- factor(dataframe_forplot$pN,
                               levels = c("pN3", "pN2", "pN1", "pN0"))

dataframe_fortest <- data.frame(t(matrix(dataframe_forplot$Count, nrow = 4)))
test <- chisq.test(dataframe_fortest)
p_value <- test$p.value

p <- ggplot(data=dataframe_forplot, aes(x=SNF_Subtype, y=Percentage, fill=pN))+
  geom_bar(stat = "identity", position = "fill", width = 0.8)+
  scale_fill_manual(values = c("#FF0000", "#FF6600", "#3399FF", "#0066CC"))+
  #scale_fill_manual(values = c("#0066CC", "#3399FF", "#FF6600", "#FF0000"))+
  scale_y_continuous(expand = expansion(mult = c(0.01,0.1)),
                     labels = scales::percent_format())+
  labs(x=NULL,y="Percentage")+
  ggtitle(paste0("p = ", as.character(round(p_value, 4))))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        panel.border = element_rect(fill=NA, colour = "black", size=2),
        axis.ticks = element_line(size = 2),
        axis.text.x=element_text(size = 16),
        axis.text.y=element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) #Supplementary Figure 8E
