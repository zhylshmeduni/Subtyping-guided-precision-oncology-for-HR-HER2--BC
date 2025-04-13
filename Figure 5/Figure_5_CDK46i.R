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

color_ORR <- c("#B24745", "gray")
names(color_ORR) <- c("Yes", "No")

clinical_dataframe <- clinical_dataframe %>%
  select(HID, SNF_Subtype, BestResponse_3) %>%
  filter(BestResponse_3 != "/") %>%
  mutate(ORR = ifelse(BestResponse_3 == "PR", "Yes", "No"))

dataframe_forplot <- as.data.frame(table(clinical_dataframe$SNF_Subtype,
                                         clinical_dataframe$ORR))
colnames(dataframe_forplot) <- c("SNF_Subtype", "ORR", "Count")
####Stat
dataframe_fortest <- as.data.frame(matrix(dataframe_forplot$Count, nrow = 3))
set.seed(2022)
test <- fisher.test(dataframe_fortest)
p_value <- test$p.value
####
dataframe_forplot <- plyr::ddply(dataframe_forplot, "SNF_Subtype", transform,
                                 Percentage = Count / sum(Count))
dataframe_forplot$SNF_Subtype <- factor(dataframe_forplot$SNF_Subtype,
                                        levels = c("SNF1", "SNF2", "SNF3", "SNF4"))
dataframe_forplot$ORR <- factor(dataframe_forplot$ORR,
                                levels = c("No", "Yes"))
####BarPlot
p <- ggplot(data = dataframe_forplot, aes(x=SNF_Subtype, y=Percentage, fill=ORR))+
  geom_bar(stat = "identity", position = "fill", width = 0.75)+
  scale_fill_manual(values = color_ORR)+
  scale_y_continuous(expand = expansion(mult = c(0.01,0.1)),
                     labels = scales::percent_format())+
  labs(x=NULL,y="Response to CDK4/6")+
  ggtitle(paste0("p = ", round(p_value, 3)))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(),
        legend.position = "right")
