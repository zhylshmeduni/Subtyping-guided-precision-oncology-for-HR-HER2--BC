library(openxlsx)
library(survminer) 
library(survival) 

rm(list = ls())
graphics.off()

load("FUSCC_SNF_RW_Validation_20250401.RData")

rm(list = ls()[which(ls() != "matched_dataframe")])

fit2 <- survfit(Surv(DFS_Interval, DFS_Status) ~ Match_Group, data = matched_dataframe)
summary(fit2)
p<- ggsurvplot(fit2, data = PFS,
               surv.median.line = "hv",
               risk.table = TRUE,
               pval = TRUE,
               #add.all = TRUE,
               #conf.int = TRUE,
               palette = "npg",
               #xlim = c(0, 21), ylim = c(0, 1),
               ylab = "Progression-free survival"
) #Figure 6B
print(p)