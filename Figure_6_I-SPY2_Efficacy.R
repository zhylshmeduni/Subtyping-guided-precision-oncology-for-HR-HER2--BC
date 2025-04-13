library(ggplot2)
library(openxlsx)

rm(list = ls())
graphics.off()

load("FUSCC_SNF_RW_Validation_20250401.RData")

other.pd1 <- I_SPY2_Data
other.pd1$Treatment_group[other.pd1$Treatment_group != "anti-PD1"&other.pd1$Treatment_group != "Ctr"] <- "others"
dataframe <- matrix(table(other.pd1$Treatment_group,other.pd1$pCR),nrow = 3, ncol = 2, byrow = FALSE,
                    dimnames = list(c( "anti-PD1" ,"Ctr","others"), c("nonpCR","pCR")))
result <- fisher.test(dataframe,simulate.p.value = TRUE)
p_value <- format(result$p.value, digits = 4)

row_percentage <- dataframe/rowSums(dataframe)
data_df <- as.data.frame(row_percentage)
data_df$Treatment<- rownames(data_df)
data_long <- tidyr::gather(data_df, key = "response", value = "percentage", -Treatment)
data_long[,1]<- factor(data_long[,1],levels=c("Ctr","others","anti-PD1"))  
data_long[,2]<- factor(data_long[,2],levels=c("nonpCR","pCR") )
p <- ggplot(data_long,aes(x=Treatment,y=percentage,fill=response))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = c("#c7e9c0","#53AD4A"))+
  geom_text(aes(label = paste0("p =", p_value)),
            x = 1.5, y = 1.02,  size = 4) #Figure 6A
