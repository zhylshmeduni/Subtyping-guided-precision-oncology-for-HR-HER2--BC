library(GSVA)
library(GSEABase)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pheatmap)

####Cell_Cycle_Scores
rm(list=ls())
graphics.off()
load("CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20221010.Rdata")

rm(list = ls()[which(ls() != "SNF_Cluster" & ls() != "exp.fpkm.TT")])

pathway = getGmt("./genesets_Cell_Cycle.gmt")

gsva_matrix1 <- gsva(as.matrix(exp.fpkm.TT), pathway, method='ssgsea')

rownames(gsva_matrix1) <- c("Cell_Cycle_Scores")
gsva_matrix1 <- gsva_matrix1[,names(SNF_Cluster)]
gsva_matrix1 <- rbind(gsva_matrix1, SNF_Cluster)
gsva_matrix1 <- as.data.frame(t(gsva_matrix1))
gsva_matrix1 <- mutate(gsva_matrix1, SNF_Subtype = paste0("SNF", SNF_Cluster))
colnames(gsva_matrix1)[1] <- c("Cell_Cycle_Scores")

comparisons <- list(
  c("SNF1", "SNF2"),
  c("SNF1", "SNF3"),
  c("SNF1", "SNF4"),
  c("SNF2", "SNF3"),
  c("SNF2", "SNF4"),
  c("SNF3", "SNF4")
)

ggviolin(gsva_matrix1, x = "SNF_Subtype", y = "Cell_Cycle_Scores", fill = "SNF_Subtype",
         palette = c("steelblue", "forestgreen", "orange", "firebrick"),
         add = "boxplot",
         add.params = list(fill = "white", width = 0.2)) +
  stat_compare_means(comparisons = comparisons, method = "t.test") +  # 两两比较
  stat_compare_means(method = "kruskal.test") +  # 整体检验
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 12, color = c("steelblue", "forestgreen", "orange", "firebrick")),
    legend.position = "none"
  ) +
  labs(
    y = "Cell_Cycle_Scores",
    x = NULL,
    title = "G1/S"
  )

####HRD_Scores
rm(list=ls())
graphics.off()
load("CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20221010.Rdata")

rm(list = ls()[which(ls() != "SNF_Cluster" & ls() != "luminal" &
                     ls() != "germline" & ls() != "maf.df")])

luminal <- filter(luminal, HRD != "NA")
luminal$HRD <- as.numeric(luminal$HRD)
luminal <- luminal[,c("PatientCode", "HRD")]

maf.df <- filter(maf.df, Hugo_Symbol %in% c("BRCA1", "BRCA2"))

luminal$WhetherSomatic <- luminal$PatientCode %in% maf.df$Tumor_Sample_Barcode
luminal$WhetherGermline <- luminal$PatientCode %in% germline$PatientCode
luminal <- mutate(luminal, BRCA = ifelse(WhetherSomatic == TRUE | WhetherGermline == TRUE, "BRCA1/2-MT", "BRCA1/2-WT"))

SNF_Cluster <- SNF_Cluster[intersect(names(SNF_Cluster), luminal$PatientCode)]

luminal <- luminal[names(SNF_Cluster),]
luminal$SNF_Subtype <- SNF_Cluster
luminal <- mutate(luminal, SNF_Subtype = paste0("SNF", SNF_Subtype))

df <- luminal
df$group_label <- paste(df$BRCA, ifelse(df$SNF_Subtype == "SNF3", "(SNF3)", "(non-SNF3)"))
df <- mutate(df, SNF_Subtype = ifelse(SNF_Subtype == "SNF3", "SNF3", "non-SNF3"))

df$group_label <- factor(df$group_label, levels = c(
  "BRCA1/2-MT (non-SNF3)", "BRCA1/2-MT (SNF3)", 
  "BRCA1/2-WT (non-SNF3)", "BRCA1/2-WT (SNF3)"
))

fill_colors <- c("non-SNF3" = "gray70", "SNF3" = "#f5a623")

p <- ggplot(df, aes(x = group_label, y = HRD, fill = SNF_Subtype)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = fill_colors) +
  labs(x = "", y = "HRD", title = "HRD score") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  ) +
  geom_signif(comparisons = list(
    c("BRCA1/2-MT (non-SNF3)", "BRCA1/2-MT (SNF3)"),
    c("BRCA1/2-WT (non-SNF3)", "BRCA1/2-WT (SNF3)"),
    c("BRCA1/2-MT (non-SNF3)", "BRCA1/2-WT (non-SNF3)"),
    c("BRCA1/2-MT (SNF3)", "BRCA1/2-WT (SNF3)")
  ),
  step_increase = 0.12,
  map_signif_level = TRUE,
  tip_length = 0.01)

print(p)
