library(GSVA)
library(GSEABase)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pheatmap)

####SET_Index
rm(list = ls())
graphics.off()
load("CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20221010.Rdata")

rm(list = ls()[which(ls() != "SNF_Cluster" & ls() != "exp.fpkm.TT")])

set_genes <- c(
  "SLC39A6", "ESR1", "NAT1", "AZGP1", "SCUBE2", "DNAJC12",
  "QDPR", "STC2", "MRPS30", "MAPT", "CA12", "CD3D",
  "CD2", "NPY1R", "ABAT", "ADCY1", "PDZK1", "KCNE4"
)

control_genes <- c(
  "LDHA", "ATP5J2", "VDAC2", "DARS", "UGP2",
  "UBE2Z", "AK2", "WIPF2", "APPBP2", "TRIM2"
)

exp.fpkm.TT <- exp.fpkm.TT[c(set_genes, control_genes),]
exp.fpkm.TT <- as.data.frame(t(exp.fpkm.TT))
exp.fpkm.TT <- exp.fpkm.TT[names(SNF_Cluster),]
exp.fpkm.TT$SET_Index <- (rowSums(exp.fpkm.TT[, 1:18]) / 18) - ((rowSums(exp.fpkm.TT[, (ncol(exp.fpkm.TT)-9):ncol(exp.fpkm.TT)]) / 10) + 2)
exp.fpkm.TT$SNF_Subtype <- SNF_Cluster
exp.fpkm.TT <- mutate(exp.fpkm.TT, SNF_Subtype = paste0("SNF", SNF_Subtype))
exp.fpkm.TT$SNF_Subtype <- factor(exp.fpkm.TT$SNF_Subtype, levels = c("SNF1", "SNF2", "SNF3", "SNF4"))

comparisons <- list(
  c("SNF1", "SNF2"),
  c("SNF1", "SNF3"),
  c("SNF1", "SNF4"),
  c("SNF2", "SNF3"),
  c("SNF2", "SNF4"),
  c("SNF3", "SNF4")
)

ggviolin(exp.fpkm.TT, x = "SNF_Subtype", y = "SET_Index", fill = "SNF_Subtype",
         palette = c("steelblue", "forestgreen", "orange", "firebrick"),
         add = "boxplot",
         add.params = list(fill = "white", width = 0.2)) +
  stat_compare_means(comparisons = comparisons, method = "t.test") +
  stat_compare_means(method = "kruskal.test", label.y = 210) +
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
    y = "Sensitivity to endocrine therapy\n(SET index)",
    x = NULL,
    title = "SET index across SNF clusters"
  ) #Figrue 4B

####Endo_Sen&Resist_Heatmap
rm(list = ls())
graphics.off()
load("CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20221010.Rdata")

pathway = getGmt("./genesets_endo_sen&resist.gmt")

gsva_matrix1 <- gsva(as.matrix(exp.fpkm.TT), pathway, method='ssgsea')  
rownames(gsva_matrix1) = c("Aromatase_inhibitor_responsive_module",
                           "Tamoxifen_or_SERD_responsive_module",
                           "RTK_pathway","EMT_pathway","WNT_pathway","CellCycle_pathway","PI3K_pathway")

annot_color = list(SNF = c("SNF1" = color[1],"SNF2" = color[2],"SNF3" = color[3],"SNF4" = color[4]))
annot_col = data.frame(row.names = colnames(gsva_matrix1),
                       SNF = paste0("SNF",SNF_Cluster[colnames(gsva_matrix1)]) )

pheatmap(gsva_matrix1[,names(sort(SNF_Cluster)) ],cluster_cols = F,cluster_rows = F,
         scale = "row",
         breaks = unique(c(seq(-1.6,1.6,length = 100))),
         annotation_col = annot_col,
         annotation_colors = annot_color,
         show_colnames = F,
         treeheight_row = 0,
         filename = "drug_pathway_SNF.pdf",
         height = 3,width = 10) #Figure 4C