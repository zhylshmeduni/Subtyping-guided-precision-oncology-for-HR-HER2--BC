library(GSVA)
library(GSEABase)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pheatmap)

####RTK_panel
rm(list = ls())
graphics.off()
load("CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20221010.Rdata")

rm(list = ls()[which(ls() != "SNF_Cluster" & ls() != "exp.fpkm.TT")])

RTK_GENEsymbol <- c(
  "EGFR", "ERBB3", "ERBB4", "KDR", "FLT1", "FLT4", "PDGFRA", "PDGFRB", 
  "FGFR1", "FGFR2", "FGFR3", "FGFR4", "MET", "FLT3", "KIT", "SRC", 
  "ERBB2", "ABL1", "MAP2K2", "MAP2K1", "BTK", "RET", "NTRK1", "NTRK2", 
  "NTRK3", "IGF1R", "AXL", "ALK", "SYK", "CSF1R", "PTK2", "MERTK", 
  "ROS1", "TYRO3", "RON", "TEK", "PTK2B", 
  "DYRK1A", "DYRK1B", "DYRK2", "DYRK3", "DYRK4", 
  "EPHA1", "EPHA2", "EPHA3", "EPHA4", "EPHA5", "EPHA6", "EPHA7", 
  "EPHA8", "EPHA10", "EPHB1", "EPHB2", "EPHB3", "EPHB4", "EPHB6", 
  "TNK2", "DDR1", "DDR2", "FER"
)

pval_to_signif <- function(p) {
  if (is.na(p)) return("")
  else if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}

expr_rtk <- exp.fpkm.TT[rownames(exp.fpkm.TT) %in% RTK_GENEsymbol, ]
rm(exp.fpkm.TT)
expr_rtk <- expr_rtk[,names(SNF_Cluster)]
expr_rtk_log <- log2(expr_rtk + 1)
expr_rtk <- t(scale(t(expr_rtk)))

pvals <- apply(expr_rtk, 1, function(gene_expr) {
  df <- data.frame(expr = gene_expr, group = SNF_Cluster)
  kruskal.test(expr ~ group, data = df)$p.value
})

signif_labels <- sapply(pvals, pval_to_signif)

new_rownames <- paste0(rownames(expr_rtk), " ", signif_labels)

group_levels <- levels(as.factor(SNF_Cluster))
avg_expr_by_group <- sapply(group_levels, function(g) {
  rowMeans(expr_rtk[, SNF_Cluster == g, drop = FALSE])
})
avg_expr_by_group <- t(avg_expr_by_group)
rownames(avg_expr_by_group) <- c("SNF1", "SNF2", "SNF3", "SNF4")

labels_matrix <- round(avg_expr_by_group, 2)

avg_expr_rotated <- avg_expr_by_group
avg_expr_rotated <- t(avg_expr_rotated)

labels_rotated <- t(labels_matrix)

pheatmap(
  mat = avg_expr_rotated,
  display_numbers = labels_rotated,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_number = 10,
  main = "RTK Gene Expression",
  angle_col = 0
)
