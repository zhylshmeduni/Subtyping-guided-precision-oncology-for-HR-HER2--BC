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

load("CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20221010.Rdata")
rm(list=ls()[which(ls() != "color")])

color_SNF <- color
names(color_SNF) <- c("SNF1", "SNF2", "SNF3", "SNF4")
rm(color)

load("FUSCC_SNF_RW_Validation_20250401.RData")

wanted_somatic_panel <- c("BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2")
wanted_germline_panel <- c("BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2")

rownames(clinical_dataframe) <- clinical_dataframe$Analysis_ID
Luminal_WES_dataframe <- Luminal_WES_dataframe[intersect(wanted_somatic_panel, rownames(Luminal_WES_dataframe)),
                                               intersect(colnames(Luminal_WES_dataframe), clinical_dataframe$Analysis_ID)]

Luminal_WES_dataframe[Luminal_WES_dataframe == "Frame_Shift_Del"] <- "Frameshift"
Luminal_WES_dataframe[Luminal_WES_dataframe == "Frame_Shift_Ins"] <- "Frameshift"
Luminal_WES_dataframe[Luminal_WES_dataframe == "In_Frame_Del"] <- "Inframe"
Luminal_WES_dataframe[Luminal_WES_dataframe == "In_Frame_Ins"] <- "Inframe"
Luminal_WES_dataframe[Luminal_WES_dataframe == "Missense_Mutation"] <- "Missense"
Luminal_WES_dataframe[Luminal_WES_dataframe == "Nonsense_Mutation"] <- "Nonsense"
Luminal_WES_dataframe[Luminal_WES_dataframe == "Nonstop_Mutation"] <- "Nonstop"
Luminal_WES_dataframe[Luminal_WES_dataframe == "Splice_Site"] <- "Splicing"
Luminal_WES_dataframe[Luminal_WES_dataframe == "Multi_Hit"] <- "Multi_Hit"

nondected <- setdiff(clinical_dataframe$Analysis_ID, colnames(Luminal_WES_dataframe))
supplementary <- as.data.frame(matrix("", nrow = nrow(Luminal_WES_dataframe), ncol = length(nondected)))
colnames(supplementary) <- nondected
rownames(supplementary) <- rownames(Luminal_WES_dataframe)
Luminal_WES_dataframe <- cbind(Luminal_WES_dataframe, supplementary)

clinical_dataframe <- clinical_dataframe[colnames(Luminal_WES_dataframe),]

####Raw_Germline_Data
raw_dataframe <- filter(raw_dataframe, Tumor_Sample_Barcode %in% colnames(Luminal_WES_dataframe))
Luminal_Germ_dataframe <- as.data.frame(matrix("", nrow = 7, ncol = 60))
rownames(Luminal_Germ_dataframe) <- c("BRCA1", "PALB2", "BRCA2_Germline", "BRCA1_Germline", "PALB2_Germline", "ATM_Germline", "CHEK2_Germline")
colnames(Luminal_Germ_dataframe) <- colnames(Luminal_WES_dataframe)

for (i in 1:nrow(raw_dataframe)){
  Luminal_Germ_dataframe[paste0(raw_dataframe$Hugo_Symbol[i], "_Germline"), raw_dataframe$Tumor_Sample_Barcode[i]] <- raw_dataframe$Variant_Classification[i]
}

Luminal_Germ_dataframe[Luminal_Germ_dataframe == "Missense_Mutation"] <- "Missense"
Luminal_Germ_dataframe[Luminal_Germ_dataframe == "Nonsense_Mutation"] <- "Nonsense"
Luminal_WES_dataframe <- rbind(Luminal_WES_dataframe, Luminal_Germ_dataframe)

gene_category <- c(rep("HRR",5), rep("Germline",5))
gene_category <- factor(gene_category, levels = c("HRR", "Germline"))

###Heatmap legend
col <- c("Splicing" = "#6b4694", "Missense" = "#4d93cd", "Nonsense" = "#e63a3a",
         "Frameshift" = '#fab71b', "Inframe" = "#ecd71e", "Nonstop" = "#ab5b9e",
         "Translation_Start_Site" = '#8BBF8D', "Multi_Hit" = "#018b38")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h, 
              gp = gpar(fill = "#eaeaea", col = "white"))},
  Splicing = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"),
              gp = gpar(fill = col["Splicing"], col = NA))
  },
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"),
              gp = gpar(fill = col["Missense"], col = NA))
  },
  Nonsense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"),
              gp = gpar(fill = col["Nonsense"], col = NA))
  },
  Frameshift = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"),
              gp = gpar(fill = col["Frameshift"], col = NA))
  },
  Inframe = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"),
              gp = gpar(fill = col["Inframe"], col = NA))
  },
  Nonstop = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"),
              gp = gpar(fill = col["Nonstop"], col = NA))
  },
  Translation_Start_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"),
              gp = gpar(fill = col["Translation_Start_Site"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"),
              gp = gpar(fill = col["Multi_Hit"], col = NA))})

heatmap_legend <- list(title = "SNV", at = c("Splicing", "Missense", "Nonsense", "Frameshift", "Inframe", "Nonstop", "Translation_Start_Site", "Multi_Hit"),
                       labels = c("Splicing", "Missense", "Nonsense", "Frameshift", "Inframe", "Nonstop", "Translation_Start_Site", "Multi_Hit"))

####Stat
dataframe_fortest <- Luminal_WES_dataframe
dataframe_fortest[dataframe_fortest == ""] <- "WT"
dataframe_fortest[dataframe_fortest != "WT"] <- "Mut"
dataframe_fortest <- rbind(dataframe_fortest, clinical_dataframe$SNF_Subtype)
rownames(dataframe_fortest)[nrow(dataframe_fortest)] <- "SNF_Subtype"
dataframe_fortest <- as.data.frame(t(dataframe_fortest))
dataframe_fortest <- mutate(dataframe_fortest, SNF_Subtype = ifelse(SNF_Subtype == "SNF3", "SNF3", "non-SNF3"))

Significance <- c()
set.seed(2022)

for (j in 1:nrow(Luminal_WES_dataframe)){
  list_fortest <- as.data.frame(table(dataframe_fortest[, rownames(Luminal_WES_dataframe)[j]],
                                      dataframe_fortest[, "SNF_Subtype"]))
  colnames(list_fortest) <- c("SNV", "Group", "Count")
  list_fortest <- matrix(list_fortest$Count, nrow = 2)
  p_value <- chisq.test(list_fortest)$p.value
  
  if (p_value >= 0.05){
    Significance <- c(Significance, "ns")
  } else if (p_value < 0.05 & p_value >= 0.01){
    Significance <- c(Significance, "*")
  } else if (p_value < 0.01 & p_value >= 0.001){
    Significance <- c(Significance, "**")
  } else if (p_value < 0.001){Significance <- c(Significance, "***")}
}

oncoPrint(Luminal_WES_dataframe[,clinical_dataframe$Analysis_ID[clinical_dataframe$SNF_Subtype == "SNF3"]],
          alter_fun = alter_fun, alter_fun_is_vectorized = TRUE, col = col, show_pct = TRUE,
          top_annotation = HeatmapAnnotation(SNF_Subtype = clinical_dataframe$SNF_Subtype[which(clinical_dataframe$SNF_Subtype == "SNF3")],
                                             annotation_name_side = "left",
                                             annotation_name_gp = gpar(fontsize = 10),
                                             col = list(SNF_Subtype = color_SNF,
                                                        Olaparib_IC50 = pal_H),
                                             gp = gpar(col = "black")),
          left_annotation = rowAnnotation(Significance = anno_text(Significance, just = "center")),
          right_annotation = NULL,
          row_split = gene_category,
          remove_empty_rows = FALSE,
          remove_empty_columns = FALSE,
          row_names_side = "left",
          pct_side = "right")+
  oncoPrint(Luminal_WES_dataframe[,clinical_dataframe$Analysis_ID[clinical_dataframe$SNF_Subtype == "SNF1"]],
            alter_fun = alter_fun, alter_fun_is_vectorized = TRUE, col = col, show_pct = TRUE,
            top_annotation = HeatmapAnnotation(SNF_Subtype = clinical_dataframe$SNF_Subtype[which(clinical_dataframe$SNF_Subtype == "SNF1")],
                                               annotation_name_side = "left",
                                               annotation_name_gp = gpar(fontsize = 10),
                                               col = list(SNF_Subtype = color_SNF,
                                                          Olaparib_IC50 = pal_H),
                                               gp = gpar(col = "black")),
            right_annotation = NULL,
            row_split = gene_category,
            remove_empty_rows = FALSE,
            remove_empty_columns = FALSE,
            show_row_names = FALSE,
            row_names_side = "left",
            pct_side = "right")+
  oncoPrint(Luminal_WES_dataframe[,clinical_dataframe$Analysis_ID[clinical_dataframe$SNF_Subtype == "SNF2"]],
            alter_fun = alter_fun, alter_fun_is_vectorized = TRUE, col = col, show_pct = TRUE,
            top_annotation = HeatmapAnnotation(SNF_Subtype = clinical_dataframe$SNF_Subtype[which(clinical_dataframe$SNF_Subtype == "SNF2")],
                                               annotation_name_side = "left",
                                               annotation_name_gp = gpar(fontsize = 10),
                                               col = list(SNF_Subtype = color_SNF,
                                                          Olaparib_IC50 = pal_H),
                                               gp = gpar(col = "black")),
            right_annotation = NULL,
            row_split = gene_category,
            remove_empty_rows = FALSE,
            remove_empty_columns = FALSE,
            show_row_names = FALSE,
            row_names_side = "left",
            pct_side = "right")+
  oncoPrint(Luminal_WES_dataframe[,clinical_dataframe$Analysis_ID[clinical_dataframe$SNF_Subtype == "SNF4"]],
            alter_fun = alter_fun, alter_fun_is_vectorized = TRUE, col = col, show_pct = TRUE,
            top_annotation = HeatmapAnnotation(SNF_Subtype = clinical_dataframe$SNF_Subtype[which(clinical_dataframe$SNF_Subtype == "SNF4")],
                                               annotation_name_side = "left",
                                               annotation_name_gp = gpar(fontsize = 10),
                                               col = list(SNF_Subtype = color_SNF,
                                                          Olaparib_IC50 = pal_H),
                                               gp = gpar(col = "black")),
            right_annotation = NULL,
            row_split = gene_category,
            remove_empty_rows = FALSE,
            remove_empty_columns = FALSE,
            show_row_names = FALSE,
            row_names_side = "left",
            pct_side = "right")