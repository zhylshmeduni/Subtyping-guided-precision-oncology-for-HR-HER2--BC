library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)

rm(list = ls())
graphics.off()

##Loading_Geneset
geneset.PDGF <- read.gmt("./gmt/REACTOME_SIGNALING_BY_PDGF.v2023.1.Hs.gmt")
geneset.PDGF <- geneset.PDGF[is.na(geneset.PDGF$gene)==F,]
geneset.VEGF <- read.gmt("./gmt/REACTOME_SIGNALING_BY_VEGF.v2023.1.Hs.gmt")
geneset.VEGF <- geneset.VEGF[is.na(geneset.VEGF$gene)==F,]

##PDGF-Sorafenib
geneList.SO <- read.csv("./GSEA/selected/BL-SO/analysis/BL_SO_REACTOME_PDGF.Gsea.1679760863824/1.csv")
geneList.SO <- geneList.SO[order(geneList.SO$SCORE, decreasing = T),]
geneList.SO_list <- geneList.SO$SCORE
names(geneList.SO_list) <- geneList.SO$Name
SO.PDGF <-GSEA(geneList.SO_list,TERM2GENE=geneset.PDGF,pvalueCutoff=1)

##PDGFR-Famitinib
geneList.FA <- read.csv("./GSEA/selected/BL-FA/analysis/BL_FA_REACTOME_PDGFR.Gsea.1679760494717/1.csv")
geneList.FA <- geneList.FA[order(geneList.FA$SCORE, decreasing = T),]
geneList.FA_list <- geneList.FA$SCORE
names(geneList.FA_list) <- geneList.FA$Name
FA.PDGF <-GSEA(geneList.FA_list,TERM2GENE=geneset.PDGF,pvalueCutoff=1)

####
pdf("./GSEA/selected/FA.PDGF_GSEA.pdf",width = 10,height = 6)
gseaplot2(FA.PDGF,1,color="red",pvalue_table = T,title="",base_size=10,ES_geom="line") ####Figure 4L
dev.off()
pdf("./GSEA/selected/SO.PDGF_GSEA.pdf",width = 8,height = 6)
gseaplot2(SO.PDGF,1,color="red",pvalue_table = T,title="",base_size=10,ES_geom="line") ####Figure 4M
dev.off()

#####
#gseaplot2(FA.PDGF,"REACTOME_SIGNALING_BY_PDGF",color="red",pvalue_table = T,title="",base_size=10,ES_geom="line") 