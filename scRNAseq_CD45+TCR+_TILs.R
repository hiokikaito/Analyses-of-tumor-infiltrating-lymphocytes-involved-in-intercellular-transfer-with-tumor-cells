#################################################
# Name: scRNAseq_CD45+TCR+_TILs.R
# 
# Description: R codes used for the analysis of scRNAseq data on all CD45+ TCRb+ sorted Tumor Infiltrating Lymphocytes collected from B16 ZsGreen+ tumors.
#
# Input: scRNAseq data in 10x Genomics file format.
#
# Date updated: Nov. 10, 2023
# 
# Author: Kaito Hioki 
# Contact: khioki@umass.edu
#################################################

#install.packages("dplyr")
library(dplyr)
#install.packages('Seurat')
library(Seurat)
#install.packages('patchwork')
library(patchwork)

#install.packages("remotes")
library(remotes)

library(MAST)
# BiocManager::install("DESeq2")
library(DESeq2)
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

# remotes::install_github("carmonalab/UCell", ref="v1.3")
library(UCell)
# remotes::install_github("carmonalab/scGate")
# library(scGate)
# remotes::install_github("carmonalab/ProjecTILs")
library(ProjecTILs)
# # BiocManager::install("dittoSeq")
# library(dittoSeq)

library(tidyr)
library(purrr)
library(magrittr)
library(ggplot2)
library(tidyseurat)


#Load the ProjecTILs reference atlas
ref <- load.reference.map("ref_TILAtlas_mouse_v1.rds")

Idents(ref, cells = WhichCells(ref, idents = c("Tfh"))) <- "Other"
ref$functional.cluster <- Idents(ref)

celltypes <- c("Other","Treg","Th1","CD4_NaiveLike","CD8_NaiveLike","CD8_EarlyActiv","CD8_EffectorMemory","CD8_Tpex","CD8_Tex")
for (i in celltypes){
  Idents(ref, cells = WhichCells(ref, idents = i)) <- i
}

#Open the reference atlas
refCols <- c("CD8_Tex"="#edbe2a","CD8_Tpex"="#A58AFF","CD8_EffectorMemory"="#53B400","CD8_EarlyActiv"="#F8766D","CD8_NaiveLike"="#00B6EB","CD4_NaiveLike"="#b57007","Th1"="#87f6a5","Treg"="#e812dd","Other"="gray")
DimPlot(ref,label = F, cols = refCols)
DimPlot(ref,label = T, cols = refCols)


#expression distribution of marker genes across reference subtypes
markers <- c("Cd4","Cd8a","Ccr7","Tcf7","Pdcd1","Havcr2","Tox","Izumo1r","Cxcr6","Xcl1","Gzmb","Gzmk","Ifng","Foxp3")
VlnPlot(ref, features = markers, stack = T, flip = T, fill.by = "ident", cols=refCols, assay = "RNA")

DotPlot(ref, features = markers, assay = "RNA", cols = c("gray", "red")) + RotatedAxis()







########################
#Load sample data
########################
test <- read.sc.query("datafiles/", type = "10x", min.cells = 3, min.features = 50)



######Assign group by barcode. -1 is ZsHigh, -2 is ZsLow.
sub <- (colnames(test))
df <- data.frame(sub)
df$group <- ("")
for (row in 1:nrow(df)) {
  if ( endsWith(df[row,1], "-1") ) {
    df[row,2] = 1
  }
  else{
    df[row,2] = 2
  }
}
test@meta.data$orig.ident <- df[,2]
df <- NULL
table(test$orig.ident)



########################
#Quality control
########################
summary(test$nCount_RNA)
summary(test$nFeature_RNA)
test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^mt-")
test[["percent.cd3"]] <- PercentageFeatureSet(test, features = c("Cd3e", "Cd3g", "Cd3d"))
test[["percent.FOXP3"]] <- PercentageFeatureSet(test, features = c("Foxp3"))


FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(test, feature1="nFeature_RNA", feature2="percent.mt")

VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cd3"), ncol = 4)
VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cd3"), ncol = 4, pt.size = 0.1, group.by = "orig.ident")


test <- subset(test, subset = nFeature_RNA < 6000 & percent.mt < 15 & percent.cd3 > 0.055)





########################
#Project sample dataset onto ProjecTILs reference plot
########################
####Projection algorithm
# BiocManager::install("BiocParallel")
query.projected <- make.projection(test, ref=ref)
test <- NULL
Hs2Mm.convert.table <- NULL



#pre-filters T cells using scGate
palette <- refCols
plot.projection(ref, query.projected, cols = refCols)
table(query.projected$orig.ident)




###Predict cell states in query
query.projected <- cellstate.predict(ref=ref, query=query.projected)

Idents(query.projected, cells = WhichCells(query.projected, idents = c("Tfh"))) <- "Other"
celltypes <- c("Other","Treg","Th1","CD4_NaiveLike","CD8_NaiveLike","CD8_EarlyActiv","CD8_EffectorMemory","CD8_Tpex","CD8_Tex")
for (i in celltypes){
  Idents(query.projected, cells = WhichCells(query.projected, idents = i)) <- i
}

query.projected$functional.cluster <- Idents(query.projected)
table(query.projected$functional.cluster)



plot.projection(ref, query.projected, cols = refCols) + labs(title = NULL)
plot.projection(ref, cols = refCols) + labs(title = NULL)





###make a violin plot of markers in our data (not the reference)
markers <- c("Cd4","Cd8a","Ccr7","Tcf7","Izumo1r","Cxcr6","Xcl1","Gzmb","Gzmk","Ifng","Foxp3")
VlnPlot(query.projected, features = markers, stack = T, flip = T, fill.by = "ident", cols=refCols, assay = "RNA")


T_markers <- c("Cd8a","Cd4","Cd44")
Ex_markers <- c("Pdcd1","Havcr2","Tox")
pEx_markers <- c("Ifng","Xcl1")
EM_markers <- c("Cxcr6","Eomes","Gzmb","Gzma","Gzmk")
EA_markers <- c()
NL8_markers <- c("Sell","Dapl1","Tcf7","Lef1","Ccr7")
NL4_markers <- c()
Th1_markers <- c("Tbx21","Cxcr3","Stat1","Stat4")
Treg_markers <- c("Foxp3","Izumo1r","Icos","Ikzf2")
markers <- c(T_markers,Ex_markers,pEx_markers,EM_markers,EA_markers,NL8_markers,NL4_markers,Th1_markers,Treg_markers)


# gd_markers <-c("Maf","Ccr2","Icos")
# EfCD4_markers <- c("Klrg1","Bhlhe40","Cxcr3","Ifng")
# NCD8_markers <- c("Sell","Dapl1","Tcf7","Lef1","Ccr7")
# EMCD8_markers <- c("Cxcr6","Eomes","Gzmb","Gzma","Gzmk")
# NK_markers <- c("Klrd1","Klrc2","Gata3","Il2rb","Runx3")
# Other_markers <- c("Foxp3","Il7r","Id2")
# markers <- c(T_markers, gd_markers, EfCD4_markers, NCD8_markers, EMCD8_markers, NK_markers, Other_markers)


DotPlot(query.projected, features = markers, assay = "RNA",dot.scale = 10) + 
  scale_color_gradient2(low="#4575b4", mid = "#ffffbf",  high="#d73027") +
  RotatedAxis()


#reverse the order of celltypes
rev_celltypes <- rev(celltypes)
for (i in rev_celltypes){
  Idents(query.projected, cells = WhichCells(query.projected, idents = i)) <- i
}



query.projected <- ScaleData(query.projected, features = rownames(query.projected))
DoHeatmap(subset(query.projected, downsample = 100), features = markers, size = 1, assay = "RNA")

head(query.projected@assays$RNA@scale.data)
DoHeatmap(query.projected, features = VariableFeatures(query.projected)[1:100], size = 4,
          angle = 90)






########################
#cell counts per subtype
########################
sub<- table(query.projected$functional.cluster)
sub

#display counts as percentages
plot.statepred.composition(ref, query.projected,metric = "Percent", cols = refCols)
#radar plots; distribution of marker genes for each cell subtype. compare reference to query
plot.states.radar(ref, query=query.projected, min.cells=30)


summary(query.projected$orig.ident)








########################
#subset the total query to ZsG- and ZsG+ subgroups
########################

zshigh <- subset(query.projected, subset= orig.ident == 1)
zslow <- subset(query.projected, subset= orig.ident == 2)


###visualize projection of ZsG- or ZsG+ subgroups onto the reference UMAP
plot.projection(ref, zshigh, cols = refCols)+ labs(title = NULL)
plot.projection(ref, zslow, cols = refCols) + labs(title = NULL)




########################
#Calculate DEGs between ZsG- and ZsG+ cells within the same subtype
########################

### finds degs between groups, make volcano plots, write csv files
celltypes <- c("all", "CD8_Tex", "CD8_Tpex", "CD8_EffectorMemory", "CD8_EarlyActiv", "CD8_NaiveLike", "Th1", "Treg")

dir <- "output/"
for (i in 1:8){
  type <- celltypes[i]

  #want to use DESeq2 because group sizes are relatively small for scRNAseq data, but because every gene contains at least one zero, can't calculate geometric means
  #stick with the wilcox test
  discriminantGenes <- find.discriminant.genes(ref=ref, query.control=zslow,
                                               query=zshigh, state=type, test="wilcox", genes.use="all",
                                               logfc.threshold = 0, min.pct = 0)
  
  filename <- paste(dir, type, "_all_DEGs.csv", sep="")
  write.csv(discriminantGenes, file = filename, row.names=T)
}






##################
#2022/11/10
#Make a violin plot comparing expression of Helios/Ikzf2 in Tregs, ZsLow vs ZsHigh
##################
sub <- subset(query.projected, subset = functional.cluster == "Treg")
sub$group <- ifelse(sub$orig.ident == 1, "ZsHigh", "ZsLow")
Idents(sub) <- sub$group

#Plot ZsLow first, then ZsHigh
#https://www.biostars.org/p/446743/
subplot <- VlnPlot(sub, features = "Ikzf2", flip = T, split.by = "group", assay = "RNA", cols = c("#118002", "gray"), pt.size = 1, add.noise = F)
subplot$data$ident <- factor(x = subplot$data$ident, levels = c("ZsLow", "ZsHigh"))
subplot


#Plot of Tbet/Tbx21
subplot <- VlnPlot(sub, features = "Tbx21", flip = T, split.by = "group", assay = "RNA", cols = c("#118002", "gray"), pt.size = 1, add.noise = F)
subplot$data$ident <- factor(x = subplot$data$ident, levels = c("ZsLow", "ZsHigh"))
subplot






